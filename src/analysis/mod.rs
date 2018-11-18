extern crate ndarray;
use ndarray::{Array2, Zip};

extern crate serde_json;

extern crate serde_derive;
use serde_derive::{Deserialize, Serialize};

use std::collections::HashMap;
use std::f64;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::ops::{Add, AddAssign};
use std::path::Path;

use super::{Base, Sequence, SequenceExt, Startcodon};

#[derive(Default, Debug)]
pub struct Evaluation {
    false_positive: usize,
    true_positive: usize,
}

impl Add for Evaluation {
    type Output = Evaluation;

    fn add(mut self, other: Evaluation) -> Evaluation {
        self.false_positive += other.false_positive;
        self.true_positive += other.true_positive;
        self
    }
}

impl AddAssign for Evaluation {
    fn add_assign(&mut self, other: Evaluation) {
        self.false_positive += other.false_positive;
        self.true_positive += other.true_positive;
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PWM {
    matrix: Array2<f64>,
    offset: i32,
    pseudocount: f64,
    background_distribution: [f64; 4],
}

impl PWM {
    /// This function calculates a position weight matrix for the given sequences.
    /// The length of the matrix can be adjusted, as well as the offset for the last position
    /// of the matrix with respect to the TIS site. E.g. when the position containing the first base
    /// of the Startcodon should be included in the PWM, you can pass an offset of 1.
    pub fn new(
        sequences: &Vec<Sequence>,
        tis_position: usize,
        pwm_offset: i32,
        length: usize,
        pseudocount: f64,
        background_distribution: [f64; 4],
    ) -> PWM {
        // the pwm has a dim of (4, length) because there are 4 Bases
        let mut matrix = Array2::<f64>::zeros((4, length));
        // add the pseudocount to every element in the matrix
        matrix += pseudocount;
        let start = (tis_position - length) as i32 + pwm_offset;
        let end = tis_position as i32 + pwm_offset;
        if start < 0 {
            panic!(
                "Combination of tis_position,\
                 pwm_offset and length results in starting index less than 0."
            )
        }

        for seq in sequences {
            for (idx, base) in seq[start as usize..end as usize].into_iter().enumerate() {
                matrix[[*base as usize, idx]] += 1.;
            }
        }
        matrix /= sequences.len() as f64 + 4. * pseudocount;

        for (row_index, row) in matrix.genrows_mut().into_iter().enumerate() {
            // maybe use map_inplace here
            Zip::from(row).apply(|el| *el = el.log2() - background_distribution[row_index].log2())
        }

        PWM {
            matrix,
            offset: pwm_offset,
            pseudocount,
            background_distribution,
        }
    }

    pub fn len(&self) -> usize {
        self.matrix.shape()[1]
    }

    pub fn store(&self, path: &str) -> io::Result<bool> {
        let path = Path::new(path);
        let mut file = File::create(path)?;
        file.write_all(serde_json::to_string(self)?.as_bytes())?;
        Ok(true)
    }

    pub fn load<'a>(path: &str) -> Result<Self, &str> {
        let file = match File::open(path) {
            Err(_) => return Err("Error reading pwm file"),
            Ok(data) => data,
        };

        match serde_json::from_reader(file) {
            Err(_) => Err("Error parsing pwm file"),
            Ok(pwm) => Ok(pwm),
        }
    }

    pub fn score(&self, sequence_slice: &[Base]) -> Result<f64, &str> {
        if sequence_slice.len() != self.len() {
            eprintln!(
                "pwm len: {}, slice len: {}",
                self.len(),
                sequence_slice.len()
            );
            return Err("sequence_slice must be same length as pwm matrix");
        }

        let score = sequence_slice
            .into_iter()
            .zip(self.matrix.gencolumns())
            .fold(0_f64, |score, (base, pwm_column)| {
                score + pwm_column[*base as usize]
            });
        Ok(score)
    }

    pub fn score_label_sequence(
        &self,
        sequence: &Sequence,
        tis_position: usize,
    ) -> Vec<(bool, f64)> {
        let mut scores = Vec::new();

        for i in self.len()..sequence.len() {
            if let None = sequence.contains_startcodon_at(i) {
                continue;
            }
            let start = (i - self.len()) as i32 + self.offset;
            let end = i as i32 + self.offset;
            if start < 0 || end >= sequence.len() as i32 {
                continue;
            }
            let label = tis_position == i;
            // the result should never contain the error case
            let score = self.score(&sequence[start as usize..end as usize]).unwrap();
            scores.push((label, score));
        }
        scores
    }

    pub fn score_label_sequences(
        &self,
        sequences: &Vec<Sequence>,
        tis_position: usize,
    ) -> Vec<(bool, f64)> {
        sequences.iter().fold(Vec::new(), |mut v, seq| {
            v.append(&mut self.score_label_sequence(seq, tis_position));
            v
        })
    }

    pub fn eval_sequence(
        &self,
        sequence: &Sequence,
        threshold: f64,
        tis_position: usize,
    ) -> Evaluation {
        self.score_label_sequence(sequence, tis_position)
            .into_iter()
            .fold(Evaluation::default(), |mut eval, (label, score)| {
                if label {
                    if score >= threshold {
                        eval.true_positive += 1;
                    }
                } else {
                    if score >= threshold {
                        eval.false_positive += 1;
                    }
                }
                eval
            })
    }

    pub fn eval_sequences(
        &self,
        sequences: &Vec<Sequence>,
        threshold: f64,
        tis_position: usize,
    ) -> Evaluation {
        sequences
            .iter()
            .fold(Evaluation::default(), |mut eval, seq| {
                eval += self.eval_sequence(seq, threshold, tis_position);
                eval
            })
    }


    /// Algorithm adapted from http://www.cs.ru.nl/~tomh/onderwijs/dm/dm_files/ROC101.pdf
    pub fn calc_roc(
        &self,
        sequences: &Vec<Sequence>,
        tis_position: usize,
    ) -> Result<Vec<(f64, f64)>, &str> {
        let mut label_score_pairs = self.score_label_sequences(sequences, tis_position);
        let positive_count = sequences.len();
        let negative_count = label_score_pairs.len() as i32 - positive_count as i32;
        if positive_count == 0 || negative_count <= 0 {
            return Err("Positive and negative count must be greater than 0.");
        }
        label_score_pairs
            .sort_by(|(_, score_fst), (_, score_snd)| score_fst.partial_cmp(score_snd).unwrap());
        let mut false_positive = 0;
        let mut true_positive = 0;
        let mut score_prev = f64::NEG_INFINITY;
        let mut roc_points = Vec::new();
        for (label, score) in label_score_pairs {
            if score != score_prev {
                roc_points.push((
                    false_positive as f64 / negative_count as f64,
                    true_positive as f64 / positive_count as f64,
                ));
                score_prev = score;
            }
            if label {
                true_positive += 1;
            } else {
                false_positive += 1;
            }
            roc_points.push((
                false_positive as f64 / negative_count as f64,
                true_positive as f64 / positive_count as f64,
            ));
        }
        Ok(roc_points)
    }
}

pub fn calc_score_threshold_for_sensitivity(
    pwm: PWM,
    sequences: &Vec<Sequence>,
    tis_position: usize,
    sensitivity: f64,
    max_iterations: usize,
) -> (f64, f64, Evaluation) {
    if sensitivity <= 0. || sensitivity > 1. {
        panic!("Sensitivity must be within (0, 1] .")
    }
    if max_iterations == 0 {
        panic!("max_iterations must be bigger than 0")
    }

    let mut low_score_bound = 0.;
    let mut high_score_bound = 0.;

    for column in pwm.matrix.gencolumns() {
        if let Some((min, max)) = min_max(column.iter()) {
            low_score_bound += min;
            high_score_bound += max;
        }
    }

    if low_score_bound == high_score_bound {
        panic!("Error calculating score bounds, low == high. Possibly the pwm is empty.")
    }

    let mut threshold = (high_score_bound - low_score_bound) / 2.;
    let mut curr_sensitivity = 0.;
    let mut eval = Evaluation::default();
    for i in 0..max_iterations {
        println!("step {} of {}", i, max_iterations);
        // TODO instead of always recalculating the scores,
        // they should be cached and compared against the new threshold
        eval = pwm.eval_sequences(sequences, threshold, tis_position);
        curr_sensitivity = eval.true_positive as f64 / sequences.len() as f64;
        if curr_sensitivity == sensitivity {
            break;
        } else if curr_sensitivity > sensitivity {
            low_score_bound = threshold;
            threshold += (high_score_bound - low_score_bound) / 2.;
        } else {
            high_score_bound = threshold;
            threshold -= (high_score_bound - low_score_bound) / 2.;
        }
    }
    (threshold, curr_sensitivity, eval)
}

pub fn count_startcodon_variants(sequences: &Vec<Sequence>) -> HashMap<Startcodon, u32> {
    sequences.iter().fold(HashMap::new(), |map, seq| {
        seq.count_startcodon_variants(map)
    })
}

fn min_max<T>(mut it: T) -> Option<(T::Item, T::Item)>
where
    T: Iterator,
    T::Item: PartialOrd + Copy,
{
    let head = it.next()?;
    let (mut min, mut max) = (head, head);
    for el in it {
        if el < min {
            min = el;
        } else if el > max {
            max = el;
        }
    }
    Some((min, max))
}
