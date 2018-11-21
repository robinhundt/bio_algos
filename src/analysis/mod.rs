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

    pub fn store(&self, path: &str) -> io::Result<()> {
        let path = Path::new(path);
        let mut file = File::create(path)?;
        file.write_all(serde_json::to_string(self)?.as_bytes())?;
        Ok(())
    }

    pub fn load(path: &str) -> Result<Self, &str> {
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
                if i == tis_position {
                    panic!("SHOULDNT HAPPEN! {:?}", sequence)
                }
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
            .sort_by(|(_, score_fst), (_, score_snd)| score_snd.partial_cmp(score_fst).unwrap());
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
        }
        roc_points.push((
            false_positive as f64 / negative_count as f64,
            true_positive as f64 / positive_count as f64,
        ));
        Ok(roc_points)
    }
}

pub fn calc_auc(roc_points: &Vec<(f64, f64)>) -> f64 {
    roc_points
        .into_iter()
        .zip(roc_points[1..].into_iter())
        .fold(0., |acc, ((x1, y1), (x2, y2))| {
            acc + trapezoid_area(*y1, *y2, (x1 - x2).abs())
        })
}

pub fn calc_background_model(
    sequences: &Vec<Sequence>,
    tis_position: usize,
    window_size: usize,
) -> [f64; 4] {
    let mut bg_model = [0.; 4];
    for seq in sequences {
        for i in window_size..seq.len() {
            if seq.contains_startcodon_at(i).is_none() || i == tis_position {
                continue;
            }
            for base in seq[i - window_size..i].into_iter() {
                bg_model[*base as usize] += 1.;
            }
        }
    }
    let sum_occurences: f64 = bg_model.iter().sum();
    bg_model
        .iter_mut()
        .for_each(|count| *count /= sum_occurences);
    bg_model
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

    let mut label_score_pairs = pwm.score_label_sequences(sequences, tis_position);
    label_score_pairs
        .sort_by(|(_, score_fst), (_, score_snd)| score_snd.partial_cmp(score_fst).unwrap());
    let mut eval = Evaluation::default();
    let mut threshold = 0.;
    let mut achieved_sensitivity = 0.;
    for (label, score) in label_score_pairs {
        threshold = score;
        if label {
            eval.true_positive += 1;
        } else {
            eval.false_positive += 1;
        }
        achieved_sensitivity = eval.true_positive as f64 / sequences.len() as f64;
        if achieved_sensitivity >= sensitivity {
            break;
        }
    }
    (threshold, achieved_sensitivity, eval)
}

pub fn count_startcodon_variants(sequences: &Vec<Sequence>) -> HashMap<Startcodon, u32> {
    sequences.iter().fold(HashMap::new(), |map, seq| {
        seq.count_startcodon_variants(map)
    })
}

fn trapezoid_area(base1: f64, base2: f64, height: f64) -> f64 {
    height * (base1 + base2) / 2.
}

#[cfg(test)]
mod tests {
    use self::Base::*;
    use super::*;
    #[test]
    fn auc_calculates_zero() {
        let roc = vec![(0., 0.), (1., 0.), (1., 1.)];
        let auc = calc_auc(&roc);
        assert_eq!(auc, 0.);
    }
    #[test]
    fn auc_calculates_one() {
        let roc = vec![(0., 0.), (0., 1.), (1., 1.)];
        let auc = calc_auc(&roc);
        assert_eq!(auc, 1.);
    }
    #[test]
    fn auc_calculates_0_5() {
        let roc = vec![(0., 0.), (1., 1.)];
        let auc = calc_auc(&roc);
        assert_eq!(auc, 0.5);
    }
    #[test]
    fn auc_calculates_0_625() {
        let roc = vec![
            (0., 0.),
            (0., 0.25),
            (0.25, 0.25),
            (0.5, 0.5),
            (0.5, 0.75),
            (0.75, 1.),
            (1., 1.),
        ];
        let auc = calc_auc(&roc);
        assert_eq!(auc, 0.625);
    }

    #[test]
    fn calc_background_model_uniform() {
        let seqs = vec![
            vec![A, C, G, T, A, T, G, A, A, G, T, G],
            vec![A, C, G, T, A, T, G, A, A, G, T, G],
        ];
        let bg_model = calc_background_model(&seqs, 9, 4);
        assert_eq!(bg_model, [0.25; 4]);
    }
}
