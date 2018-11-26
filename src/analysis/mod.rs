//! This module contains contains methods and a position weight matrix
//! implementation useful for detecting translataion initiation sites.  
extern crate ndarray;
use ndarray::{Array2};

extern crate serde_json;

extern crate serde_derive;
use serde_derive::{Deserialize, Serialize};

use std::collections::HashMap;
use std::error::Error;
use std::f64;
use std::fmt;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::ops::{Add, AddAssign};
use std::path::Path;

use super::{Base, Sequence, SequenceExt, Startcodon};

/// Holds the results of evaluating a pwm with threshold on sequences.
#[derive(Default, Debug)]
pub struct Evaluation {
    false_positive: usize,
    true_positive: usize,
    negative_count: usize,
    positive_count: usize,
}

impl Evaluation {
    /// Returns the true positive rate for the Evaluation.
    pub fn tpr(&self) -> f64 {
        self.true_positive as f64 / self.positive_count as f64
    }

    /// Returns the false positive rate for the Evaluation.
    pub fn fpr(&self) -> f64 {
        self.false_positive as f64 / self.negative_count as f64
    }
}

impl fmt::Display for Evaluation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "FP: {}, TP: {}, FPR: {}, TPR: {}",
            self.false_positive,
            self.true_positive,
            self.fpr(),
            self.tpr()
        )
    }
}

/// Implement addition on Evalution.
impl Add for Evaluation {
    type Output = Evaluation;

    fn add(mut self, other: Evaluation) -> Evaluation {
        self.false_positive += other.false_positive;
        self.true_positive += other.true_positive;
        self.negative_count += other.negative_count;
        self.positive_count += other.negative_count;
        self
    }
}

/// Implement addition and assignment ( += ) on Evalutaion.
impl AddAssign for Evaluation {
    fn add_assign(&mut self, other: Evaluation) {
        self.false_positive += other.false_positive;
        self.true_positive += other.true_positive;
        self.negative_count += other.negative_count;
        self.positive_count += other.negative_count;
    }
}

/// A PWM (position weight matrix) ca be used to identify possible
/// translattion initiation (tis) sites. This is achieved by training
/// the matrix on a set of training sequences which have been aligned
/// such that their actual tis site is always at the same position in
/// the sequence.  
/// As hypterparamters, the length of the slice, the offset from the
/// tis site, a pseudocount and a background distribution can be
/// specified. The background model can also be estimated
/// with the calc_auc() function of this module.
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
        // Every column of the matrix will contain number of sequences + 4 times
        // the pseudocount elements. Thus we can divide the matrix elementwise
        // by that number, in order to get a position probability
        matrix /= sequences.len() as f64 + 4. * pseudocount;

        for (row_index, row) in matrix.genrows_mut().into_iter().enumerate() {
            for el in row {
                *el = el.log2() - background_distribution[row_index].log2()
            }
        }

        PWM {
            matrix,
            offset: pwm_offset,
            pseudocount,
            background_distribution,
        }
    }

    /// Returns the length of the PWM.
    pub fn len(&self) -> usize {
        self.matrix.shape()[1]
    }

    /// Stores the PWM as json at the specified path.
    pub fn store(&self, path: &str) -> io::Result<()> {
        let path = Path::new(path);
        let mut file = File::create(path)?;
        file.write_all(serde_json::to_string(self)?.as_bytes())?;
        Ok(())
    }

    /// Loads a json serialized PWM from the specified path. 
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

    /// Calculates a score for the specied sequence slice
    pub fn score(&self, sequence_slice: &[Base]) -> Result<f64, &str> {
        if sequence_slice.len() != self.len() {
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

    /// Returns a Vector of tuples where the first element designates
    /// the true label of a potential TIS and the second element the 
    /// score of this site. 
    pub fn score_label_sequence(
        &self,
        sequence: &Sequence,
        tis_position: usize,
    ) -> Vec<(bool, f64)> {
        let mut scores = Vec::new();

        for i in self.len()..sequence.len() {
            if let None = sequence.contains_startcodon_at(i) {
                if i == tis_position {
                    panic!(
                        "Invalid input data! Sequence does not contain startcodon at {}\n {:?}",
                        tis_position, sequence
                    )
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

    /// Returns a Vector of tuples where the first element designates
    /// the true label of a potential TIS and the second element the 
    /// score of this site. 
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

    /// Evaluates a single sequence. The returned Evaluation contains
    /// the true positive count, the false positive count, as well as
    /// the total number of actual and false TIS candidates.  
    /// It also offers methods to calculate the TPR and FPR.
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
                    eval.positive_count += 1;
                    if score >= threshold {
                        eval.true_positive += 1;
                    }
                } else {
                    eval.negative_count += 1;
                    if score >= threshold {
                        eval.false_positive += 1;
                    }
                }
                eval
            })
    }

    /// Evaluates a multiple sequences. The returned Evaluation contains
    /// the true positive count, the false positive count, as well as
    /// the total number of actual and false TIS candidates.  
    /// It also offers methods to calculate the TPR and FPR.
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

    /// This is method is used to calculate the points for a receiver operator curve (ROC).
    /// The Vector contained in the Result holds tuples where the first element is the
    /// x and the second element the z value of the respecive ROC point.
    /// 
    /// Algorithm adapted from http://www.cs.ru.nl/~tomh/onderwijs/dm/dm_files/ROC101.pdf .
    pub fn calc_roc(
        &self,
        sequences: &Vec<Sequence>,
        tis_position: usize,
    ) -> Result<Vec<(f64, f64)>, &str> {
        // Calculate a vector that contains a (true label, score) tuple for every
        // possible TIS in the given sequences.
        let mut label_score_pairs = self.score_label_sequences(sequences, tis_position);
        // Each sequence contains exactly one actual TIS
        let positive_count = sequences.len();
        let negative_count = label_score_pairs.len() as i32 - positive_count as i32;
        if positive_count == 0 || negative_count <= 0 {
            return Err("Positive and negative count must be greater than 0.");
        }
        // sort the label, score pairs by score descending
        label_score_pairs
            .sort_by(|(_, score_fst), (_, score_snd)| score_snd.partial_cmp(score_fst).unwrap());
        let mut false_positive = 0;
        let mut true_positive = 0;
        let mut score_prev = f64::NEG_INFINITY;
        let mut roc_points = Vec::new();
        for (label, score) in label_score_pairs {
            // only push a new curve point if the score is different than the last one
            // this prevents unjustified curves and will result in a diagonal
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

/// Calculates the are under curve (AUC) for a Vector of ROC points.  
/// The tuples are (x_val, y_val).
pub fn calc_auc(roc_points: &Vec<(f64, f64)>) -> f64 {
    roc_points
        .into_iter()
        .zip(roc_points[1..].into_iter())
        .fold(0., |acc, ((x1, y1), (x2, y2))| {
            acc + trapezoid_area(*y1, *y2, (x1 - x2).abs())
        })
}

/// Calculates the AUC for a range of pseudocounts and outpouts a vector of
/// (pseudocount, AUC) values.
pub fn eval_auc_for_pseudocount(
    sequences_train: &Vec<Sequence>,
    sequences_test: &Vec<Sequence>,
    pseudocount_start: f64,
    pseudocount_end: f64,
    pseudocount_step: f64,
    tis_position: usize,
    pwm_offset: i32,
    length: usize,
    background_distribution: [f64; 4],
) -> Result<Vec<(f64, f64)>, Box<Error>> {
    let mut pseudocount = pseudocount_start;
    let mut count_auc_points: Vec<(f64, f64)> = Vec::new();
    while pseudocount <= pseudocount_end {
        let pwm = PWM::new(
            sequences_train,
            tis_position,
            pwm_offset,
            length,
            pseudocount,
            background_distribution,
        );
        let roc = pwm.calc_roc(sequences_test, tis_position)?;
        count_auc_points.push((pseudocount, calc_auc(&roc)));
        pseudocount += pseudocount_step;
    }
    Ok(count_auc_points)
}

/// Estimates the background model of the different bases for false
/// TIS candidates in the given sequences.
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

/// Given a PWM, sequences and a tis_position calculate the threshold needed
/// to achieve the specified sensitivity.  
/// Returns the needed threshold and an Evaluation object.
pub fn calc_score_threshold_for_sensitivity(
    pwm: PWM,
    sequences: &Vec<Sequence>,
    tis_position: usize,
    sensitivity: f64,
) -> (f64, Evaluation) {
    if sensitivity <= 0. || sensitivity > 1. {
        panic!("Sensitivity must be within (0, 1] .")
    }

    let mut label_score_pairs = pwm.score_label_sequences(sequences, tis_position);
    label_score_pairs
        .sort_by(|(_, score_fst), (_, score_snd)| score_snd.partial_cmp(score_fst).unwrap());
    let mut eval = Evaluation::default();
    let mut threshold = 0.;
    eval.positive_count = sequences.len();
    eval.negative_count = label_score_pairs.len() - eval.positive_count;

    for (label, score) in label_score_pairs {
        threshold = score;
        if label {
            eval.true_positive += 1;
        } else {
            eval.false_positive += 1;
        }
        if eval.tpr() >= sensitivity {
            break;
        }
    }

    (threshold, eval)
}

/// Count the number of the different startcodon variants in the given sequences.
pub fn count_startcodon_variants(sequences: &Vec<Sequence>) -> HashMap<Startcodon, u32> {
    sequences.iter().fold(HashMap::new(), |map, seq| {
        seq.count_startcodon_variants(map)
    })
}

/// Area of a trapezoid.
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
