//! This crate contains algorithms and datastructures for the fulfillment of the
//! the first assignment of the w*ourse at the
//! university of Göttingen.  
//! The crate is divided this base module, an analysis module and a util module.  
//! In the base module some types and methods on them have been implemented (Base, Startcodon, Sequence).
//! The Analysis module contains an implementation of a position weight matrix, as well as
//! methods that can be used to evaluate a pwm on aligned sequences.  
//! The util module contains some utility functions for reading sequence data from a file.
//!
use std::collections::HashMap;

use self::Base::*;

pub mod analysis;
pub mod util;

/// Enumeration representing the individual bases of a genome sequence.
#[derive(Debug, Clone, Copy)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

/// Enumeration representing the different possible Startcodons which are at the start of
/// translation initiation sites.
#[derive(PartialEq, Eq, Hash, Debug)]
pub enum Startcodon {
    ATG,
    GTG,
    TTG,
}

/// A type alias for Sequence to a Vector of bases. Note that a Vector is a
/// heap allocated dynamically resizing array.
pub type Sequence = Vec<Base>;

/// Since Sequence is only a type alias, methods can not be directly implemented in it.
/// To circumvent this, a trait SequenceExt is defined which contains useful methods on
/// Sequences like parse().
pub trait SequenceExt {
    /// Parses the input string into a Sequence. Returns an error in case the string
    /// contains characters other than A, C, G or T.
    fn parse(input: &str) -> Result<Sequence, &'static str> {
        input
            .chars()
            .map(|c| match c {
                'A' => Ok(A),
                'C' => Ok(C),
                'G' => Ok(G),
                'T' => Ok(T),
                _ => return Err("Illegal character: "),
            })
            .collect()
    }
    /// Returns whether the sequence contains a Startcodon at the specified position.
    fn contains_startcodon_at(&self, index: usize) -> Option<Startcodon>;
    // Counts all the possible startcodon variants in the sequence.
    fn count_startcodon_variants(&self, map: HashMap<Startcodon, u32>) -> HashMap<Startcodon, u32>;
}

impl SequenceExt for Sequence {
    fn contains_startcodon_at(&self, index: usize) -> Option<Startcodon> {
        if index + 3 > self.len() {
            return None;
        }
        match self[index..index + 3] {
            [A, T, G] => Some(Startcodon::ATG),
            [G, T, G] => Some(Startcodon::GTG),
            [T, T, G] => Some(Startcodon::TTG),
            _ => None,
        }
    }

    fn count_startcodon_variants(&self, map: HashMap<Startcodon, u32>) -> HashMap<Startcodon, u32> {
        // iterate and enumerate all bases in sequence, then reduce to a HashMap which
        // counts the occurences of the startcodons
        self.iter().enumerate().fold(map, |mut map, (i, _)| {
            if let Some(codon) = self.contains_startcodon_at(i) {
                let count = map.entry(codon).or_insert(0);
                *count += 1;
            }
            map
        })
    }
}
