use std::collections::HashMap;

use self::Base::*;

pub mod analysis;
pub mod util;

#[derive(Debug, Clone, Copy)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

#[derive(PartialEq, Eq, Hash, Debug)]
pub enum Startcodon {
    ATG,
    GTG,
    TTG,
}

pub type Sequence = Vec<Base>;

pub trait SequenceExt {
    fn parse(input: &str) -> Result<Sequence, &'static str>;
    fn contains_startcodon_at(&self, index: usize) -> Option<Startcodon>;
    fn count_startcodon_variants(&self, map: HashMap<Startcodon, u32>) -> HashMap<Startcodon, u32>;
}

impl SequenceExt for Sequence {
    fn parse(input: &str) -> Result<Sequence, &'static str> {
        input
            .chars()
            .map(|c| match c {
                'A' => Ok(A),
                'C' => Ok(C),
                'G' => Ok(G),
                'T' => Ok(T),
                _ => return Err("Illegal character: "),
            }).collect()
    }

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
        self.iter().enumerate().fold(map, |mut map, (i, _)| {
            if let Some(codon) = self.contains_startcodon_at(i) {
                let count = map.entry(codon).or_insert(0);
                *count += 1;
            }
            map
        })
    }
}
