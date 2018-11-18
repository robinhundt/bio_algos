use std::{fs, io, io::BufRead, process};

use super::{Sequence, SequenceExt};

pub fn read_sequences_from_file(path: &str) -> io::Result<Vec<Result<Sequence, &'static str>>> {
    let f = fs::File::open(path)?;
    let f = io::BufReader::new(f);
    Ok(f.lines()
        .map(|line_res| match line_res {
            Ok(line) => Sequence::parse(line.as_str()),
            Err(_) => Err("Encountered io error when reading sequence file"),
        }).collect())
}

pub fn read_sequences_or_exit(path: &str) -> Vec<Sequence> {
    let input = read_sequences_from_file(path);

    let sequences = input.unwrap_or_else(|e| {
        eprintln!("Encountered error while reading sequence data: {:?}", e);
        process::exit(1);
    });
    sequences
        .into_iter()
        .filter(|el| {
            if el.is_err() {
                eprintln!("Error in sequence {:?}", el);
                return false;
            }
            true
        }).map(|el| el.unwrap())
        .collect()
}
