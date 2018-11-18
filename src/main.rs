extern crate clap;
use clap::{App, Arg, ArgMatches, SubCommand};

use std::*;

extern crate bio_algos;
use bio_algos::{analysis::PWM, *};

fn main() {
    let matches = App::new("Bio-Algorithms")
        .version("0.1.0")
        .author("Robin William Hundt")
        .about("A program implementing some bio algorithmic stuff.")
        .subcommand(
            SubCommand::with_name("count_startcodons")
                .about("Counts all possible startcodon variants in the input")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                ),
        ).subcommand(
            SubCommand::with_name("calc_pwm")
                .about("Calculates a position weight matrix for the input")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                ).arg(
                    Arg::with_name("output")
                        .help("Where to store the pwm")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .value_name("PATH"),
                ),
        ).subcommand(
            SubCommand::with_name("calc_threshold")
                .about("Calculates a score threshold give a sensitivity and an optional pwm matrix")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                ).arg(
                    Arg::with_name("pwm")
                        .help("Position weight matrix for given sequence data")
                        .long("pwm")
                        .takes_value(true)
                        .value_name("PATH"),
                ),
        ).subcommand(
            SubCommand::with_name("calc_roc")
                .about("Calculates a ROC curve given the inout sequences and the pwm matrix")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                ).arg(
                    Arg::with_name("pwm")
                        .help("Position weight matrix for given sequence data")
                        .long("pwm")
                        .required(true)
                        .takes_value(true)
                        .value_name("PATH"),
                ),
        ).get_matches();

    if let Some(matches) = matches.subcommand_matches("count_startcodons") {
        exec_count_startcodons(matches);
    } else if let Some(matches) = matches.subcommand_matches("calc_pwm") {
        exec_calc_pwm(matches);
    } else if let Some(matches) = matches.subcommand_matches("calc_threshold") {
        exec_calc_threshold(matches);
    }
}

fn exec_count_startcodons(matches: &ArgMatches) {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let codon_count_map = analysis::count_startcodon_variants(&sequences);
    println!("Results:\n{:?}", codon_count_map)
}

fn exec_calc_pwm(matches: &ArgMatches) {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let pwm = analysis::PWM::new(&sequences, 100, 0, 30, 1., [0.25, 0.25, 0.25, 0.25]);
    if let Some(path) = matches.value_of("output") {
        match pwm.store(path) {
            Err(msg) => eprintln!("{}", msg),
            _ => println!("Stored pwm at: {}", path),
        };
    } else {
        println!("{:?}", pwm)
    }
}

fn exec_calc_threshold(matches: &ArgMatches) {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let pwm = match matches.value_of("pwm") {
        Some(path) => match PWM::load(path) {
            Err(msg) => {
                eprintln!("{}", msg);
                process::exit(1);
            }
            Ok(pwm) => pwm,
        },
        None => PWM::new(&sequences, 100, 0, 30, 1., [0.25, 0.25, 0.25, 0.25]),
    };
    let tis_position = 100;
    let sensitivity = 0.5;
    let max_iterations = 100;
    let (threshold, achieved_sensitivity, eval) = analysis::calc_score_threshold_for_sensitivity(
        pwm,
        &sequences,
        tis_position,
        sensitivity,
        max_iterations,
    );
    println!(
        "Acvieved sensitivity of {} with threshold: {}",
        achieved_sensitivity, threshold
    );
    println!("{:?}", eval);
}
