extern crate clap;
use clap::{App, Arg, ArgMatches, SubCommand};

extern crate serde_json;

use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::*;

extern crate bio_algos;
use bio_algos::{analysis::PWM, *};

/// This program allows using different algoryihms/datastructures offered by
/// the bio_algos crate on the command line.
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
        )
        .subcommand(
            SubCommand::with_name("calc_pwm")
                .about("Calculates a position weight matrix for the input")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                )
                .arg(
                    Arg::with_name("output")
                        .help("Where to store the pwm")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("tis_position")
                        .help("Position of the aligned translation initiation site")
                        .long("tis")
                        .takes_value(true)
                        .value_name("POS")
                        .default_value("100"),
                )
                .arg(
                    Arg::with_name("pseudocount")
                        .help("Pseudocount to use for the pwm")
                        .long("pseudocount")
                        .takes_value(true)
                        .default_value("1.0")
                        .value_name("COUNT"),
                )
                .arg(
                    Arg::with_name("length")
                        .help("Length of the pwm")
                        .long("length")
                        .short("L")
                        .takes_value(true)
                        .default_value("30")
                        .value_name("LENGTH"),
                )
                .arg(
                    Arg::with_name("offset")
                        .help(
                            "Offset of the pwm window. A positive value will move the pwm \
                        to the right, over the tis site, a negative value in the other direction",
                        )
                        .next_line_help(true)
                        .long("offset")
                        .takes_value(true)
                        .default_value("0")
                        .value_name("COUNT"),
                )
                .arg(
                    Arg::with_name("background_model")
                        .help("Specify a background model to use for the pwm \
                        [default: estimate background distribution of wrong sites]")
                        .long("background")
                        .value_names(&["A", "C", "G", "T"])
                ),
        )
        .subcommand(
            SubCommand::with_name("calc_threshold")
                .about("Calculates a score threshold give a sensitivity and an optional pwm matrix")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                )
                .arg(
                    Arg::with_name("pwm")
                        .help("Position weight matrix for given sequence data")
                        .long("pwm")
                        .takes_value(true)
                        .value_name("PATH"),
                ),
        )
        .subcommand(
            SubCommand::with_name("calc_roc")
                .about("Calculates a ROC curve given the inout sequences and the pwm matrix")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                )
                .arg(
                    Arg::with_name("pwm")
                        .help("Position weight matrix for given sequence data")
                        .long("pwm")
                        .required(true)
                        .takes_value(true)
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("output")
                        .help("Where to store the roc data")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("auc")
                        .help("Output the AUC for the curce")
                        .long("auc"),
                ),
        )
        .subcommand(
            SubCommand::with_name("eval_sequences")
                .about("Evaluates the input sequences given the provided pwm and a threshold")
                .arg(
                    Arg::with_name("INPUT")
                        .help("Input sequence data")
                        .required(true),
                )
                .arg(
                    Arg::with_name("pwm")
                        .help("Position weight matrix for given sequence data")
                        .long("pwm")
                        .required(true)
                        .takes_value(true)
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("threshold")
                        .help("Threshold to use for evalutaion")
                        .short("t")
                        .required(true)
                        .takes_value(true),
                ),
        )
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("count_startcodons") {
        exec_count_startcodons(matches);
    } else if let Some(matches) = matches.subcommand_matches("calc_pwm") {
        if let Err(err) = exec_calc_pwm(matches) {
            print_and_exit(err);
        }
    } else if let Some(matches) = matches.subcommand_matches("calc_threshold") {
        exec_calc_threshold(matches);
    } else if let Some(matches) = matches.subcommand_matches("calc_roc") {
        if let Err(err) = exec_calc_roc(matches) {
            print_and_exit(err);
        }
    } else if let Some(matches) = matches.subcommand_matches("eval_sequences") {
        if let Err(err) = exec_eval_sequences(matches) {
            print_and_exit(err);
        }
    }
}

fn exec_count_startcodons(matches: &ArgMatches) {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let codon_count_map = analysis::count_startcodon_variants(&sequences);
    println!("Results:\n{:?}", codon_count_map)
}

fn exec_calc_pwm(matches: &ArgMatches) -> Result<(), Box<Error>> {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let tis_position: usize = matches.value_of("tis_position").unwrap().parse()?;
    let offset: i32 = matches.value_of("offset").unwrap().parse()?;
    let pseudocount: f64 = matches.value_of("pseudocount").unwrap().parse()?;
    let length: usize = matches.value_of("length").unwrap().parse()?;
    let background_model: [f64; 4] = match matches.values_of("background_model") {
        Some(values) => {
            let bg_in: Vec<f64> = values.map(|s| s.parse::<f64>()).collect::<Result<_, _>>()?;
            [bg_in[0], bg_in[1], bg_in[2], bg_in[3]]
        }
        None => analysis::calc_background_model(&sequences, tis_position, length),
    };
    let pwm = analysis::PWM::new(
        &sequences,
        tis_position,
        offset,
        length,
        pseudocount,
        background_model,
    );

    if let Some(path) = matches.value_of("output") {
        match pwm.store(path) {
            Err(msg) => eprintln!("{}", msg),
            _ => println!("Stored pwm at: {}", path),
        };
    } else {
        println!("{:?}", pwm)
    }
    Ok(())
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

fn exec_calc_roc<'a>(matches: &'a ArgMatches) -> Result<(), &'a str> {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let tis_position = 100;
    let pwm_path = matches.value_of("pwm").unwrap();
    let pwm = PWM::load(pwm_path)?;
    let roc = match pwm.calc_roc(&sequences, tis_position) {
        Ok(roc) => roc,
        Err(_) => return Err("Error calculating roc"),
    };
    match matches.value_of("output") {
        Some(path) => {
            let path = Path::new(path);
            let mut file = match File::create(path) {
                Ok(file) => file,
                _ => return Err("Error opening file."),
            };
            let json = match serde_json::to_string(&roc) {
                Ok(json) => json,
                _ => return Err("Error converting roc data to json."),
            };
            if let Err(_) = file.write_all(json.as_bytes()) {
                return Err("Unable to write json data.");
            }
        }
        None => println!("ROC:\n{:?}", roc),
    };
    if matches.is_present("auc") {
        println!("AUC: {}", analysis::calc_auc(&roc));
    }
    Ok(())
}

fn exec_eval_sequences<'a>(matches: &'a ArgMatches) -> Result<(), &'a str> {
    let tis_position = 100;
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let threshold: f64 = match matches.value_of("threshold").unwrap().parse() {
        Ok(val) => val,
        _ => return Err("Error parsing threshold input"),
    };
    let pwm_path = matches.value_of("pwm").unwrap();
    let pwm = PWM::load(pwm_path)?;
    println!(
        "{:?}",
        pwm.eval_sequences(&sequences, threshold, tis_position)
    );
    Ok(())
}

fn print_and_exit<T: Display>(err: T) -> ! {
    eprintln!("Aborting!\n{}", err);
    process::exit(1);
}
