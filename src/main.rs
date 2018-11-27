extern crate bio_algos;
extern crate clap;
extern crate serde_json;

use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::*;

use clap::{App, Arg, ArgMatches, SubCommand};

use bio_algos::{analysis::PWM, *};

/// This program allows using different algoryihms/datastructures offered by
/// the bio_algos crate on the command line.
fn main() {
    let input_arg = Arg::with_name("INPUT")
        .help("Input sequence data")
        .required(true);
    let output_arg = Arg::with_name("output")
        .help("Where to store the output")
        .short("o")
        .long("output")
        .takes_value(true)
        .value_name("PATH");

    let tis_position_arg = Arg::with_name("tis_position")
        .help("Position of the aligned TIS")
        .long("tis")
        .takes_value(true)
        .value_name("POS")
        .default_value("100");

    let pseudocount_arg = Arg::with_name("pseudocount")
        .help("Pseudocount to use for the PWM")
        .long("pseudocount")
        .takes_value(true)
        .default_value("1.0")
        .value_name("COUNT");

    let length_arg = Arg::with_name("length")
        .help("Length of the PWM")
        .long("length")
        .short("L")
        .takes_value(true)
        .default_value("30")
        .value_name("LENGTH");

    let offset_arg = Arg::with_name("offset")
        .help(
            "Offset of the PWM window. A positive value will move the PWM \
             to the right, over the tis site, a negative value in the other direction",
        )
        .next_line_help(true)
        .long("offset")
        .takes_value(true)
        .default_value("0")
        .value_name("COUNT");

    let background_model_arg = Arg::with_name("background_model")
        .help(
            "Specify a background model to use for the PWM \
             [default: estimate background distribution of wrong sites]",
        )
        .long("background")
        .value_names(&["A", "C", "G", "T"]);

    let pwm_arg = Arg::with_name("pwm")
        .help("Position weight matrix for given sequence data")
        .long("pwm")
        .takes_value(true)
        .value_name("PATH");

    let matches = App::new("Bio-Algorithms")
        .version("0.1.0")
        .author("Robin William Hundt")
        .about(
            "A program implementing algorithms and datastructures mostly related to finding TIS",
        )
        .subcommand(
            SubCommand::with_name("count_startcodons")
                .about("Counts all possible startcodon variants in the input")
                .arg(input_arg.clone()),
        )
        .subcommand(
            SubCommand::with_name("calc_pwm")
                .about("Calculates a position weight matrix for the input")
                .arg(input_arg.clone())
                .arg(output_arg.clone())
                .arg(tis_position_arg.clone())
                .arg(pseudocount_arg.clone())
                .arg(length_arg.clone())
                .arg(offset_arg.clone())
                .arg(background_model_arg.clone()),
        )
        .subcommand(
            SubCommand::with_name("calc_threshold")
                .about("Calculates a score threshold given a sensitivity and an optional PWM")
                .arg(input_arg.clone())
                .arg(pwm_arg.clone())
                .arg(tis_position_arg.clone())
                .arg(
                    Arg::with_name("sensitivity")
                        .default_value("0.50")
                        .long("sensitivity")
                        .help("Sensitivity to achieve with threshold"),
                ),
        )
        .subcommand(
            SubCommand::with_name("calc_roc")
                .about("Calculates a ROC curve given the inoput sequences and the PWM")
                .arg(input_arg.clone())
                .arg(pwm_arg.clone())
                .arg(output_arg.clone())
                .arg(
                    Arg::with_name("auc")
                        .help("Output the AUC for the curve")
                        .long("auc"),
                ),
        )
        .subcommand(
            SubCommand::with_name("eval_sequences")
                .arg(input_arg.clone())
                .arg(pwm_arg.clone())
                .about("Evaluates the input sequences given the provided PWM and a threshold")
                .arg(
                    Arg::with_name("threshold")
                        .help("Threshold to use for evalutaion")
                        .short("t")
                        .required(true)
                        .takes_value(true),
                ),
        )
        .subcommand(
            SubCommand::with_name("pseudocount_auc_data")
                .about(
                    "This command will evaluate a range of pseudocounts \
                     and output the AUC of resultig ROC curves",
                )
                .arg(output_arg.clone())
                .arg(tis_position_arg.clone())
                .arg(length_arg.clone())
                .arg(offset_arg.clone())
                .arg(background_model_arg.clone())
                .arg(
                    Arg::with_name("INPUT_TRAIN")
                        .help("Training input sequence data")
                        .required(true),
                )
                .arg(Arg::with_name("INPUT_TEST").help("Testing input sequence data"))
                .arg(
                    Arg::with_name("pseudocount_start")
                        .help("Pseudocount start value to use in the evaluation")
                        .long("count_start")
                        .takes_value(true)
                        .default_value("1.0")
                        .value_name("START"),
                )
                .arg(
                    Arg::with_name("pseudocount_end")
                        .help("Pseudocount end value to use in the evaluation")
                        .long("count_end")
                        .takes_value(true)
                        .default_value("500")
                        .value_name("END"),
                )
                .arg(
                    Arg::with_name("pseudocount_step")
                        .help("Pseudocount step size to use")
                        .long("count_step")
                        .takes_value(true)
                        .default_value("1.0")
                        .value_name("STEP"),
                ),
        )
        .get_matches();

    // match each subcommand to it's respective function
    let execution_result = match matches.subcommand() {
        ("count_startcodons", Some(matches)) => exec_count_startcodons(matches),
        ("calc_pwm", Some(matches)) => exec_calc_pwm(matches),
        ("calc_threshold", Some(matches)) => exec_calc_threshold(matches),
        ("calc_roc", Some(matches)) => exec_calc_roc(matches),
        ("eval_sequences", Some(matches)) => exec_eval_sequences(matches),
        ("pseudocount_auc_data", Some(matches)) => exec_count_auc_pairs(matches),
        (_, _) => Ok(()),
    };

    if let Err(err) = execution_result {
        print_and_exit(err);
    }
}

// The following functions are merely using functionality provided
// by the util and analysis modules.

fn exec_count_startcodons(matches: &ArgMatches) -> Result<(), Box<Error>> {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let codon_count_map = analysis::count_startcodon_variants(&sequences);
    println!("Results:\n{:?}", codon_count_map);
    Ok(())
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
            _ => println!("Stored PWM at: {}", path),
        };
    } else {
        println!("{:?}", pwm)
    }
    Ok(())
}

fn exec_calc_threshold(matches: &ArgMatches) -> Result<(), Box<Error>> {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let pwm = match matches.value_of("PWM") {
        Some(path) => PWM::load(path)?,
        None => PWM::new(&sequences, 100, 0, 30, 1., [0.25, 0.25, 0.25, 0.25]),
    };
    let tis_position: usize = matches.value_of("tis_position").unwrap().parse()?;
    let sensitivity: f64 = matches.value_of("sensitivity").unwrap().parse()?;
    let (threshold, eval) =
        analysis::calc_score_threshold_for_sensitivity(pwm, &sequences, tis_position, sensitivity);
    println!(
        "Acvieved sensitivity of {} with threshold: {}",
        eval.tpr(),
        threshold
    );
    println!("{}", eval);
    Ok(())
}

fn exec_calc_roc<'a>(matches: &'a ArgMatches) -> Result<(), Box<Error>> {
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let tis_position = 100;
    let pwm_path = matches.value_of("pwm").unwrap();
    let pwm = PWM::load(pwm_path)?;
    let roc = pwm.calc_roc(&sequences, tis_position)?;
    match matches.value_of("output") {
        Some(path) => {
            let path = Path::new(path);
            let mut file = File::create(path)?;
            let json = serde_json::to_string(&roc)?;
            file.write_all(json.as_bytes())?;
        }
        None => println!("ROC:\n{:?}", roc),
    };
    if matches.is_present("auc") {
        println!("AUC: {}", analysis::calc_auc(&roc));
    }
    Ok(())
}

fn exec_eval_sequences<'a>(matches: &'a ArgMatches) -> Result<(), Box<Error>> {
    let tis_position = 100;
    let input = matches.value_of("INPUT").unwrap();
    let sequences = util::read_sequences_or_exit(input);
    let threshold: f64 = matches.value_of("threshold").unwrap().parse()?;
    let pwm_path = matches.value_of("pwm").unwrap();
    let pwm = PWM::load(pwm_path)?;
    println!(
        "{}",
        pwm.eval_sequences(&sequences, threshold, tis_position)
    );
    Ok(())
}

fn exec_count_auc_pairs(matches: &ArgMatches) -> Result<(), Box<Error>> {
    let input_train = matches.value_of("INPUT_TRAIN").unwrap();
    let sequences_train = util::read_sequences_or_exit(input_train);
    let sequence_test_input = match matches.value_of("INPUT_TEST") {
        Some(path) => Some(util::read_sequences_or_exit(path)),
        None => None,
    };
    let sequences_test_ref = match sequence_test_input {
        Some(ref seq) => seq,
        None => &sequences_train,
    };
    let tis_position: usize = matches.value_of("tis_position").unwrap().parse()?;
    let offset: i32 = matches.value_of("offset").unwrap().parse()?;
    let length: usize = matches.value_of("length").unwrap().parse()?;
    let pseudocount_start: f64 = matches.value_of("pseudocount_start").unwrap().parse()?;
    let pseudocount_end: f64 = matches.value_of("pseudocount_end").unwrap().parse()?;
    let pseudocount_step: f64 = matches.value_of("pseudocount_step").unwrap().parse()?;
    let background_model: [f64; 4] = match matches.values_of("background_model") {
        Some(values) => {
            let bg_in: Vec<f64> = values.map(|s| s.parse::<f64>()).collect::<Result<_, _>>()?;
            [bg_in[0], bg_in[1], bg_in[2], bg_in[3]]
        }
        None => analysis::calc_background_model(&sequences_train, tis_position, length),
    };
    let count_auc_points = analysis::eval_auc_for_pseudocount(
        &sequences_train,
        sequences_test_ref,
        pseudocount_start,
        pseudocount_end,
        pseudocount_step,
        tis_position,
        offset,
        length,
        background_model,
    )?;
    match matches.value_of("output") {
        Some(path) => {
            let path = Path::new(path);
            let mut file = File::create(path)?;
            let json = serde_json::to_string(&count_auc_points)?;
            file.write_all(json.as_bytes())?;
        }
        None => println!("(count, AUC):\n{:?}", count_auc_points),
    };
    Ok(())
}

fn print_and_exit<T: Display>(err: T) -> ! {
    eprintln!("Aborting!\n{}", err);
    process::exit(1);
}
