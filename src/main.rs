extern crate getopts;
extern crate regex;

extern crate serde;
extern crate serde_json;
extern crate serde_derive;

use std::ffi::OsString;
use std::{env, process};
use getopts::Options;
use std::thread;
use std::thread::JoinHandle;
use std::process::Command;
use regex::Regex;

use std::fs::File;
use std::io::{Write, BufReader, BufRead, BufWriter};
use std::collections::HashMap;

extern crate domain_process;

use domain_process::{segmasker_parse, types::*};
use domain_process::interpro_parse::parse;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn make_tmhmm_thread(protein_file_name: &str)
                     -> JoinHandle<HashMap<String, Vec<TMMatch>>>
{
    let protein_file_name_ostring: OsString = protein_file_name.into();
    let re = Regex::new(r"(?i)(\S+)\s+tmhmm\S+\s+tmhelix\s+(\d+)\s+(\d+)").unwrap();

    thread::spawn(move || {
        let mut ret = HashMap::new();
        let tmhmm = Command::new("tmhmm")
            .arg(protein_file_name_ostring)
            .stdout(process::Stdio::piped())
            .spawn()
            .expect("tmhmm command failed to start");
        let buf_reader = BufReader::new(tmhmm.stdout.unwrap());
        'LINE: for line_result in buf_reader.lines() {
            let line = line_result.unwrap();
            if line.starts_with("#") {
                continue 'LINE;
            }

            let re_result = re.captures(&line);

            if let Some(captures) = re_result {
                let uniprot_id = captures.get(1).unwrap().as_str();
                let start = captures.get(2).unwrap().as_str().parse::<usize>().unwrap();
                let end = captures.get(3).unwrap().as_str().parse::<usize>().unwrap();
                ret.entry(String::from(uniprot_id))
                    .or_insert(vec![])
                    .push(TMMatch {
                        start,
                        end,
                    });
            }
        }
        ret
    })
}



fn make_segmasker_thread(protein_file_name: &str)
        -> JoinHandle<HashMap<String, Vec<Location>>>
{
    let protein_file_name_ostring: OsString = protein_file_name.into();

    thread::spawn(move || {

        let tmhmm = Command::new("segmasker")
            .arg("-in")
            .arg(protein_file_name_ostring)
            .stdout(process::Stdio::piped())
            .spawn()
            .expect("segmasker command failed to start");
        let mut buf_reader = BufReader::new(tmhmm.stdout.unwrap());
        segmasker_parse::parse(&mut buf_reader)
    })
}


/// Parse the InterPro XML and run TMHMM to create a JSON file for the PomBase
/// front end to display.
fn main() -> Result<(), std::io::Error> {
    print!("{} v{}\n", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.optflag("h", "help", "print this help message");
    opts.optopt("p", "postgresql-connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.optopt("i", "input-file",
                "Input InterPro XML file", "FILE");
    opts.optopt("o", "output-file",
                "Output JSON file", "FILE");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let program = args[0].clone();

    if matches.opt_present("help") {
        print_usage(&program, opts);
        process::exit(0);
    }

    if !matches.opt_present("postgresql-connection-string") {
        print!("no -p|--postgresql-connection-string option\n");
        print_usage(&program, opts);
        process::exit(1);
    }

    let protein_filename = matches.opt_str("p").unwrap();
    let input_filename = matches.opt_str("i").unwrap();
    let output_filename = matches.opt_str("o").unwrap();

    let (interproscan_version, mut domains_by_id) = parse(&input_filename);

    let tmhmm_handle = make_tmhmm_thread(&protein_filename);

    let tmhmm_matches =
        tmhmm_handle.join().expect("Failed to get TMHMM results");

    for (protein_id, domain_match) in tmhmm_matches {
        let gene_uniquename = protein_id.replace(".1:pep", "");
        domains_by_id.entry(gene_uniquename.clone())
            .or_insert(GeneMatches {
                gene_uniquename,
                interpro_matches: vec![],
                segmasker_matches: vec![],
                tmhmm_matches: vec![],
            })
            .tmhmm_matches.extend(domain_match.into_iter());
    }
 
    let segmasker_handle = make_segmasker_thread(&protein_filename);

    let segmasker_matches = segmasker_handle.join().expect("Failed to run segmasker");

    for (gene_uniquename, locations) in segmasker_matches {

        domains_by_id.entry(gene_uniquename.clone())
            .or_insert(GeneMatches {
                gene_uniquename,
                interpro_matches: vec![],
                segmasker_matches: vec![],
                tmhmm_matches: vec![],
            })
            .segmasker_matches.extend(locations.into_iter());
    }

    let domain_data = DomainData {
        interproscan_version,
        domains_by_id,
    };

    let s = serde_json::to_string(&domain_data).unwrap();
    let f = File::create(output_filename).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    Ok(())
}
