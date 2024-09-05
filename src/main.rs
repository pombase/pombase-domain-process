extern crate getopts;
extern crate postgres;
extern crate regex;

extern crate serde;
extern crate serde_json;
extern crate serde_derive;

extern crate tempfile;

use std::{env, process};
use getopts::Options;
use std::thread;
use std::thread::JoinHandle;
use std::process::Command;
use regex::Regex;

use tempfile::NamedTempFile;

use std::fs::File;
use std::path::Path;
use std::io::{Write, BufReader, BufRead, BufWriter};
use std::collections::{HashMap,HashSet};

use postgres::{Client, NoTls};

extern crate domain_process;

use domain_process::types::*;
use domain_process::interpro_parse::parse;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

type UniProtId = String;
fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

/// Return a HashSet of all the UniProt IDs of genes in Chado
fn get_chado_uniprot_ids(conn: &mut Client) -> HashSet<String> {
    let mut return_set = HashSet::new();

    for row in &conn.query("select value from featureprop p join cvterm t on p.type_id = t.cvterm_id where name = 'uniprot_identifier'", &[]).unwrap() {
        let uniprot_id: String = row.get(0);
        return_set.insert(uniprot_id);
    }

    return_set
}

/// Return a HashMap where the keys are the uniquenames of protein
/// coding genes and the values are the corresponding protein
/// sequence
fn get_chado_prot_sequences(conn: &mut Client) -> HashMap<String, String> {
    let mut ret = HashMap::new();

    for row in &conn.query("
SELECT uniprot_id_prop.value AS uniprot_id, prot.residues
FROM feature g
JOIN featureprop uniprot_id_prop ON g.feature_id = uniprot_id_prop.feature_id
JOIN feature_relationship gt_rel ON gt_rel.object_id = g.feature_id
JOIN feature tr ON gt_rel.subject_id = tr.feature_id
JOIN feature_relationship tp_rel ON tr.feature_id = tp_rel.object_id
JOIN feature prot ON tp_rel.subject_id = prot.feature_id
WHERE gt_rel.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'part_of')
  AND tp_rel.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'derives_from')
  AND prot.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'polypeptide')
  AND uniprot_id_prop.type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'uniprot_identifier')
", &[]).unwrap() {
        let gene_uniquename: String = row.get(0);
        let residues: String = row.get(1);
        ret.insert(gene_uniquename, residues);
    }

    ret
}


/// Write a FASTA file of sequences to the given Path
fn write_prot_sequences(path: &Path, seq_map: &HashMap<String, String>) {
    let f = File::create(path).expect("Can't create protein sequence file");
    let mut writer = BufWriter::new(&f);
    for (gene_uniquename, residues) in seq_map {
        writer.write_fmt(format_args!(">{}\n", gene_uniquename)).unwrap();
        writer.write_fmt(format_args!("{}\n", residues)).unwrap();
    }
}


/// Fist connect to Postgres to read the UniProt IDs and protein sequences
/// of all genes, then write them to a FASTA file.
/// Next spawn a thread for running TMHMM on that file.
/// We parse the results to make a map of TMMatches for each UniProt ID.
fn make_tmhmm_thread(conn: &mut Client) -> JoinHandle<HashMap<UniProtId, Vec<TMMatch>>> {
    let seq_map = get_chado_prot_sequences(conn);
    let tmpfile = NamedTempFile::new().unwrap();
    write_prot_sequences(tmpfile.path(), &seq_map);

    let re = Regex::new(r"(?i)(\S+)\s+tmhmm\S+\s+tmhelix\s+(\d+)\s+(\d+)").unwrap();

    thread::spawn(move || {
        let mut ret = HashMap::new();
        let tmhmm = Command::new("tmhmm")
            .arg(tmpfile.path().as_os_str())
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
                        start: start,
                        end: end,
                    });
            }
        }
        ret
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

    let connection_string = matches.opt_str("p").unwrap();
    let input_filename = matches.opt_str("i").unwrap();
    let output_filename = matches.opt_str("o").unwrap();

    let mut conn = Client::connect(connection_string.as_str(), NoTls).unwrap();

    let tmhmm_handle = make_tmhmm_thread(&mut conn);

    let chado_uniprot_ids = get_chado_uniprot_ids(&mut conn);
    let mut domain_data = parse(&chado_uniprot_ids, &input_filename);

    let tmhmm_matches =
        tmhmm_handle.join().expect("Failed to get TMHMM results");

    for (uniprot_id, domain_match) in tmhmm_matches {
        domain_data.domains_by_id.entry(uniprot_id.clone())
            .or_insert(UniProtResult {
                uniprot_id: uniprot_id.clone(),
                interpro_matches: vec![],
                tmhmm_matches: vec![],
            })
            .tmhmm_matches.extend(domain_match.into_iter());
    }

    let s = serde_json::to_string(&domain_data).unwrap();
    let f = File::create(output_filename).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    Ok(())
}
