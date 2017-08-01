extern crate getopts;
extern crate postgres;
extern crate xml;
extern crate regex;

extern crate serde;
extern crate serde_json;
#[macro_use]
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
use std::collections::hash_map::HashMap;
use std::collections::HashSet;

use xml::reader::{EventReader, XmlEvent};
use xml::attribute::OwnedAttribute;
use xml::ParserConfig;

use postgres::{Connection, TlsMode};

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

type UniProtId = String;

#[derive(Serialize, Debug, Clone)]
struct Location {
    start: usize,
    end: usize,
    score: f32,
}

#[derive(Serialize, Debug, Clone)]
struct InterproMatch {
    id: String,
    dbname: String,
    name: String,
    evidence: String,
    interpro_id: String,
    interpro_name: String,
    interpro_type: String,
    locations: Vec<Location>,
}

#[derive(Serialize, Debug, Clone)]
struct TMMatch {
    start: usize,
    end: usize,
}

#[derive(Serialize, Debug, Clone)]
struct UniProtResult {
    uniprot_id: String,
    interpro_matches: Vec<InterproMatch>,
    tmhmm_matches: Vec<TMMatch>,
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn get_chado_uniprot_ids(conn: &Connection) -> HashSet<String> {
    let mut return_set = HashSet::new();

    for row in &conn.query("select value from featureprop p join cvterm t on p.type_id = t.cvterm_id where name = 'uniprot_identifier'", &[]).unwrap() {
        let uniprot_id: String = row.get(0);
        return_set.insert(uniprot_id);
    }

    return_set
}

fn get_chado_prot_sequences(conn: &Connection) -> HashMap<String, String> {
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

fn write_prot_sequences(path: &Path, seq_map: &HashMap<String, String>) {
    let f = File::create(path).expect("Can't create protein sequence file");
    let mut writer = BufWriter::new(&f);
    for (gene_uniquename, residues) in seq_map {
        writer.write_fmt(format_args!(">{}\n", gene_uniquename)).unwrap();
        writer.write_fmt(format_args!("{}\n", residues)).unwrap();
    }
}

fn make_tmhmm_thread(conn: &Connection) -> JoinHandle<HashMap<UniProtId, Vec<TMMatch>>> {
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

fn add_ipr_to_match(domain_match: &mut InterproMatch, attributes: &Vec<OwnedAttribute>) {
    let mut interpro_id = None;
    let mut interpro_name = None;
    let mut interpro_type = None;

    for attr in attributes {
        match &attr.name.local_name as &str {
            "id" => {
                interpro_id = Some(attr.value.clone());
            },
            "name" => {
                interpro_name = Some(attr.value.clone());
            },
            "type" => {
                interpro_type = Some(attr.value.clone());
            },
            _ => (),
        }
    }

    domain_match.interpro_id = interpro_id.unwrap();
    domain_match.interpro_name = interpro_name.unwrap();
    domain_match.interpro_type = interpro_type.unwrap();
}

fn add_lcn_to_match(domain_match: &mut InterproMatch, attributes: &Vec<OwnedAttribute>) {
    let mut start = None;
    let mut end = None;
    let mut score = None;

    for attr in attributes {
        match &attr.name.local_name as &str {
            "start" => {
                start = Some(attr.value.parse::<usize>().unwrap());
            },
            "end" => {
                end = Some(attr.value.parse::<usize>().unwrap());
            },
            "score" => {
                score = Some(attr.value.parse::<f32>().unwrap());
            },
            _ => (),
        }
    }

    domain_match.locations
        .push(Location {
            start: start.unwrap(),
            end: end.unwrap(),
            score: score.unwrap(),
        });
}

fn make_uniprot_result(chado_uniprot_ids: &HashSet<String>,
                       attributes: &Vec<OwnedAttribute>) -> Option<UniProtResult>
{
    let mut uniprot_id = None;

    for attr in attributes {
        if attr.name.local_name == "id" {
            uniprot_id = Some(attr.value.clone());
            break;
        }
    }

    if chado_uniprot_ids.get(uniprot_id.as_ref().unwrap()).is_some() {
        Some(UniProtResult {
            uniprot_id: uniprot_id.unwrap(),
            interpro_matches: vec![],
            tmhmm_matches: vec![],
        })
    } else {
        None
    }
}

fn make_match(attributes: &Vec<OwnedAttribute>) -> Option<InterproMatch> {
    let mut id = None;
    let mut name = None;
    let mut dbname = None;
    let mut status = None;
    let mut evidence = None;

    for attr in attributes {
        match &attr.name.local_name as &str {
            "id" => {
                id = Some(attr.value.clone());
            },
            "name" => {
                name = Some(attr.value.clone());
            },
            "dbname" => {
                dbname = Some(attr.value.clone());
            },
            "status" => {
                status = Some(attr.value.clone());
            },
            "evd" => {
                evidence = Some(attr.value.clone());
            },
            _ => {
                panic!("unknown attribute name: {}", attr.name.local_name);
            }
        }
    }

    if status.is_some() {
        Some(InterproMatch {
            id: id.unwrap(),
            dbname: dbname.unwrap(),
            name: name.unwrap(),
            evidence: evidence.unwrap(),
            interpro_id: "".into(),
            interpro_name: "".into(),
            interpro_type: "".into(),
            locations: vec![],
        })
    } else {
        None
    }
}

fn parse(chado_uniprot_ids: &HashSet<String>, filename: &str)
         -> HashMap<String, UniProtResult>
{
    let file = File::open(filename).unwrap();
    let file = BufReader::new(file);
    let parser = EventReader::new_with_config(file, ParserConfig {
        trim_whitespace: true,
        cdata_to_characters: true,
        ..Default::default()
    });

    let mut results = HashMap::new();

    let mut uniprot_result = None;
    let mut domain_match = None;

    for e in parser {
        match e.unwrap() {
            XmlEvent::StartElement { ref name, ref attributes, .. } => {
                match &*name.local_name as &str {
                    "protein" => {
                        uniprot_result = make_uniprot_result(chado_uniprot_ids, attributes);
                    },
                    "match" => {
                        domain_match = make_match(attributes);
                    },
                    "ipr" => {
                        add_ipr_to_match(domain_match.as_mut().unwrap(), attributes);
                    },
                    "lcn" => {
                        add_lcn_to_match(domain_match.as_mut().unwrap(), attributes);
                    },
                    _ => (),
                }
            },
            XmlEvent::EndElement { ref name, .. } => {
                match &*name.local_name as &str {
                    "protein" => {
                        if let Some(finished_uniprot_result) = uniprot_result {
                            results.insert(finished_uniprot_result.uniprot_id.clone(),
                                           finished_uniprot_result);
                            uniprot_result = None;
                        }
                    },
                    "match" => {
                        if let Some(ref mut uniprot_result) = uniprot_result {
                            let matches = &mut uniprot_result.interpro_matches;
                            matches.push(domain_match.unwrap());
                            domain_match = None;
                        }
                    },
                    _ => (),
                }
            },
            _ => ()
        }
    }

    results
}

fn main() {
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

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();

    let tmhmm_handle = make_tmhmm_thread(&conn);

    let chado_uniprot_ids = get_chado_uniprot_ids(&conn);
    let mut results = parse(&chado_uniprot_ids, &input_filename);

    let tmhmm_matches =
        tmhmm_handle.join().expect("Failed to get TMHMM results");

    for (uniprot_id, domain_match) in tmhmm_matches {
        results.entry(uniprot_id.clone())
            .or_insert(UniProtResult {
                uniprot_id: uniprot_id.clone(),
                interpro_matches: vec![],
                tmhmm_matches: vec![],
            })
            .tmhmm_matches.extend(domain_match.into_iter());
    }

    let s = serde_json::to_string(&results).unwrap();
    let f = File::create(output_filename).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");
}
