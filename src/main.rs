extern crate getopts;
extern crate postgres;
extern crate xml;

extern crate serde;
extern crate serde_json;
#[macro_use]
extern crate serde_derive;

use std::{env, process};
use getopts::Options;

use std::fs::File;
use std::io::{Write, BufReader, BufWriter};
use std::collections::hash_map::HashMap;
use std::collections::HashSet;

use xml::reader::{EventReader, XmlEvent};
use xml::attribute::OwnedAttribute;
use xml::ParserConfig;

use postgres::{Connection, TlsMode};

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Serialize, Debug, Clone)]
struct Location {
    start: usize,
    end: usize,
    score: f32,
}

#[derive(Serialize, Debug, Clone)]
struct Match {
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
struct UniprotResult {
    uniprot_id: String,
    interpro_matches: Vec<Match>,
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

fn add_ipr_to_match(interpro_match: &mut Match, attributes: &Vec<OwnedAttribute>) {
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

    interpro_match.interpro_id = interpro_id.unwrap();
    interpro_match.interpro_name = interpro_name.unwrap();
    interpro_match.interpro_type = interpro_type.unwrap();
}

fn add_lcn_to_match(interpro_match: &mut Match, attributes: &Vec<OwnedAttribute>) {
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

    interpro_match.locations
        .push(Location {
            start: start.unwrap(),
            end: end.unwrap(),
            score: score.unwrap(),
        });
}

fn make_uniprot_result(chado_uniprot_ids: &HashSet<String>,
                       attributes: &Vec<OwnedAttribute>) -> Option<UniprotResult>
{
    let mut uniprot_id = None;

    for attr in attributes {
        if attr.name.local_name == "id" {
            uniprot_id = Some(attr.value.clone());
            break;
        }
    }

    if chado_uniprot_ids.get(uniprot_id.as_ref().unwrap()).is_some() {
        Some(UniprotResult {
            uniprot_id: uniprot_id.unwrap(),
            interpro_matches: vec![],
        })
    } else {
        None
    }
}

fn make_match(attributes: &Vec<OwnedAttribute>) -> Option<Match> {
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
        Some(Match {
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
         -> HashMap<String, UniprotResult>
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
    let mut interpro_match = None;

    for e in parser {
        match e.unwrap() {
            XmlEvent::StartElement { ref name, ref attributes, .. } => {
                match &*name.local_name as &str {
                    "protein" => {
                        uniprot_result = make_uniprot_result(chado_uniprot_ids, attributes);
                    },
                    "match" => {
                        interpro_match = make_match(attributes);
                    },
                    "ipr" => {
                        add_ipr_to_match(interpro_match.as_mut().unwrap(), attributes);
                    },
                    "lcn" => {
                        add_lcn_to_match(interpro_match.as_mut().unwrap(), attributes);
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
                            let interpro_matches =
                                &mut uniprot_result.interpro_matches;
                            interpro_matches.push(interpro_match.unwrap());
                            interpro_match = None;
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

    let chado_uniprot_ids = get_chado_uniprot_ids(&conn);

    let results = parse(&chado_uniprot_ids, &input_filename);

    let s = serde_json::to_string(&results).unwrap();
    let f = File::create(output_filename).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");
}
