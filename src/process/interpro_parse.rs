extern crate xml;

use self::xml::reader::{EventReader, XmlEvent};
use self::xml::attribute::OwnedAttribute;
use self::xml::ParserConfig;

use std::collections::{HashMap,HashSet};
use std::fs::File;
use std::io::BufReader;

use types::*;

/// Add the InterPro name, ID and type from the current <match> element to
/// the InterProMatch
fn add_ipr_to_match(interpro_match: &mut InterProMatch, attributes: &Vec<OwnedAttribute>) {
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


/// Add LocationScores to the match
fn add_lcn_to_match(interpro_match: &mut InterProMatch, attributes: &Vec<OwnedAttribute>) {
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
        .push(LocationScore {
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


/// Use the attributes to make an InterProMatch object.  The InterPro name,
/// ID and type and the LocationScore objects will be added later.
fn make_match(attributes: &Vec<OwnedAttribute>) -> Option<InterProMatch> {
    let mut id = None;
    let mut name = None;
    let mut dbname = None;
    let mut status = None;
    let mut evidence = None;
    let mut model = None;

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
            "model" => {
                model = Some(attr.value.clone());
            }
            _ => {
                panic!("unknown attribute name: {}", attr.name.local_name);
            }
        }
    }

    if status.is_some() {
        Some(InterProMatch {
            id: id.unwrap(),
            dbname: dbname.unwrap(),
            name: name.unwrap(),
            model: model,
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


/// Parse a <dbinfo> element and return the version attribute if
/// the dbname is "INTERPRO"
fn get_version(attributes: &Vec<OwnedAttribute>) -> Option<String> {
    let mut dbname = None;
    let mut version = None;

    for attr in attributes {
        match &attr.name.local_name as &str {
            "dbname" => dbname = Some(attr.value.clone()),
            "version" => version = Some(attr.value.clone()),
            _ => (),
        }
    }

    if let Some(ref dbname) = dbname {
        if dbname != "INTERPRO" {
            return None;
        }
    } else {
        panic!("missing dbname attribute in <dbinfo> element");
    }

    if let Some(version) = version {
        Some(version)
    } else {
        panic!("missing version attribute in <dbinfo> element");
    }
}


/// Parse an InterPro XML file.  Return a map from UniProt ID to struct
/// containing its InterProMatches.
pub fn parse(chado_uniprot_ids: &HashSet<String>, filename: &str)
         -> DomainData
{
    let file = File::open(filename).unwrap();
    let file = BufReader::new(file);
    let parser = EventReader::new_with_config(file, ParserConfig {
        trim_whitespace: true,
        cdata_to_characters: true,
        ..Default::default()
    });

    let mut domains_by_id = HashMap::new();

    let mut uniprot_result = None;
    let mut interpro_match = None;
    let mut interpro_version = None;

    for e in parser {
        match e.unwrap() {
            XmlEvent::StartElement { ref name, ref attributes, .. } => {
                match &name.local_name as &str {
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
                    "dbinfo" => {
                        if interpro_version.is_none() {
                            if let Some(parsed_interpro_version) = get_version(attributes)
                            {
                                interpro_version = Some(parsed_interpro_version);
                            }
                        }
                    },
                    _ => (),
                }
            },
            XmlEvent::EndElement { ref name, .. } => {
                match &name.local_name as &str {
                    "protein" => {
                        if let Some(finished_uniprot_result) = uniprot_result {
                            domains_by_id.insert(finished_uniprot_result.uniprot_id.clone(),
                                           finished_uniprot_result);
                            uniprot_result = None;
                        }
                    },
                    "match" => {
                        if let Some(ref mut uniprot_result) = uniprot_result {
                            let matches = &mut uniprot_result.interpro_matches;
                            matches.push(interpro_match.unwrap());
                            interpro_match = None;
                        }
                    },
                    _ => (),
                }
            },
            _ => ()
        }
    }

    let interpro_version =
        interpro_version.expect("failing to find a <dbinfo> element where dbname=\"INTERPRO\"");

    DomainData {
        interpro_version,
        domains_by_id,
    }
}
