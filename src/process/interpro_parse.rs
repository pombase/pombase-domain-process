use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::cmp::Ordering;

use crate::types::{GeneMatches, InterProMatch, Location};
use crate::util::merge_locations;

#[derive(Debug, Deserialize)]
pub struct InterProScanOutput {
    #[serde(rename = "interproscan-version")]
    pub interproscan_version: String,
    pub results: Vec<InterProScanResult>,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanResult {
    pub matches: Vec<InterProScanMatch>,
    pub xref: Vec<InterProScanXref>,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanMatch {
    pub signature: InterProScanSignature,
    pub locations: Vec<InterProScanLocation>,
    #[serde(rename = "model-ac")]
    pub model_ac: String,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanSignature {
    pub accession: String,
    pub name: Option<String>,
    pub description: Option<String>,
    #[serde(rename = "signatureLibraryRelease")]
    pub library_release: InterProScanSignatureLibraryRelease,
    pub entry: Option<InterProScanEntry>,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanLocationFragment {
    pub start: usize,
    pub end: usize,
    #[serde(rename = "dc-status")]
    pub dc_status: String,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanLocation {
    pub start: usize,
    pub end: usize,
    #[serde(rename = "location-fragments")]
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub location_fragments: Vec<InterProScanLocationFragment>,
    #[serde(rename = "sequence-feature")]
    pub sequence_feature: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanSignatureLibraryRelease {
    pub library: String,
    pub version: String,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanXref {
    pub id: String,
}

#[derive(Debug, Deserialize)]
pub struct InterProScanEntry {
    pub accession: String,
    pub name: String,
    pub description: String,
    #[serde(rename = "type")]
    pub entry_type: String,
}

pub type VersionString = String;


/// Parse an InterPro TSV file.  Return a map from UniProt ID to struct
/// containing its InterProMatches.
pub fn parse(filename: &str)
         -> (VersionString, HashMap<String, GeneMatches>)
{
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(err) => {
            panic!("Failed to read {}: {}\n", filename, err)
        }
    };

    let reader = BufReader::new(file);

    let mut interproscan_output: InterProScanOutput =
        match serde_json::from_reader(reader) {
            Ok(config) => config,
            Err(err) => {
                panic!("failed to parse {}: {}", filename, err)
            },
        };

    for result in interproscan_output.results.iter_mut() {
        for interpro_match in result.matches.iter_mut() {
            if let Some(ref name) = interpro_match.signature.name {
                if name == "" {
                    interpro_match.signature.name = None
                }
            }
        }
    }

    let mut gene_match_map: HashMap<String, HashMap<String, InterProMatch>> = HashMap::new();

    for result in interproscan_output.results.into_iter() {
        let gene_uniquename = result.xref.get(0).unwrap().id.replace(".1:pep", "");

        let mut match_map: HashMap<String, InterProMatch> = HashMap::new();

        for interpro_match in result.matches.into_iter() {
            let signature = &interpro_match.signature;
            let library = &signature.library_release.library.replace("MOBIDB_LITE", "MOBIDB");

            for loc in interpro_match.locations {

                let sequence_feature_str =
                    if let Some(ref sequence_feature) = loc.sequence_feature {
                        if sequence_feature.len() > 0 {
                            format!("-{}", sequence_feature.replace(" ", "-"))
                        } else {
                            if library == "MOBIDB" {
                                "-Disorder".into()
                            } else {
                                "".into()
                            }
                        }
                    } else {
                        "".into()
                    };

                let match_id = format!("{}{}", signature.accession, sequence_feature_str);

                let fragment_locs =
                    if loc.location_fragments.len() > 0 {
                        loc.location_fragments
                            .iter()
                            .map(|frag| Location {
                                start: frag.start,
                                end: frag.end,
                            })
                            .collect()
                    } else {
                        vec![Location {
                            start: loc.start,
                            end: loc.end,
                        }]
                    };

                match_map
                    .entry(match_id.clone())
                    .or_insert_with(|| {
                        let dbname = format!("{}{}", library,
                                             sequence_feature_str);
                        let interpro_id =
                            if let Some(ref entry) = signature.entry {
                                Some(entry.accession.clone())
                            } else {
                                None
                            };
                        let interpro_name =
                            if let Some(ref entry) = signature.entry {
                                Some(entry.name.clone())
                            } else {
                                None
                            };
                        let interpro_description =
                            if let Some(ref entry) = signature.entry {
                                Some(entry.description.clone())
                            } else {
                                None
                            };
                        InterProMatch {
                            id: match_id.clone(),
                            dbname,
                            name: signature.name.clone(),
                            description: signature.description.clone(),
                            interpro_id,
                            interpro_name,
                            interpro_description,
                            match_start: usize::MAX,
                            match_end: 0,
                            locations: vec![],
                        }
                    })
                    .locations
                    .extend(fragment_locs.into_iter());
            }
        }

        for interpro_match in match_map.values_mut() {
            for loc in interpro_match.locations.iter() {
                if loc.start < interpro_match.match_start {
                    interpro_match.match_start = loc.start;
                }
                if loc.end > interpro_match.match_end {
                    interpro_match.match_end = loc.end;
                }
            }
        }

        gene_match_map.insert(gene_uniquename, match_map);
    }

    let mut results = HashMap::new();

    for (gene_uniquename, domains_by_id) in gene_match_map.into_iter() {
        for mut interpro_match in domains_by_id.into_values() {
            interpro_match.locations.sort();
            merge_locations(&mut interpro_match.locations);
            results
                .entry(gene_uniquename.clone())
                .or_insert_with(|| GeneMatches {
                    gene_uniquename: gene_uniquename.clone(),
                    interpro_matches: vec![],
                    segmasker_matches: vec![],
                    tmhmm_matches: vec![],
                })
                .interpro_matches
                .push(interpro_match);
        }

        if let Some(ref mut gene_matches) = results.get_mut(&gene_uniquename) {
            gene_matches.interpro_matches
                .sort_by(|a, b| {
                    let dbname_cmp = a.dbname.cmp(&b.dbname);
                    if dbname_cmp == Ordering::Equal {
                        Ordering::Equal
                    } else {
                        if a.dbname.to_ascii_lowercase().starts_with("pfam") {
                            Ordering::Less
                        } else {
                            if b.dbname.to_ascii_lowercase().starts_with("pfam") {
                                Ordering::Greater
                            } else {
                                dbname_cmp
                            }
                        }
                    }
                });
        }
    }

    (interproscan_output.interproscan_version, results)
}
