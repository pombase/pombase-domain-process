use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use crate::types::{GeneMatches, InterProMatch, Location};

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
pub struct InterProScanLocation {
    pub start: usize,
    pub end: usize,
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

    let interproscan_output: InterProScanOutput =
        match serde_json::from_reader(reader) {
            Ok(config) => config,
            Err(err) => {
                panic!("failed to parse {}: {}", filename, err)
            },
        };

    let mut match_map: HashMap<String, HashMap<String, InterProMatch>> = HashMap::new();

    for result in interproscan_output.results.into_iter() {
        for interpro_match in result.matches.into_iter() {

            let gene_uniquename = result.xref.get(0).unwrap().id.replace(".1:pep", "");

            let first_location = &interpro_match.locations[0];
            let sequence_feature_str =
                if let Some(ref sequence_feature) = first_location.sequence_feature {
                    if sequence_feature.len() > 0 {
                        format!("-{}", sequence_feature.replace(" ", "-"))
                    } else {
                        "".into()
                    }
                } else {
                    "".into()
                };

            let match_id = format!("{}{}", interpro_match.signature.accession,
                                   sequence_feature_str);

            let locations: Vec<_> = interpro_match.locations.iter()
                .map(|loc| {
                    Location {
                        start: loc.start,
                        end: loc.end,
                    }
                })
                .collect();

            match_map
                .entry(gene_uniquename)
                .or_insert_with(HashMap::new)
                .entry(match_id.clone())
                .or_insert_with(|| {
                    let signature = &interpro_match.signature;
                    let dbname = format!("{}{}", signature.library_release.library,
                                         sequence_feature_str);
                    let interpro_id =
                        if let Some(ref entry) = signature.entry {
                            entry.accession.clone()
                        } else {
                            "".into()
                        };
                    let interpro_name =
                        if let Some(ref entry) = signature.entry {
                            entry.name.clone()
                        } else {
                            "".into()
                        };
                    let interpro_description =
                        if let Some(ref entry) = signature.entry {
                            entry.description.clone()
                        } else {
                            "".into()
                        };
                    InterProMatch {
                        id: match_id.clone(),
                        dbname,
                        name: signature.name.clone(),
                        description: signature.description.clone(),
                        interpro_id,
                        interpro_name,
                        interpro_description,
                        locations: vec![],
                    }
                })
                .locations
                .extend(locations.clone().into_iter());
        }

    }

    let mut results = HashMap::new();

    for (gene_uniquename, domains_by_id) in match_map.into_iter() {
        for mut interpro_match in domains_by_id.into_values() {
            interpro_match.locations.sort();
            results
                .entry(gene_uniquename.clone())
                .or_insert_with(|| GeneMatches {
                    gene_uniquename: gene_uniquename.clone(),
                    interpro_matches: vec![],
                    tmhmm_matches: vec![],
                })
                .interpro_matches
                .push(interpro_match);
        }

        if let Some(ref mut gene_matches) = results.get_mut(&gene_uniquename) {
            gene_matches.interpro_matches
                .sort_by(|a, b| { a.dbname.cmp(&b.dbname) });
        }
    }

    (interproscan_output.interproscan_version, results)
}
