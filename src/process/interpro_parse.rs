use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use csv;

use crate::types::{GeneMatches, InterProMatch, Location};

#[derive(Debug, Deserialize)]
pub struct InterProtRow {
    pub protein_id: String,
    pub md5: String,
    pub seq_len: String,
    pub analysis: String,
    pub signature_accession: String,
    pub signature_description: String,
    pub start: String,
    pub end: String,
    pub score: String,
    pub status: String,
    pub date: String,
    pub interpro_accession: String,
    pub interpro_description: String,
    pub go_annotations: String,
    pub pathways_annotations: String,
}


/// Parse an InterPro TSV file.  Return a map from UniProt ID to struct
/// containing its InterProMatches.
pub fn parse(filename: &str)
         -> HashMap<String, GeneMatches>
{
    let file = File::open(filename).unwrap();
    let file = BufReader::new(file);

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(Box::new(file));

    let headers = csv::StringRecord::from(vec![
        "protein_id",
        "md5",
        "seq_len",
        "analysis",
        "signature_accession",
        "signature_description",
        "start",
        "end",
        "score",
        "status",
        "date",
        "interpro_accession",
        "interpro_description",
        "go_annotations",
        "pathways_annotations",
    ]);

    let mut match_map = HashMap::new();

    for record in reader.records() {
        let mut record: InterProtRow = record.unwrap().deserialize(Some(&headers)).unwrap();

        if record.signature_description == "-" {
           record.signature_description = "".into();
        }
        if record.interpro_accession == "-" {
           record.interpro_accession = "".into();
        }
        if record.interpro_description == "-" {
           record.interpro_description = "".into();
        }

        let gene_uniquename = record.protein_id.replace(".1:pep", "");

        let match_id = record.signature_accession.clone();

        match_map
            .entry(gene_uniquename)
            .or_insert_with(HashMap::new)
            .entry(match_id)
            .or_insert_with(|| InterProMatch {
                id: record.signature_accession.clone(),
                dbname: record.analysis.clone(),
                name: record.signature_description.clone(),
                interpro_id: record.interpro_accession.clone(),
                interpro_name: record.interpro_description.clone(),
                locations: vec![],
            })
            .locations
            .push(Location {
                start: record.start.parse::<usize>().unwrap(),
                end: record.end.parse::<usize>().unwrap(),
            });
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

    results
}
