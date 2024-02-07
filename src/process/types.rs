
use std::collections::HashMap;

#[derive(Serialize, Debug, Clone)]
pub struct LocationScore {
    pub start: usize,
    pub end: usize,
    pub score: f32,
}

#[derive(Serialize, Debug, Clone)]
pub struct InterProMatch {
    pub id: String,
    pub dbname: String,
    pub name: String,
    pub model: Option<String>,
    pub evidence: String,
    pub interpro_id: String,
    pub interpro_name: String,
    pub interpro_type: String,
    pub locations: Vec<LocationScore>,
}

#[derive(Serialize, Debug, Clone)]
pub struct TMMatch {
    pub start: usize,
    pub end: usize,
}

#[derive(Serialize, Debug, Clone)]
pub struct UniProtResult {
    pub uniprot_id: String,
    pub interpro_matches: Vec<InterProMatch>,
    pub tmhmm_matches: Vec<TMMatch>,
}

#[derive(Serialize, Debug, Clone)]
pub struct DomainData {
    pub interpro_version: String,
    pub domains_by_id: HashMap<String, UniProtResult>,
}
