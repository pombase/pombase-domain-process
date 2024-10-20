
use std::{cmp::Ordering, collections::HashMap, hash::{Hash, Hasher}};

#[derive(Serialize, Debug, Clone)]
pub struct Location {
    pub start: usize,
    pub end: usize,
}

impl PartialEq for Location {
    fn eq(&self, other: &Location) -> bool {
        self.start == other.start && self.end == other.end
    }
}
impl Eq for Location { }
impl Ord for Location {
    fn cmp(&self, other: &Location) -> Ordering {
        let start_cmp = self.start.cmp(&other.start);

        if start_cmp == Ordering::Equal {
            self.end.cmp(&other.end)
        } else {
            start_cmp
        }
    }
}
 impl PartialOrd for Location {
     fn partial_cmp(&self, other: &Location) -> Option<Ordering> {
         Some(self.cmp(other))
     }
 }
 impl Hash for Location {
     fn hash<H: Hasher>(&self, state: &mut H) {
         self.start.hash(state);
         self.end.hash(state);
     }
 }


#[derive(Serialize, Debug, Clone)]
pub struct InterProMatch {
    pub id: String,
    pub dbname: String,
    pub name: Option<String>,
    pub description: Option<String>,
    pub interpro_id: String,
    pub interpro_name: String,
    pub interpro_description: String,
    pub match_start: usize,
    pub match_end: usize,
    pub locations: Vec<Location>,
}

#[derive(Serialize, Debug, Clone)]
pub struct TMMatch {
    pub start: usize,
    pub end: usize,
}

#[derive(Serialize, Debug, Clone)]
pub struct GeneMatches {
    pub gene_uniquename: String,
    pub interpro_matches: Vec<InterProMatch>,
    pub segmasker_matches: Vec<Location>,
    pub tmhmm_matches: Vec<TMMatch>,
}

#[derive(Serialize, Debug, Clone)]
pub struct DomainData {
    pub interproscan_version: String,
    pub domains_by_id: HashMap<String, GeneMatches>,
}
