use std::collections::HashMap;
use std::io::BufRead;

use regex::Regex;

use crate::types::Location;

// merge locations/ranges that abut or overlap
pub fn merge_locations(locations: &mut Vec<Location>)
{
    if locations.len() <= 1 {
        return;
    }

    locations.sort_by_key(|a| a.start);

    for i in (1..locations.len()).rev() {
        let (this_start, this_end) = {
            let this = locations.get(i).unwrap();
            (this.start, this.end)
        };
        let prev = locations.get_mut(i-1).unwrap();

        if this_start <= prev.end + 1 {
            if this_end > prev.end {
                prev.end = this_end;
            }
            locations.remove(i);
        }
    }
}

pub fn parse(buf_reader: &mut dyn BufRead)
    -> HashMap<String, Vec<Location>>
{
    let mut ret = HashMap::new();

    let gene_re = Regex::new(r"^>(\S+)\.1:pep\s.*").unwrap();

    let mut current_gene = String::from("");

    for line_result in buf_reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with(">") {
            let re_result = gene_re.captures(&line);
            let captures = re_result.unwrap();
            current_gene = captures.get(1).unwrap().as_str().to_owned();
        } else {
            let line_parts: Vec<_> = line.split(" ").collect();
            if *line_parts.get(1).unwrap() != "-" {
                panic!("can't parse line from segmasker: {}", line);
            }
            let start = line_parts.get(0).unwrap().parse::<usize>().unwrap();
            let end = line_parts.get(2).unwrap().parse::<usize>().unwrap();
            ret.entry(current_gene.clone())
                .or_insert(vec![])
                .push(Location {
                    start,
                    end,
                });
        }
    }
    ret
}
