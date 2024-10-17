use std::collections::HashMap;
use std::io::BufRead;

use regex::Regex;

use crate::types::Location;


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
