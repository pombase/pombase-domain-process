extern crate domain_process;

use std::fs::File;
use std::io::BufReader;

use domain_process::interpro_parse;
use domain_process::segmasker;
use domain_process::segmasker::merge_locations;
use domain_process::types::Location;

#[test]
fn test_parse() {

    let (interproscan_version, matches) = interpro_parse::parse("tests/small_matches.json");

    assert_eq!(interproscan_version, "5.69-101.0");

    let spac13g6_15c = matches.get("SPAC13G6.15c").unwrap();
    assert_eq!(spac13g6_15c.gene_uniquename, "SPAC13G6.15c");

    let mobidb_match = spac13g6_15c.interpro_matches.get(0).unwrap();
    assert_eq!(mobidb_match.dbname, "MOBIDB-Disorder");
    assert_eq!(mobidb_match.id, "mobidb-lite-Disorder");

    let mobidb_locations = &mobidb_match.locations;
    assert_eq!(mobidb_locations.get(0).unwrap().end, 163);

    let panther_match = spac13g6_15c.interpro_matches.get(1).unwrap();
    assert_eq!(panther_match.dbname, "PANTHER");
    assert_eq!(panther_match.id, "PTHR10300");

    assert_eq!(panther_match.locations[0].end, 153);
}

#[test]
fn test_parse_segmasker() {
    let file = File::open("tests/small_segmasker_output.txt").unwrap();
    let mut reader = BufReader::new(file);
    let results = segmasker::parse(&mut reader);

    let spac1250_07 = results.get("SPAC1250.07").unwrap();
    assert_eq!(spac1250_07.len(), 2);
    let second_loc = spac1250_07.get(1).unwrap();
    assert_eq!(second_loc.start, 139);
    assert_eq!(second_loc.end, 155);
}

#[test]
fn test_segmasker_location_merge() {
    let mut locations = vec![
        Location {
            start: 10,
            end: 20
        },
        Location {
            start: 30,
            end: 40
        },
        Location {
            start: 41,
            end: 50
        },
        Location {
            start: 60,
            end: 70
        },
        Location {
            start: 71,
            end: 80
        },
        Location {
            start: 75,
            end: 76
        },
        Location {
            start: 90,
            end: 100
        },
        Location {
            start: 95,
            end: 105
        },
        Location {
            start: 100,
            end: 110
        },
    ];
 
    merge_locations(&mut locations);

    assert_eq!(locations.len(), 4);
}