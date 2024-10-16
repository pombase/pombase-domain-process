extern crate domain_process;

use domain_process::interpro_parse::parse;

#[test]
fn test_parse() {

    let (interproscan_version, matches) = parse("tests/small_matches.json");

    assert_eq!(interproscan_version, "5.69-101.0");

    let spac13g6_15c = matches.get("SPAC13G6.15c").unwrap();
    assert_eq!(spac13g6_15c.gene_uniquename, "SPAC13G6.15c");

    let mobidb_match = spac13g6_15c.interpro_matches.get(0).unwrap();
    assert_eq!(mobidb_match.dbname, "MOBIDB_LITE");
    assert_eq!(mobidb_match.id, "mobidb-lite");

    let mobidb_locations = &mobidb_match.locations;
    assert_eq!(mobidb_locations.get(0).unwrap().end, 163);

    let panther_match = spac13g6_15c.interpro_matches.get(1).unwrap();
    assert_eq!(panther_match.dbname, "PANTHER");
    assert_eq!(panther_match.id, "PTHR10300");

    assert_eq!(panther_match.locations[0].end, 153);
}
