extern crate domain_process;

use domain_process::interpro_parse::parse;

#[test]
fn test_parse() {

    let result = parse("tests/small_matches.tsv");

    let spbc56f2_10c = result.get("SPBC56F2.10c").unwrap();
    assert_eq!(spbc56f2_10c.gene_uniquename, "SPBC56F2.10c");

    for m in &spbc56f2_10c.interpro_matches {
        eprintln!("{} {}", m.id, m.dbname);
    }

    let interpro_match = spbc56f2_10c.interpro_matches.get(0).unwrap();
    assert_eq!(interpro_match.dbname, "CDD");
    assert_eq!(interpro_match.id, "cd04188");

    let prositeprofiles_locations = &interpro_match.locations;
    assert_eq!(prositeprofiles_locations.get(0).unwrap().end, 287);
}
