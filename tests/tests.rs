use std::collections::HashSet;

extern crate domain_process;

use domain_process::interpro_parse::parse;

#[test]
fn test_xml_parse() {
    let mut ids = HashSet::new();

    ids.insert("A0A001".into());

    let result = parse(&ids, "tests/small_matches.tidy.xml");

    assert_eq!(result.interpro_version, "63.0");

    let a0a001 = result.domains_by_id.get("A0A001").unwrap();
    assert_eq!(a0a001.uniprot_id, "A0A001");

    let gene3d_match = a0a001.interpro_matches.get(0).unwrap();
    assert_eq!(gene3d_match.id, "G3DSA:1.20.1560.10");
    assert_eq!(gene3d_match.dbname, "GENE3D");

    let gene3d_locations = &gene3d_match.locations;
    assert_eq!(gene3d_locations.get(0).unwrap().end, 301);
}
