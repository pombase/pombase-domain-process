# [PomBase](/pombase) code for processing domains

This program processes the `match_complete.xml.gz` from InterPro
and also runs [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)
to generate a JSON of domain information.

The latest InterPro file is available from: https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/

UniProt IDs for pombe proteins are queried from PostgreSQL.  Those IDs are
used to filter the InterPro file.

Protein sequences are queried from PostgreSQL and are passed to TMHMM.
We run TMHMM in a separate thread while the InterPro XML is parsed and
processed.

## Running

Run with:

    PATH=$PATH_TO_TMHMM_EXE:$PATH /var/pomcur/bin/pombase-interpro \
        -p "postgres://<username>:<password>@localhost/<dbname>" \
        -i <(gzip -d < match_complete.xml.gz) -o pombe_domain_results.json


## Status

![Tests](https://github.com/pombase/pombase-domain-process/workflows/Tests/badge.svg)
