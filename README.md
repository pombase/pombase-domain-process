# PomBase code for processing domains from InterProScan

This program processes the JSON format output of `InterProScan` to
generate an simplified JSON file of domain information.

It also runs [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)
and `segmasker` and includes the results in the JSON output.

## Running

Run with:

    PATH=$PATH_TO_TMHMM_EXE:$PATH /var/pomcur/bin/pombase-domain-process \
        -p pombe_peptide.fa -i interproscan_output.json
        -o pombe_domain_results.json

## Status

![Tests](https://github.com/pombase/pombase-domain-process/workflows/Tests/badge.svg)
