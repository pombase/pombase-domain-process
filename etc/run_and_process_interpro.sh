#!/bin/sh -

# Usage:
#   cd <interproscan_directory>
#   run_and_process_interpro.sh

python3 setup.py -f interproscan.properties

perl -pne 's/(precalculated.match.lookup.service.proxy.(host|port)=).*/$1/' interproscan.properties > interproscan.properties.new &&
    mv interproscan.properties.new interproscan.properties

curl https://curation.pombase.org/dumps/latest_build/fasta/feature_sequences/peptide.fa.gz |
    gzip -d | perl -pne 's/\*$//' | perl -pne 's/\*/X/g' > pombe_peptide.fa

PATH=/usr/local/jdk-14.0.1/bin:$PATH nice -19 ./interproscan.sh -i pombe_peptide.fa -f json

PATH=/usr/local/tmhmm-2.0c/bin:$PATH nice -19 /var/pomcur/bin/pombase-domain-process -p pombe_peptide.fa \
   -i pombe_peptide.fa.json -o pombe_domain_results.json
