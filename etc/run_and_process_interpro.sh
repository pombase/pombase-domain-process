#!/bin/sh -

# Usage:
#   cd <interproscan_directory>
#   run_and_process_interpro.sh <interproscan_version>

if [ "$1" = "" ]
then
    echo "Provide an InterProScan version as argument" 1>&2
    exit 1
fi

INTERPRO_SCAN_VERSION=$1

curl https://curation.pombase.org/dumps/latest_build/fasta/feature_sequences/peptide.fa.gz |
    gzip -d | perl -pne 's/\*$//' | perl -pne 's/\*/X/g' > pombe_peptide.fa

curl https://www.japonicusdb.org/data/genome_sequence_and_features/feature_sequences/peptide.fa.gz |
    gzip -d | perl -pne 's/\*$//' | perl -pne 's/\*/X/g' > japonicus_peptide.fa

nextflow run ebi-pf-team/interproscan6 -r $INTERPROSCAN_VERSION -profile docker --datadir data --interpro latest --input pombe_peptide.fa --max-workers 1 -c /data/pombase/interproscan6/licensed.conf

# Run only DeepTMHMM:
#nextflow run ebi-pf-team/interproscan6 -r $INTERPROSCAN_VERSION -profile docker --datadir data --interpro latest --input pombe_peptide.fa --max-workers 1 -c /data/pombase/interproscan6/licensed.conf --applications deeptmhmm

nextflow run ebi-pf-team/interproscan6 -r $INTERPROSCAN_VERSION -profile docker --datadir data --interpro latest --input japonicus_peptide.fa --max-workers 1 -c /data/pombase/interproscan6/licensed.conf

PATH=/usr/local/tmhmm-2.0c/bin:$PATH nice -19 /var/pomcur/bin/pombase-domain-process \
   --run-tmhmm -p pombe_peptide.fa \
   --extra-input-file /data/pombase/interproscan6/deeptmhmm_pombe_peptide.fa.json \
   -i pombe_peptide.fa.json -o pombe_domain_results.json

PATH=/usr/local/tmhmm-2.0c/bin:$PATH nice -19 /var/pomcur/bin/pombase-domain-process \
   --run-tmhmm -p japonicus_peptide.fa \
   -i japonicus_peptide.fa.json -o japonicus_domain_results.json
