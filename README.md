Run with:

/var/pomcur/bin/pombase-interpro -p "postgres://<username>:<password>@localhost/<dbname>" \
    -i <(gzip -d < match_complete.xml.gz) -o pombe_interpro_results.json

in: /var/pomcur/sources/interpro/
