SELECT *
FROM chr_db.locus l, chr_db.cds c, chr_db.accession a
WHERE a.locus_id=l.id
AND c.locus_id=l.id
AND a.accession_num ='AB001091';


SELECT *
FROM chr_db.locus l, chr_db.cds c, chr_db.accession a
WHERE a.locus_id=l.id
AND c.locus_id=l.id
AND l.chr_location ='16q24'
AND a.latest_version='T';


SELECT *
FROM chr_db.locus l, chr_db.cds c, chr_db.accession a
WHERE a.locus_id=l.id
AND c.locus_id=l.id
AND a.latest_version='T'
AND c.product_name='H-cadherin';