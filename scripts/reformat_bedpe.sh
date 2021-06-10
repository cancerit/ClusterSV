sed '/^#/d' BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.gridss_somatic_filtered.bedpe | \
 awk -F '\t' 'BEGIN {OFS="\t"}; {for(i=0; i<=10; i++) print $i}' > BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.bedpe

sed '/^#/d' BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.gridss_somatic_filtered.bedpe | \
 awk -F '\t' 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.bedpe

Rscript clustering_index.R test_clusterSV/sample_data/gridss_SV/BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.bedpe

sed '/^#/d' /projects/rmorin/projects/gambl-repos/gambl-hwinata/results/gambl/gridss-1.0/99-outputs/bedpe/genome--hg38/BLGSP-71-32-00701-01A-01E--BLGSP-71-32-00701-14A-01E--matched.gridss_somatic_filtered.bedpe | \
 awk -F '\t' 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > BLGSP-71-32-00701-01A-01E--BLGSP-71-32-00701-14A-01E--matched.bedpe

sed '/^#/d' BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.gridss_somatic_filtered.bedpe | \
 sed -e 's/chr//g' |\
 awk -F '\t' 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' \
 > BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.bedpe


# reformat output
sed -i '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tCLUSTER_ID\tNUM_SV\tFP_ID_LOW\tFP_ID_HIGH\tFP_COORDS_LOW\tFP_COORDS_HIGH\tP_VAL' test.bedpe
