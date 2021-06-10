#! bin/bash 

Rscript R/run_cluster_sv.R \
-bedpe test_clusterSV/sample_data/gridss_SV/BLGSP-71-32-00701-01A-01E--BLGSP-71-32-00701-14A-01E--matched.bedpe \
-chr references/hg38.chrom_sizes \
-cen_telo references/hg38_centromere_and_telomere_coords.txt \
-out test_clusterSV/results/
