# USAGE: Rscript main.R -bedpe <input.bedpe> -chr <chrom.sizes> -cen_telo <centromere_telomere_coords> -out <output.dir> -n <n.threads>

library(gtools)
library(doParallel)
library(R.utils)

# Set up command line params
required_params <- list('chr', 'cen_telo', 'out', 'bedpe', 'n')

default_params <- list(
    chr = 'hg38.chrom_sizes',
    cen_telo = 'hg38_centromere_and_telomere_coords.txt',

    n = 1
    );

ARGS <- commandArgs(
    trailingOnly = FALSE,
    asValues = TRUE,
    defaults = default_params
    );

# load functions
source(file.path(dirname(ARGS$file),'clustering_index.R'));

####################################################################################################
#     C L U S T E R - S V
####################################################################################################

if (all(names(default_params) %in% required_params)) {
    chr_sizes <- get_chr_sizes(ARGS$chr);
    genome_size <- sum(as.numeric(chr_sizes));
    centromere_telomere_coords <- get_centromere_telomere_coords(ARGS$cen_telo);
    registerDoParallel(cores = ARGS$n);

    file <- ARGS$bedpe;
    if (!file.exists(file)) {
        stop(sprintf("Input file '%s' does not exist!", file))
        };

    if (is.null(ARGS$out)) {
        output_dir <- sub('.bedpe$', '', file)
    } else {
        output_dir <- ARGS$out
    }
} else {
    cat('No input file provided.\n');
    cat('Usage:\n');
    cat('    run_cluster_sv.R -bedpe <input.bedpe> [ARGS] \n');
    cat('    run_cluster_sv.R -bedpe <input.bedpe> -chr <chrom.sizes> -cen_telo <centromere_telomere_coords> -out <output.dir> -n <n_threads>
\n');
    cat('\n');
    cat('Arguments:\n');
    cat('    input.bedpe (-bedpe) - Input BEDPE file, without header line.\n');
    cat('    chrom.sizes (-chr) - tab-delimited file with list of chromosome sizes \n');
    cat('    centromere_telomere_coords (-cen_telo) - tab-delimited file with list of chromosome sizes \n');
    cat('    outpit.dir (-out) - Optional output directory (default: input file directory) \n');
    cat('    n_threads (-n) - Optional INT for number of threads to use. (default: 1)\n')
};



if (file.info(file)$size > 0) {
    d = read.table(file, header = F, sep = "\t", stringsAsFactors = F, fill = T, comment = "")
} else {
    warning("Input file is empty.\n")
    out_table = data.frame(NULL)
    out_mat = matrix(NA, nrow = 1)
};

res = compute_clusters_and_footprints(d);
out_mat = res$m;
clustering_fdr_cutoff = res$clustering_fdr_cutoff;
by_ct_last_fdr = sapply(
    1:max(res$clust_and_fps$clust),
    function(idx) {
        if (sum(res$clust_and_fps$clust == idx) == 1) {
            NA
        } else {
            temp = out_mat[res$clust_and_fps$clust == idx, res$clust_and_fps$clust == idx]
            max(c(temp)[c(temp) < clustering_fdr_cutoff], na.rm = T)
        }
    }
);

clustering_fdr = by_ct_last_fdr[res$clust_and_fps$clust];
out_table = data.frame(
    res$clust_and_fps,
    clustering_fdr
);

# Write clusters and footprints info the the clustering matrix
write.table(
    x = out_table,
    file = file.path(output_dir, 'sv_clusters_and_footprints.tsv'),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
);
write.table(
    x = as.data.frame(out_mat),
    file = file.path(output_dir, 'sv_distance_pvals.tsv'),
    row.names = F,
    col.names = F,
    quote = F,
    sep = "\t"
);

