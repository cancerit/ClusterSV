source('clustering_index_HW.R')
#
# THE MAIN FUNCTIONS
#
if (length(commandArgs(T)) > 0) {
    file = commandArgs(T)[1]
    body = sub(".+/", "", sub(".bedpe$", "", file))
    output_file = sprintf("%s.sv_clusters_and_footprints", file)

    if (!file.exists(file)) {
        stop(sprintf("Input file '%s' does not exist!", file))
    }
    else {
        if (file.info(file)$size > 0) {
            d = read.table(file, header = F, sep = "\t", stringsAsFactors = F, fill = T, comment = "")

            # Set up number of cores
            if (length(commandArgs(T)) > 1) {
                num_cores = as.numeric(commandArgs(T)[2])
            }
            else {
                num_cores = 1
            }
            registerDoParallel(cores = num_cores)

            res = compute_clusters_and_footprints(d)
            out_mat = res$m
            clustering_fdr_cutoff = res$clustering_fdr_cutoff
            by_ct_last_fdr = sapply(
                1:max(res$clust_and_fps$clust),
                function(idx) {
                    if (sum(res$clust_and_fps$clust == idx) == 1) {
                        NA
                    }
                    else {
                        temp = out_mat[res$clust_and_fps$clust == idx, res$clust_and_fps$clust == idx]
                        max(c(temp)[c(temp) < clustering_fdr_cutoff], na.rm = T)
                    }
                }
            )
            clustering_fdr = by_ct_last_fdr[res$clust_and_fps$clust]
            out_table = data.frame(
                res$clust_and_fps,
                clustering_fdr
            )
        }
        else {
            warning("Input file is empty.\n")
            out_table = data.frame(NULL)
            out_mat = matrix(NA, nrow = 1)
        }

        # Write clusters and footprints info the the clustering matrix
        write.table(
            out_table,
            output_file,
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
        )
        write.table(
            as.data.frame(out_mat),
            sub("sv_clusters_and_footprints$", "sv_distance_pvals", output_file),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
        )
    }
} else {
    cat("No input file provided.\n")
    cat("Usage:\n")
    cat("    clustering_index.R <input.bedpe>\n")
    cat("    clustering_index.R <input.bedpe> <n_threads>\n")
    cat("\n")
    cat("Positional arguments:\n")
    cat("    input.bedpe - Input BEDPE file - no header line.\n")
    cat("    n_threads - Optional INT for number of threads to use. (default: 1)\n")
}


