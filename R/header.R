# Header filt to set up clustering_index

library(gtools)
library(doParallel)



get_chr_sizes= function() {
    temp = read.table("hg38.chrom_sizes", header = F, sep = "\t")
    chr_sizes = temp[,2]
    names(chr_sizes) = temp[,1]
    chr_sizes
}

chr_sizes = get_chr_sizes()
genome_size = sum(as.numeric(chr_sizes))

centromere_telomere_coords = as.matrix(
    read.table(
        "hg38_centromere_and_telomere_coords.txt",
        header = T,
        sep = "\t",
        row.names = 1
    )
)
rownames(centromere_telomere_coords) = c(1:22, "X") # paste0('chr', c(1:22, "X"))
colnames(centromere_telomere_coords) = c("ptel", "cen_start", "cen", "cen_end", "qtel")

