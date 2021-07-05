ClusterSV
=========

This is the [R](https://cran.r-project.org/) script used in the PCAWG-6 SV
mechanisms project to group rearrangements into rearrangement clusters and
footprints. 


Requirements
============

* [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)


Usage
=====
This program takes a [BEDPE format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
file as input and returns a BEPDE file with additional columns appended.
Ideally the breakpoints of each SV should be at basepair resolution, i.e. the
span or low and high end should be a single base. However, the script should
work even when the breakpoints are fuzzy. 

Note that the program will crash if the input BEDPE file is empty.

Usage:
    Rscript run_cluster_sv.R -bedpe <input.bedpe>
    Rscript run_cluster_sv.R -bedpe <input.bedpe> -chr <chrom.sizes> -cen_telo <centromere_telomere_coords> -out <output.prefix> -n <n.threads>

Named arguments:
    input.bedpe (-bedpe) - Input BEDPE file, without header line.
    chrom.sizes (-chr) - tab-delimited file with list of chromosome sizes. 
        Default: `./hg38.chrom_sizes`
    centromere_telomere_coords (-cen_telo) - tab-delimited file with list of chromosome sizes. 
        Default: `./hg38_centromere_and_telomere_coords.txt`
    output.prefix (-out) - Prefic for output files. 
        Default: Input prefix
    n_threads (-n) - Optional INT for number of threads to use. 
        Default: 1


Output
======
Two files are output:

* `<output_prefix>.sv_clusters_and_footprints` - this file contains the first
  ten columns in the input BEDPE file with the following columns appended.
  * Cluster ID of the cluster of each rearrangement junction
  * Number of SVs in the current SV cluster
  * Footprint ID of the low end breakpoint
  * Footprint ID of the high end breakpoint
  * Start and end coordinates of the footprint of the low end breakpoint
  * Start and end coordinates of the footprint of the high end breakpoint
  * P-value of the last agglomeration step in this cluster
* `<output_prefix>.sv_distance_pvals` - distance P-values between each pair of
  SVs. Note that this is the raw P-value matrix, where columns correspond to
  anchor SVs (see the PCAWG-6 SV mechanisms manuscript for details). Inside the
  algorithm, the pairwise distances used for hierarchical clustering of SVs are
  derived from (M + t(M))/2, where m is the P-value matrix in this output file
  and t() is the transpose function. 


Notes
=====
The program scales quadratically with respect to the number of SVs in the,
since the distance between each pair of SVs are computed. Thus, with

SVs with at least one breakpoint within the provided telomere and centromere
boundaries are discarded. Therefore the output may have fewer SVs than the
input. This feature was originally created to ensure that the output from this
script are compatible with downstream analysis scripts in the PCAWG-6 SV
mechanisms project...

The algorithm behind this script is described in the PCAWG-6 SV mechanisms
manuscript. 


Citation
========
Coming soon hopefully...


LICENSE
========
Copyright (c) 2014-2017 Genome Research Ltd.

Author: Yilong Li <yl3.at.sanger.ac.uk>, <liyilong623.at.gmail.com>
Contributor: Helena Winata <hwinata@mednet.ucla.edu>

ClusterSV is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
