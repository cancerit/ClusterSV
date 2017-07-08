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

    Rscript clustering_index.R <input.bedpe>
    Rscript clustering_index.R <input.bedpe> <n_threads>

Options:

    input.bedpe - input BEDPE file
    n_threads - optional INT for number of threads to use. (default: 1)


Output
======
Two files are output:

* `<input_prefix>.sv_clusters_and_footprints` - this file contains the first
  ten columns in the input BEDPE file with the following columns appended.
  * Cluster ID of the cluster of each rearrangement junction
  * Number of SVs in the current SV cluster
  * Footprint ID of the low end breakpoint
  * Footprint ID of the high end breakpoint
  * Start and end coordinates of the footprint of the low end breakpoint
  * Start and end coordinates of the footprint of the high end breakpoint
  * P-value of the last agglomeration step in this cluster
* `<input_prefix>.sv_distance_pvals` - distance P-values between each pair of
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


Issues
======
Currently this script is hard-coded to only work with human reference version
hg19. It might be possible to adapt to GRCh38 by changing the chromosome length
and telomere/centromere coordinates files. 

I will try to refactor the code to accept custom chromosome size or telomere
and centromere coordinate files as inputs to the code. 


Citation
========
Coming soon hopefully...


LICENCE
========
Copyright (c) 2014-2017 Genome Research Ltd.

Author: Yilong Li <yl3@sanger.ac.uk>, <liyilong623@gmail.com>

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
