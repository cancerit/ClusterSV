# Functions implementing clustering index which roughly estimates how unlikely
# it is for two rearrangements to have occurred within the observed vicinity
# from each other. 
#
# Usage: Rscript clustering_index.R
# 

source('header.R')

# The four P-value computation functions below:
# 1. Inter-chr vs. inter-chr
# 2. Inter-chr vs. intra-chr
# 3. Intra-chr vs. inter-chr
# 4. Intra-chr vs. intra-chr

L = 1.5e9
p = function(dist) {
    if (dist < 0) {
        stop("Negative input for p()")
    }
    return(pmax(dist, 1)/L)
}

inv_p = function(p) {
    return(pmax(p*L, 1))
}

p_bkpts = function(bkpt_1, bkpt_2) {
    p(abs(bkpt_1 - bkpt_2))
}

bkpt_proximity_p_value = function(index_bkpt, chr_size, score_cutoff) {
    # Returns the fraction of the coordinates in the current chromosome
    # the produces a distance likelihood smaller than score_cutoff.

    distance_cutoff = inv_p(score_cutoff)
    min_coord = max(1, index_bkpt - distance_cutoff)
    max_coord = min(chr_size, index_bkpt + distance_cutoff)
    if (min_coord >= max_coord) {
        # print(min_coord)
        # print(max_coord)
        stop(sprintf("Failed with index_bkpt=%s; chr_size=%s; score_cutoff=%s", index_bkpt, chr_size, score_cutoff))
    } else {
        return((max_coord-min_coord)/chr_size)
    }
}

rg_paired_dist = function(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2) {
    p_l = ifelse(
        chr_l_1 == chr_l_2,
        p_bkpts(pos_l_1, pos_l_2),
        1
    )
    p_h = ifelse(
        chr_h_1 == chr_h_2,
        p_bkpts(pos_h_1, pos_h_2),
        1
    )
    p_l * p_h
}

rg_dist = function(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2) {
    pmin(
        rg_paired_dist(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2),
        rg_paired_dist(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_h_2, pos_h_2, chr_l_2, pos_l_2)
    )
}

# Four different functions for computing four different types of P-values:
# 1. inter-chr vs. inter-chr
# 2. inter-chr vs. intra-chr
# 3. intra-chr vs. inter-chr
# 4. intra-chr vs. intra-chr

inter_vs_inter_p_value = function(index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, chr_sizes) {
    p_match_at_low_chr   = 2 * (chr_sizes[[index_chr_l]]/genome_size) * (1 - chr_sizes[[index_chr_h]]/(genome_size-chr_sizes[[index_chr_l]])) * bkpt_proximity_p_value(index_bkpt_l, chr_sizes[[index_chr_l]], score_cutoff)
    p_match_at_high_chr  = 2 * (chr_sizes[[index_chr_h]]/genome_size) * (1 - chr_sizes[[index_chr_l]]/(genome_size-chr_sizes[[index_chr_h]])) * bkpt_proximity_p_value(index_bkpt_h, chr_sizes[[index_chr_h]], score_cutoff)

    # Calculating the area where two inter-chromosomal events produce a p smaller than threshold. 
    # 1: d1 > 0, d2 > 0. Calculate where d1 is when d2 = chr_len - index_bkpt_h
    max_d2 = chr_sizes[[index_chr_h]] - index_bkpt_h
    a = score_cutoff / max_d2 * L^2
    b = chr_sizes[[index_chr_l]] - index_bkpt_l
    area_sum = pmin(b, a) * max_d2
    if (a < b) { 
        area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
    }

    # 2: d1 > 0, d2 < 0
    max_d2 = index_bkpt_h - 1
    a = score_cutoff / max_d2 * L^2
    b = chr_sizes[[index_chr_l]] - index_bkpt_l
    area_sum = area_sum + pmin(b, a) * max_d2
    if (a < b) {
        area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
    }

    # 3: d1 < 0, d2 > 0
    max_d2 = chr_sizes[[index_chr_h]] - index_bkpt_h
    a = score_cutoff / max_d2 * L^2
    b = index_bkpt_l - 1
    area_sum = area_sum + pmin(a, b) * max_d2
    if (a < b) { 
        area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
    }

    # 2: d1 < 0, d2 < 0
    max_d2 = index_bkpt_h - 1
    a = score_cutoff / max_d2 * L^2
    b = index_bkpt_l - 1
    area_sum = area_sum + pmin(a, b) * max_d2
    if (a < b) {
        area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
    }

    p_match_at_both_chrs = 2 * (chr_sizes[[index_chr_l]]/genome_size) * (chr_sizes[[index_chr_h]]/(genome_size-chr_sizes[[index_chr_l]])) * (area_sum/chr_sizes[[index_chr_l]]/chr_sizes[[index_chr_h]])

    out = p_match_at_low_chr + p_match_at_high_chr + p_match_at_both_chrs

    if (out > 1 || out <= 0) {
        stop(sprintf("Failed at inter_vs_inter_p_value() with index_chr_l=%s index_bkpt_l=%s index_chr_h=%s index_bkpt_h=%s score_cutoff=%s chr_sizes=c(%s)", index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, paste(chr_sizes, collapse = ",")))
    }
    return(out)
}

inter_vs_intra_p_value = function(index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, chr_sizes) {
    # Index rearrangement is intra-chromosomal, other rearrangement is inter-chromosomal
    chr_size = chr_sizes[[index_chr]]

    max_dist = inv_p(score_cutoff)
    min_pos_l = pmax(1, index_bkpt_l - max_dist)
    max_pos_l = pmin(chr_size, index_bkpt_l + max_dist)
    min_pos_h = pmax(1, index_bkpt_h - max_dist)
    max_pos_h = pmin(chr_size, index_bkpt_h + max_dist)

    prob_one_end_in_index_chr = 2 * (chr_size/genome_size)
    if (min_pos_l <= max_pos_h && min_pos_h <= max_pos_l) {
        out = prob_one_end_in_index_chr * (max(max_pos_l, max_pos_h) - min(min_pos_l, min_pos_h) + 1) / chr_size
    }
    else {
        out = prob_one_end_in_index_chr * (max_pos_l-min_pos_l+1 + max_pos_h-min_pos_h+1) / chr_size
    }
    
    if (out > 1 || out <= 0) {
        stop(sprintf("Failure at inter_vs_intra_p_value() with index_chr=%s, index_bkpt_l=%s, index_bkpt_h=%s, score_cutoff=%s, chr_sizes=%s, OUT=%s", index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, "chr_sizes", out))
    }
    return(out)
}

# This is helper function for intra_vs_inter_p_value()
intra_vs_single_bkpt_p_value = function (index_bkpt, rg_sizes, chr_size, score_cutoff) {
    # Range of positions giving sufficient proximity
    max_dist = inv_p(score_cutoff)

    if (all(rg_sizes > chr_size)) {
        rg_sizes = chr_size
    }

    out = mean(sapply(
        rg_sizes[rg_sizes <= chr_size],
        function(size) {
            min_pos_l = pmin(pmax(1, index_bkpt - max_dist), chr_size - size)
            max_pos_l = pmin(chr_size - size, index_bkpt + max_dist)
            min_pos_h = pmax(1, index_bkpt - size - max_dist)
            max_pos_h = pmax(1, pmin(chr_size - size, index_bkpt - size + max_dist))
            # cat(paste(size, min_pos_l, max_pos_l, min_pos_h, max_pos_h, "\n"))

            # Sanity check
            if (min_pos_l > max_pos_l || min_pos_h > max_pos_h) {
                stop()
            }

            if (
                min_pos_l <= max_pos_h &&  # Low end range start vs. high end range end
                max_pos_l >= min_pos_h     # Low end range end vs. high end range start
            ) {
                (pmax(max_pos_l, max_pos_h) - pmin(min_pos_l, min_pos_l) + 1)/chr_size
            }
            else {
                (max_pos_l-min_pos_l+1 + max_pos_h-min_pos_h+1) / chr_size
            }
        }
    ))

    if (out > 1 || out < 0) {
        stop(sprintf("Failed with index_bkpt=%s; chr_size=%s; score_cutoff=%s; rg_sizes=%s", index_bkpt, chr_size, score_cutoff, paste(rg_sizes, collapse = ",")))
    }
    out
}

intra_vs_inter_p_value = function(index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, chr_sizes, rg_sizes) {
    # Index rearrangement is inter-chromosomal, other rearrangement is intra-chromosomal
    out = chr_sizes[[index_chr_l]]/genome_size * intra_vs_single_bkpt_p_value(index_bkpt_l, rg_sizes, chr_sizes[[index_chr_l]], score_cutoff) +
          chr_sizes[[index_chr_h]]/genome_size * intra_vs_single_bkpt_p_value(index_bkpt_h, rg_sizes, chr_sizes[[index_chr_h]], score_cutoff)

    if (out > 1 || out < 0) {
        stop()
    }
    out
}

quadratic_root = function(a, b, c) {
    D = b^2 - 4*a*c
    if (D < 0) {
        return(c())
    }
    else if (D == 0) {
        return(-b/2*a)
    }
    else {
        return(sort(c(-b + sqrt(D), -b - sqrt(D)) / (2 * a)))
    }
}

# This is a helper function for a single secondary breakpoint
intra_vs_single_intra_p_value = function(il, ih, score_cutoff, s, chr_size) {
    il = as.numeric(il)
    ih = as.numeric(ih)
    score_cutoff = as.numeric(score_cutoff)
    s = pmax(as.numeric(s), 1)
    chr_size = as.numeric(chr_size)

    # Is there any chance this size produces good enough score?
    if (p(0) * p(abs(ih-il-s)) == score_cutoff) {
        return(1/(chr_size - s + 1))
    }
    if (p(0) * p(abs(ih-il-s)) > score_cutoff) {
        return(0)
    }

    good_regions_start = numeric()
    good_regions_end = numeric()

    # 1. ll and hh pairing
    if (il <= ih - s) {  # This is if low ends meet first
        # 1) pos <= il. Solution for (il-pos) * (ih-pos-size) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -il-ih+s, -il*s + il*ih - score_cutoff*L^2)
        roots = round(roots[roots-1 <= il])  # -1 for some rounding problem
        if (length(roots) == 0) {
            warning("At ll, il <= ih - s, 1): no roots found???")
            roots = il
        }
        else if (length(roots) == 2) {
            stop("At ll, il <= ih - s, 1): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, pmax(roots, 1))
            good_regions_end   = c(good_regions_end,   max(c(il, roots, 1)))
        }

        # 2) il < pos <= ih-s. Solution for (pos-il) * (ih-pos-size) - score_cutoff*L^2 = 0
        roots = quadratic_root(-1, ih-s+il, -il*ih + il*s - score_cutoff*L^2)
        if (length(roots) <= 1) {
            good_regions_start = c(good_regions_start, il)
            good_regions_end   = c(good_regions_end,   ih-s)
        }
        else {
            roots = round(roots)
            good_regions_start = c(good_regions_start, il,       roots[2])
            good_regions_end   = c(good_regions_end,   roots[1], ih-s)
        }

        # 3) ih-s < pos. Solution for (pos-il) * (pos+size-ih) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -ih+s-il, il*ih - il*s - score_cutoff*L^2)
        roots = round(roots[roots+1 >= ih-s])  # +1 for some rounding problems
        if (length(roots) == 0) {
            warning("At hh, il <= ih - s, 3): no roots found???")
            roots = ih - s
        }
        else if (length(roots) == 2) {
            stop("At hh, il <= ih - s, 3): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, min(c(ih-s, roots, chr_size)))
            good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
        }
    }
    else {  # This is if high ends meet first
        # 1) pos+size <= ih. Solution for (il-pos) * (ih-pos-size) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -il-ih+s, -il*s + il*ih - score_cutoff*L^2)
        roots = round(roots[roots-1 <= ih-s])  # -1 for some rounding problem
        if (length(roots) == 0) {
            warning("At ll, il > ih - s, 1): no roots found???")
            roots = ih - s
        }
        else if (length(roots) == 2) {
            stop("At ll, il > ih - s, 1): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, pmax(1, roots))
            good_regions_end   = c(good_regions_end,   max(c(ih-s, 1, roots)))
        }

        # 2) ih-size < pos <= il. Solution for (il-pos) * (pos+size-ih) - score_cutoff*L^2 = 0
        roots = quadratic_root(-1, ih-s+il, -il*ih + il*s - score_cutoff*L^2)
        if (length(roots) <= 1) {
            good_regions_start = c(good_regions_start, ih-s)
            good_regions_end   = c(good_regions_end,   il)
        }
        else {
            round(roots)
            good_regions_start = c(good_regions_start, ih-s,   roots[2])
            good_regions_end   = c(good_regions_end,   roots[1], il)
        }

        # 3) il < pos. Solution for (pos-il) * (pos+size-ih) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -ih+s-il, il*ih - il*s - score_cutoff*L^2)
        roots = round(roots[roots+1 >= il])  # +1 for some rounding problem
        if (length(roots) == 0) {
            warning("At hh, il > ih - s, 3): no roots found???")
            roots = il
        }
        else if (length(roots) == 2) {
            stop("At hh, il > ih - s, 3): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, min(c(il+1, roots, chr_size)))
            good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
        }
    }

    # Now lh and hl pairing. Index low end will always meet with second rg high end first.
    # Need to entere here only if there's any chance some of the achieved scores
    # are more extreme than score_cutoff. 
    if (p(0) * p(abs(ih-il+s)) <= score_cutoff) {
        # 1) pos+size <= il: Solution for (il-pos-size) * (ih-pos) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -ih+s-il, il*ih - ih*s - score_cutoff*L^2)
        roots = round(roots[roots-1 <= il-s])  # -1 for some rounding problem
        if (length(roots) == 0) {
            warning("At lh, 1): no roots found???")
            roots = il - s
        }
        else if (length(roots) == 2) {
            stop("At lh, 1): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, pmax(roots, 1))
            good_regions_end   = c(good_regions_end,   max(c(il-s, roots, 1)))
        }

        # 2) pos <= ih: Solution for (pos+size-il) * (ih-pos) - score_cutoff*L^2 = 0
        roots = quadratic_root(-1, ih-s+il, -il*ih + ih*s - score_cutoff*L^2)
        roots = roots[roots+s+1 >= il & roots-1 <= ih]  # +1 and -1 for som erounding problems
        if (length(roots) <= 1) {
            good_regions_start = c(good_regions_start, pmin(il-s+1, ih))
            good_regions_end   = c(good_regions_end,   ih)
        }
        else {
            roots = round(roots)
            good_regions_start = c(good_regions_start, il-s+1,   roots[2])
            good_regions_end   = c(good_regions_end,   roots[1], ih)
        }

        # 3) pos > ih: Solution for (pos+size-il) * (pos-ih) - score_cutoff*L^2 = 0
        roots = quadratic_root(1, -ih+s-il, il*ih - ih*s - score_cutoff*L^2)
        roots = round(roots[roots+1 >= ih])  # +1 for some rounding problems
        if (length(roots) == 0) {
            warning("At hl, 3): no roots found???")
            roots = ih
        }
        else if (length(roots) == 2) {
            stop("At hl, 3): found two roots???")
        }
        else {
            good_regions_start = c(good_regions_start, min(c(roots, ih+1, chr_size)))
            good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
        }
    }

    # Make sure we don't get nasty float comparison issues.
    good_regions_start = round(good_regions_start)
    good_regions_end = round(good_regions_end)

    # Remove regions that are out of bounds
    bad_idx = good_regions_start > chr_size | good_regions_end < 1
    good_regions_start = good_regions_start[!bad_idx]
    good_regions_end   = good_regions_end[!bad_idx]
    if (any(good_regions_end < good_regions_start)) {
        # print(good_regions_start)
        # print(good_regions_end)
        # print(which(good_regions_end < good_regions_start))
        stop(sprintf("Failed after removing out of bounds regions: il=%s; ih=%s; score_cutoff=%s; s=%s; chr_size=%s", il, ih, score_cutoff, s, chr_size))
    }

    # Merge regions
    o = order(good_regions_start, good_regions_end)
    good_regions_start = good_regions_start[o]
    good_regions_end   = good_regions_end[o]
    i = 1
    while(i < length(good_regions_start)) {
        if (
            good_regions_start[i]   <= good_regions_end[i+1] &&
            good_regions_start[i+1] <= good_regions_end[i]
        ) {
            good_regions_end[i] = good_regions_end[i+1]
            good_regions_start = good_regions_start[-(i+1)]
            good_regions_end = good_regions_end[-(i+1)]
        }
        else {
            i = i + 1
        }
    }

    out = sum(pmin(good_regions_end, chr_size) - pmax(1, good_regions_start) + 1) / chr_size
    if (out > 1 || out < 0) {
        # print(good_regions_start)
        # print(good_regions_end)
        # print(good_regions_end - good_regions_start + 1)
        # print(sum(good_regions_end - good_regions_start + 1))
        # print(chr_size)
        # print(sum(good_regions_end - good_regions_start + 1) > chr_size)
        stop(sprintf("Failed after summing regions: out=%s; il=%s; ih=%s; score_cutoff=%s; s=%s; chr_size=%s", out, il, ih, score_cutoff, s, chr_size))
    }

    return(out)
}

intra_vs_intra_p_value = function(index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, chr_size, rg_sizes) {
    if (all(rg_sizes > chr_size)) {
        # return(0)
        rg_sizes = chr_size
    }

    res = sapply(
        rg_sizes[rg_sizes <= chr_size],
        function(size) {
            # print(size)
            intra_vs_single_intra_p_value(index_bkpt_l, index_bkpt_h, score_cutoff, size, chr_size)
        }
    )

    out = (chr_size/genome_size) * mean(res)
    if (out > 1 || out < 0) {
        stop()
    }
    return(out)
}


#
# Returns a matrix, where each entry is the expected number of FPs if the
# observed proximity P-value of each pair of SVs was used a the P-value threshold
# for the given pair. FP values are calculated by SV pair category, i.e.
#   del-del
#   del-inter
#   del-inv
#   del-td
#   inter-inter
#   inter-inv
#   inter-td
#   inv-inv
#   inv-td
#   td-td
#
get_clust_mat = function(d, is_already_clustered = rep(F, nrow(d))) {
    # Returns a matrix of expected numbers of given SV distances under
    # random breakpoint positioning. 

    n_inter = sum(d[!is_already_clustered, 1] != d[!is_already_clustered, 4])
    is_del = d[,1] == d[,4] & d[,9] == "+" & d[,10] == "-"
    is_td  = d[,1] == d[,4] & d[,9] == "-" & d[,10] == "+"
    is_inv = d[,1] == d[,4] & d[,9] == d[,10]
    rg_type = ifelse(is_del, "del",
              ifelse(is_td,  "td",
              ifelse(is_inv, "inv", "inter")))
    
    # Fix coordinate columns and perform some sanity checks
    bkpt_l = rowMeans(d[,2:3])
    bkpt_h = rowMeans(d[,5:6])
    if (any((is_del | is_td) & bkpt_l >= bkpt_h)) {
        warning(sprintf("is_td & bkpt_l >= bkpt_h at lines c(%s)", paste(which((is_del | is_td) & bkpt_l >= bkpt_h), collapse = ", ")))
        idx = (is_del | is_td) & bkpt_l >= bkpt_h
        temp = bkpt_l[idx]
        bkpt_l[idx] = bkpt_h[idx]
        bkpt_h[idx] = temp
    }
    temp_m = pmin(bkpt_l, bkpt_h)
    temp_M = pmax(bkpt_l, bkpt_h)
    bkpt_l = ifelse(d[,1] == d[,4], temp_m, bkpt_l)
    bkpt_h = ifelse(d[,1] == d[,4], temp_M, bkpt_h)

    rg_sizes = list(
        "del" = setNames(bkpt_h-bkpt_l, d[,7])[is_del & !is_already_clustered],
        "td"  = setNames(bkpt_h-bkpt_l, d[,7])[is_td & !is_already_clustered],
        "inv" = setNames(bkpt_h-bkpt_l, d[,7])[is_inv & !is_already_clustered]
    )

    # Replaced above with the multithreaded version
    m = foreach(i = 1:nrow(d), .combine = cbind) %dopar% {
        sapply(
            1:nrow(d),
            function(j) {
                # if (j == 1) cat(paste("[", i, "] ", d[i,7], "\n", sep = ""))
                # cat(paste("", j))
                # if (j == nrow(d)) cat("\n")
                if (j == 1) cat(sprintf("%s%s", i, if (i == nrow(d)) "\n" else " "))

                score_cutoff = rg_dist( d[i,1], bkpt_l[i], d[i,4], bkpt_h[i],
                                        d[j,1], bkpt_l[j], d[j,4], bkpt_h[j])

                if (i == j) {
                    out = NA
                }
                else if (rg_type[i] == "inter" && rg_type[j] == "inter") {
                    temp_n = n_inter + sum(is_already_clustered[c(i,j)])
                    out = temp_n * (temp_n-1) / 2
                    if (score_cutoff < 1) {
                        out = out * inter_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,4], bkpt_h[i], score_cutoff, chr_sizes)
                    }
                }
                else if (rg_type[i] != "inter" && rg_type[j] == "inter") {
                    out = (length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i]) * (n_inter + is_already_clustered[j])
                    if (score_cutoff < 1) {
                        out = out * inter_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes)
                    }
                }
                else if (rg_type[i] == "inter" && rg_type[j] != "inter") {
                    out = (n_inter+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]]) + is_already_clustered[j])
                    cur_rg_sizes = rg_sizes[[rg_type[j]]]
                    if (is_already_clustered[j]) {
                        cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
                    }
                    if (score_cutoff < 1) {
                        out = out * intra_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,4], bkpt_h[i], score_cutoff, chr_sizes, cur_rg_sizes)
                    }
                }
                else if (rg_type[i] != "inter" && rg_type[j] != "inter") {
                    if (rg_type[i] == rg_type[j]) {
                        temp_n = length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i] + is_already_clustered[j]
                        out = temp_n * (temp_n-1) / 2
                        cur_rg_sizes = rg_sizes[[rg_type[j]]]
                        if (!is_already_clustered[i]) cur_rg_sizes = cur_rg_sizes[-(which(cur_rg_sizes == bkpt_h[i]-bkpt_l[i])[1])]
                        if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
                    } else {
                        out = (length(rg_sizes[[rg_type[i]]])+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]])+is_already_clustered[j])
                        cur_rg_sizes = rg_sizes[[rg_type[j]]]
                        if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
                    }
                    if (score_cutoff < 1) {
                        out = out * intra_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes[[d[i,1]]], cur_rg_sizes)
                    }
                }
                else {
                    stop()
                }

                out
            }
        )
    }
    
    return(m)
}



get_footprints = function(pos, chr, cutoff = 0.01) {
    # Compute footprints given all breakpoints of a cluster

    # Edge case
    if (length(pos) == 2 || length(pos) == length(unique(chr))) {
        return(list(
            footprint_idx = 1:length(pos),
            footprint_coords = paste(pos, pos, sep = "-"),
            footprint_bounds = paste(pos, pos, sep = "-")
        ))
    }

    # Calculate distances between breakpoints
    o = order(chr, pos)
    reverse_o = match(1:length(pos), o)
    pos = pos[o]
    chr = chr[o]
    dists = diff(pos) + 1
    dists[chr[-length(chr)] != chr[-1]] = NA
    is_new_footprint = rep(0, length(dists) + 1)
    is_new_footprint[which(is.na(dists)) + 1] = 1

    # Below q-value is computed with a simple exponential distribution model
    # lambda = 1/mean(dists)
    # q = p.adjust(pexp(dists, lambda, lower.tail = F), method = "fdr")
    # footprint_idx_offset = rep(0, length(pos))
    # footprint_idx_offset[c(F, q <= cutoff)] = 1
    # footprint_idx = rep(1, length(pos)) + cumsum(footprint_idx_offset)
    
    # Below P-values are computed using likelihood ratio test instead
    while (sum(!is.na(dists)) > 2) {
        is_max = which.max(dists)
        p0 = 1/mean(dists, na.rm = T)
        p1 = 1/mean(dists[-is_max], na.rm = T)
        l0 = sum(dexp(dists, p0, log = T), na.rm = T)
        l1 = sum(dexp(dists[-is_max], p1, log = T), na.rm = T) + dexp(dists[is_max], 1/dists[is_max], log = T)
        D = 2 * (l1 - l0)
        # pval = pchisq(D, 1, lower.tail = F) * sum(!is.na(dists))
        pval = p.adjust(pchisq(D, 1, lower.tail = F), method = "fdr", n = sum(!is.na(dists)))
        if (pval < cutoff) {
            is_new_footprint[is_max+1] = 1
            dists[is_max] = NA
        } else {
            break
        }
    }

    # After computing footprints based on all distances, perform per-chromosome
    # distance tests. 
    # print(dists)
    # print(is_new_footprint)
    dists = c(dists, NA)
    is_new_footprint = c(is_new_footprint, NA)
    for (c in unique(chr)) {
        idx = chr == c
        cur_dists = dists[idx]
        cur_is_new_footprint = is_new_footprint[idx]
        print(cur_dists)
        while (sum(!is.na(cur_dists)) > 1) {
            is_max = which.max(cur_dists)
            p0 = 1/mean(cur_dists, na.rm = T)
            p1 = 1/mean(cur_dists[-is_max], na.rm = T)
            l0 = sum(dexp(cur_dists, p0, log = T), na.rm = T)
            l1 = sum(dexp(cur_dists[-is_max], p1, log = T), na.rm = T) + dexp(cur_dists[is_max], 1/cur_dists[is_max], log = T)
            D = 2 * (l1 - l0)
            # pval = pchisq(D, 1, lower.tail = F) * sum(!is.na(cur_dists))
            pval = p.adjust(pchisq(D, 1, lower.tail = F), method = "fdr", n = sum(!is.na(cur_dists)))
            print(pval)
            if (pval < cutoff) {
                cur_is_new_footprint[is_max+1] = 1
                cur_dists[is_max] = NA
            } else {
                break
            }
        }
        dists[idx] = cur_dists
        is_new_footprint[idx] = cur_is_new_footprint
    }
    dists = dists[-length(dists)]
    is_new_footprint = is_new_footprint[-length(is_new_footprint)]
    print(is_new_footprint)
    print(dists)

    footprint_idx = rep(1, length(pos)) + cumsum(is_new_footprint)
    lambda = 1/mean(dists, na.rm = T)

    # Chromosome of each footprint
    chr_of_footprint_idx = c()
    for (k in 1:max(footprint_idx)) {
        temp = unique(chr[footprint_idx == k])
        if (length(temp) != 1) { stop("Sanity check failed") }
        chr_of_footprint_idx[k] = temp
    }

    # Some helper functions
    ptel_coord = function(chr_name) centromere_telomere_coords[chr_name, "ptel"]
    qtel_coord = function(chr_name) centromere_telomere_coords[chr_name, "qtel"]
    cen_start_coord = function(chr_name) centromere_telomere_coords[chr_name, "cen_start"]
    cen_end_coord = function(chr_name) centromere_telomere_coords[chr_name, "cen_end"]
    cen_coord = function(chr_name) centromere_telomere_coords[chr_name, "cen"]
    exp_distr_lambda_hat = function(x) 1/mean(x)

    # Footprint extends to chromosome ends or centromeres?
    footprint_coords_of_idx = c()
    footprint_bounds_of_idx = c()
    for (k in 1:max(footprint_idx)) {
        idx = which(footprint_idx == k)
        footprint_chr = chr_of_footprint_idx[k]
        footprint_coords_of_idx[k] = paste(range(pos[idx]), collapse = "-")

        if (length(idx) == 1) {  # No footprint for singleton SVs
            footprint_bounds_of_idx[k] = paste(pos[idx], pos[idx], sep = "-")
        } else {
            dists = diff(pos[idx])
            lambda = 1/mean(dists)
            is_first_footprint_of_chr = k == min(footprint_idx[chr_of_footprint_idx[footprint_idx] == footprint_chr])
            #print(paste0('FOOTPRINT', footprint_chr))
            is_before_centromere = min(pos[idx]) < cen_coord(footprint_chr)
            if (is_first_footprint_of_chr && is_before_centromere) {  # Only the first footprint can start from p-telomere
                if (pexp(pos[idx][1] - ptel_coord(footprint_chr), lambda, lower.tail = F) <= cutoff) {
                    start = pos[idx][1]
                } else {
                    start = "ptel"
                }
            } else if (is_first_footprint_of_chr || (pos[min(idx)-1] < cen_coord(footprint_chr) && !is_before_centromere)) {  # Only the footprint closest to centromere can start from centromere
                if (pexp(min(pos[idx]) - cen_end_coord(footprint_chr), lambda, lower.tail = F) > cutoff) {
                    start = "cen_end"
                } else {
                    start = min(pos[idx])
                }
            } else {
                start = min(pos[idx])
            }

            is_last_footprint_of_chr = k == max(footprint_idx[chr_of_footprint_idx[footprint_idx] == footprint_chr])
            is_after_centromere = max(pos[idx]) > cen_coord(footprint_chr)
            if (is_last_footprint_of_chr && is_after_centromere) {  # Only the last footprint can end at the q-telomere
                if (pexp(qtel_coord(footprint_chr) - pos[idx][length(idx)], lambda, lower.tail = F) <= cutoff) {
                    end = pos[idx][length(idx)]
                } else {
                    end = "qtel"
                }
            } else if (is_last_footprint_of_chr || (pos[max(idx)+1] > cen_coord(footprint_chr) && !is_after_centromere)) {  # Only the footprint closest to centromere can start from centromere
                if (pexp(cen_start_coord(footprint_chr) - max(pos[idx]), lambda, lower.tail = F) > cutoff) {
                    end = "cen_start"
                } else {
                    end = max(pos[idx])
                }
            } else {
                end = max(pos[idx])
            }

            footprint_bounds_of_idx[k] = paste(start, end, sep = "-")
        }
    }

    footprint_idx = footprint_idx[reverse_o]
    return(list(
        footprint_idx = footprint_idx,
        footprint_coords = footprint_coords_of_idx[footprint_idx],
        footprint_bounds = footprint_bounds_of_idx[footprint_idx]
    ))
}


clust_rgs = function(d) {
    cutoff_1 = 0.01
    cutoff_2 = 0.05

    m = get_clust_mat(d)
    m2 = m/2 + t(m)/2
    h = hclust(as.dist(m2), method = "single")
    fdr = h$height * 10 / seq_along(h$height)  # Times ten, because there are ten different pairs of SV types between del, td, inv and translocation
    fdr_cutoff = cutoff_1
    if (fdr[1] > fdr_cutoff) {
        ct = 1:nrow(d)
    } else if (all(fdr <= fdr_cutoff)) {
        ct = rep(1, nrow(d))
    } else {
        # hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
        hc = h$height[min(max(which(fdr <= fdr_cutoff)), length(h$height))]
        ct = cutree(h, h = hc)
    }
    is_already_clustered = duplicated(ct) | duplicated(ct, fromLast = T)

    m = get_clust_mat(d, is_already_clustered)
    m2 = m/2 + t(m)/2
    h = hclust(as.dist(m2), method = "single")
    fdr = h$height * 10 / seq_along(h$height)  # Times ten, because there are ten different pairs of SV types between del, td, inv and translocation
    fdr_cutoff = cutoff_2
    if (fdr[1] > fdr_cutoff) {
        ct = 1:nrow(d)
    } else if (all(fdr <= fdr_cutoff)) {
        ct = rep(1, nrow(d))
    } else {
        # hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
        hc = h$height[min(max(which(fdr <= fdr_cutoff)), length(h$height))]
        ct = cutree(h, h = hc)
    }

    cutree_res = ct
    list(
        hclust = h,
        fdr = fdr,
        cutree = ct,
        fdr_cutoff = fdr_cutoff,
        m = m
    )
}


compute_clusters_and_footprints = function(d) {
    d[,1] = as.character(d[,1])
    d[,4] = as.character(d[,4])
    if (nrow(d) == 1) {
        ct = 1
        m = matrix(NA, nrow = 1)
        footprint_idx = sprintf("%s.chr%s.%s", c(1, 1), c(d[1,1], d[1,4]), c(1, 2))
        footprint_bounds = paste(
            c(rowMeans(d[, 2:3]), rowMeans(d[, 5:6])),
            c(rowMeans(d[, 2:3]), rowMeans(d[, 5:6])),
            sep = "-"
        )
        clustering_fdr_cutoff = 1
    }
    else {
        res = clust_rgs(d)
        ct = res$cutree
        m = res$m
        clustering_fdr_cutoff = res$fdr_cutoff

        # Get the footprint info
        pos = c(rowMeans(d[,2:3]), rowMeans(d[,5:6]))
        chrs = c(d[,1], d[,4])
        footprint_idx = rep("NA", length(pos))
        footprint_bounds = rep("NA", length(pos))
        for (k in unique(ct)) {
            idx = c(ct, ct) == k
            print(k)
            res = get_footprints(pos[idx], chrs[idx])
            footprint_idx[idx] = sprintf("%s.chr%s.%s", c(ct, ct)[idx], chrs[idx], res$footprint_idx)
            footprint_bounds[idx] = res$footprint_bounds
        }
    }

    list(
        clust_and_fps = data.frame(
            d[, 1:10],
            clust = ct,
            clust_size = sapply(ct, function(x) sum(x == ct)),
            fp_l = footprint_idx[1:nrow(d)],
            fp_h = footprint_idx[-(1:nrow(d))],
            fp_l_ends = footprint_bounds[1:nrow(d)],
            fp_h_ends = footprint_bounds[-(1:nrow(d))]
        ),
        clustering_fdr_cutoff = clustering_fdr_cutoff,
        m = m
    )
}
