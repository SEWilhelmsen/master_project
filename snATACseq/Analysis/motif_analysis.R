# Motif analysis

# https://www.bioconductor.org/packages/release/bioc/vignettes/motifTestR/inst/doc/motifAnalysis.html

data("ar_er_peaks")
ar_er_peaks



sq <- seqinfo(ar_er_peaks)

# Obtainign av et of sequences for testing
test_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ar_er_peaks)
names(test_seq) <- as.character(ar_er_peaks)

# Obtainging a list of PWMs for testing 
data("ex_pfm")
names(ex_pfm)

ex_pfm$ESR1

# Searching sequences

score_thresh <- "70%"
getPwmMatches(ex_pfm$ESR1, test_seq, min_score = score_thresh)

getPwmMatches(ex_pfm$ESR1, test_seq, min_score = score_thresh, best_only = TRUE)


bm_all <- getPwmMatches(
  ex_pfm, test_seq, min_score = score_thresh, best_only = TRUE, break_ties = "all",
  mc.cores = cores
)


countPwmMatches(ex_pfm, test_seq, min_score = score_thresh, mc.cores = cores)

# Testing for positional bias

res_pos <- testMotifPos(bm_all, mc.cores = cores)
head(res_pos)

res_abs <- testMotifPos(bm_all, abs = TRUE, mc.cores = cores) 
head(res_abs)

# Viewing matches 

topMotifs <- res_pos |>
  subset(fdr < 0.05) |>
  rownames()
A <- plotMatchPos(bm_all[topMotifs], binwidth = 10, se = FALSE)
B <- plotMatchPos(
  bm_all[topMotifs], binwidth = 10, geom = "col", use_totals = TRUE
) +
  geom_smooth(se = FALSE, show.legend = FALSE) +
  facet_wrap(~name)
A + B + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")


topMotifs <- res_abs |>
  subset(fdr < 0.05) |>
  rownames()
A <- plotMatchPos(bm_all[topMotifs], abs = TRUE, type = "heatmap") +
  scale_fill_viridis_c()
B <- plotMatchPos(
  bm_all[topMotifs], abs = TRUE, type = "cdf", geom = "line", binwidth = 5
)
A + B + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")



# Defining a set of control sequences 
data("zr75_enh")
mean(overlapsAny(ar_er_peaks, zr75_enh))

ar_er_peaks$feature <- ifelse(
  overlapsAny(ar_er_peaks, zr75_enh), "enhancer", "other"
)
chr1 <- GRanges(sq)[1]
bg_ranges <- GRangesList(
  enhancer = zr75_enh,
  other = GenomicRanges::setdiff(chr1, zr75_enh)
)

data("hg19_mask")
set.seed(305)
rm_ranges <- makeRMRanges(
  splitAsList(ar_er_peaks, ar_er_peaks$feature),
  bg_ranges, exclude = hg19_mask,
  n_iter = 100
)

rm_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rm_ranges)
mcols(rm_seq) <- mcols(rm_ranges)




# Testing for enrichment
enrich_res <- testMotifEnrich(
  ex_pfm, test_seq, rm_seq, min_score = score_thresh, model = "quasi", mc.cores = cores
)
head(enrich_res)


iter_res <- testMotifEnrich(
  ex_pfm, test_seq, rm_seq, min_score = score_thresh, mc.cores = cores, model = "iteration"
)
head(iter_res)

ex_pfm |>
  getPwmMatches(test_seq, min_score = score_thresh, mc.cores = cores) |>
  lapply(\(x) x$seq) |>
  plotOverlaps(type = "upset")


# Clustered motifs
c("ANDR", "FOXA1") |>
  lapply(
    \(x) create_motif(ex_pfm[[x]], name = x, type = "PPM")
  ) |>
  view_motifs()

cl <- clusterMotifs(ex_pfm, plot = TRUE, labels = NULL)

ex_cl <- split(ex_pfm, cl)
names(ex_cl) <- vapply(split(names(cl), cl), paste, character(1), collapse = "/")

cl_matches <- getClusterMatches(
  ex_cl, test_seq, min_score = score_thresh, best_only = TRUE
)
cl_matches


testClusterPos(cl_matches, test_seq, abs = TRUE)


testClusterEnrich(
  ex_cl, test_seq, rm_seq, min_score = score_thresh, model = "quasi", mc.cores = cores
)




