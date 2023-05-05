
test_that("multiplication works", {
    print(getwd())
    library(devtools)
    library(usethis)
    library(testthat)
    laml = yload_laml_maf(maf_ = 'wes_snv_all',gl_ = 'wes_glx',frm='export')
    laml.tnm = maftools::trinucleotideMatrix(maf = laml
                                             #, prefix = 'chr'
                                             , add = TRUE
                                             , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = laml.tnm$nmf_matrix %>% t %>% as.matrix
    # sig = yget_signatures_related_all(laml,rank = 3)
    # custom_sigs = MutationalPatterns::extract_signatures(mut_mat, rank = 3, nrun = 3, single_core = T)
    expect_true(TRUE)

})
