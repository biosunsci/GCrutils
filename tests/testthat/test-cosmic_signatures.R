test_that("multiplication works", {
    laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
    #clinical information containing survival information and histology. This is optional
    laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
    laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

    sig = yget_signatures_related_all(laml,rank = 3)

    expect_true(TRUE)

})
