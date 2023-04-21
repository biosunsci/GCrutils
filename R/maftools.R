#' Extract cosmic signatures from Maf object
#'
#' @param laml maftools::maf object
#' @param cosmic_ver c('v2','v3')
#' @param rank number of signatures extracted
#' @param nrun resample times
#' @param single_core [MutationalPatterns::extract_signatures()]
#' @param ... pass to
#'
#' @return
#' @export
#'
#' @examples
ydo_extract_cosmic_signatures_related_all = function(laml,cosmic_ver = 'v2',rank=3,nrun=3,single_core=TRUE,...){
    "require MutationalPatterns::extract_signatures, yget_known_cancer_signatures, MutationalPatterns::fit_to_signatures"
    # get mut_mat of 96 tri-nucleotide
    laml.tnm = trinucleotideMatrix(maf = laml
                                   #, prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = laml.tnm$nmf_matrix %>% t %>% as.matrix
    # get custom signatures
    custom_sigs = MutationalPatterns::extract_signatures(mut_mat, rank = rank, nrun = nrun, single_core = single_core,...)
    # get contributions
    cancer_signatures = yget_known_cancer_signatures(cosmic_ver=cosmic_ver) %>% data.matrix
    fit_res <- fit_to_signatures(data.matrix(mut_mat), cancer_signatures)
    list(mut_mat=mut_mat,extracted_n_signatures_suit=custom_sigs,contribution=fit_res$contribution,tnm=laml.tnm)
}
