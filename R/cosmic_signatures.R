yget_known_cancer_signatures = function(cosmic_ver='v2'){
    tnn_order = c('A[C>A]A','A[C>A]C','A[C>A]G','A[C>A]T','C[C>A]A','C[C>A]C','C[C>A]G','C[C>A]T','G[C>A]A','G[C>A]C','G[C>A]G','G[C>A]T','T[C>A]A','T[C>A]C','T[C>A]G','T[C>A]T','A[C>G]A','A[C>G]C','A[C>G]G','A[C>G]T','C[C>G]A','C[C>G]C','C[C>G]G','C[C>G]T','G[C>G]A','G[C>G]C','G[C>G]G','G[C>G]T','T[C>G]A','T[C>G]C','T[C>G]G','T[C>G]T','A[C>T]A','A[C>T]C','A[C>T]G','A[C>T]T','C[C>T]A','C[C>T]C','C[C>T]G','C[C>T]T','G[C>T]A','G[C>T]C','G[C>T]G','G[C>T]T','T[C>T]A','T[C>T]C','T[C>T]G','T[C>T]T','A[T>A]A','A[T>A]C','A[T>A]G','A[T>A]T','C[T>A]A','C[T>A]C','C[T>A]G','C[T>A]T','G[T>A]A','G[T>A]C','G[T>A]G','G[T>A]T','T[T>A]A','T[T>A]C','T[T>A]G','T[T>A]T','A[T>C]A','A[T>C]C','A[T>C]G','A[T>C]T','C[T>C]A','C[T>C]C','C[T>C]G','C[T>C]T','G[T>C]A','G[T>C]C','G[T>C]G','G[T>C]T','T[T>C]A','T[T>C]C','T[T>C]G','T[T>C]T','A[T>G]A','A[T>G]C','A[T>G]G','A[T>G]T','C[T>G]A','C[T>G]C','C[T>G]G','C[T>G]T','G[T>G]A','G[T>G]C','G[T>G]G','G[T>G]T','T[T>G]A','T[T>G]C','T[T>G]G','T[T>G]T')
    if (cosmic_ver=='v2'){
        # get cosimic_cancer_signatrues v2 30 sigs
        sp_url <- '/GCI/jup/A_TSF/libs/signatures_probabilities.txt'
        cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
        # Match the order of the mutation types to MutationalPatterns standard
        new_order = match(tnn_order, cancer_signatures$Somatic.Mutation.Type)
        # Reorder cancer signatures dataframe
        cancer_signatures = cancer_signatures[as.vector(new_order),]
        # Add trinucletiode changes names as row.names
        row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
        # Keep only 96 contributions of the signatures in matrix
        cancer_signatures = as.matrix(cancer_signatures[,4:33])
    }else if (cosmic_ver =='v3'){
        # get cosmic_caner_signatures v3 78 sigs
        sp_url = '/GCI/jup/A_TSF/libs/COSMIC_v3.2_SBS_GRCh37.txt'
        cancer_signatures2 = read.table(sp_url, sep = "\t", header = TRUE)

        new_order = match(tnn_order, cancer_signatures2$Type)
        # Reorder cancer signatures dataframe
        cs = cancer_signatures2[as.vector(new_order),]
        #         cs = cancer_signatures2
        # Keep only 96 contributions of the signatures in matrix
        row.names(cs) = tnn_order
        cancer_signatures = cs[,-1]
    }else{
        stop("wrong COSMIC_VER value, must in ['v2','v3']")
    }
    cancer_signatures
}


#' Title
#'
#' @param laml
#' @param cosmic_ver
#' @param rank
#' @param nrun
#' @param single_core
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
yget_signatures_related_all = function(laml,cosmic_ver = 'v2',rank=3,nrun=3,single_core=TRUE,...){
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

#' Title
#'
#' @param laml_or_tnm
#' @param rank
#' @param cosmic_ver
#' @param plot
#' @param search_rank
#' @param nTry
#'
#' @return
#' @export
#'
#' @examples
yget_maftools_signatures = function(laml_or_tnm, rank = NA, cosmic_ver = 'v2', plot=FALSE, search_rank = FALSE, nTry =6 ){
    "extract N/rank signatures of group, the input patients are treated as a whole, and calculate the cosine similarity vs COSMIC v2/v3 signatures."
    if (is.na(rank) && search_rank == FALSE ){
        print('please run with search_rank=TRUE to find the best rank first')
        return (NA)
    }
    cls = class(laml_or_tnm)
    if (length(cls)==1 && cls == 'MAF'){
        tnm = trinucleotideMatrix(maf = laml_or_tnm
                                  #, prefix = 'chr'
                                  , add = TRUE
                                  , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

    }else if (cls == 'list' ){
        mask = names(laml_or_tnm) == c('nmf_matrix','APOBEC_scores')
        if (length(mask)==2 && all(mask)==TRUE){
            tnm = laml_or_tnm
        }
    }else{
        stop('wrong type input laml_or_tnm, must be a maftools::maf or a tnm by trinucleotideMatrix(maftools::maf)')
    }
    if (is.na(rank) && search_rank == TRUE){
        laml.sign = estimateSignatures(mat = tnm, nTry = nTry)
        plotCophenetic(res = laml.sign)
        return (laml.sign)
    }
    if (cosmic_ver=='v2'){
        sig_db = 'legacy'
    }else if (cosmic_ver=='v3'){
        sig_db = 'SBS'
    }
    tnm$nmf_matrix[tnm$nmf_matrix==0] = 0.01
    laml.sig = extractSignatures(mat = tnm, n = rank)
    laml.cosm = compareSignatures(nmfRes = laml.sig, sig_db = sig_db)
    if (plot ==TRUE){
        ComplexHeatmap::Heatmap(mat = laml.cosm$cosine_similarities,col=circlize::colorRamp2(breaks = c(0, 0.5,1), colors =  c("blue", "white", "red"), space = 'RGB'),
                                cluster_rows = FALSE, column_title = "cosine similarity against validated signatures")
        maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = sig_db)
    }
    list(laml.sig,laml.cosm)
}
