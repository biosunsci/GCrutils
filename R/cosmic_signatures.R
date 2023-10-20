
#' Get pre-defined Cosmic Signatures by version
#'
#' @param cosmic_ver {v2,v3}
#'
#' @return matrix, rows are 96 tri-nucleotides mutation types, columns are signatures
#'
#' @examples
yget_known_ref_cancer_signatures = function(cosmic_ver='v2'){
    # tnn_order = c('A[C>A]A','A[C>A]C','A[C>A]G','A[C>A]T','C[C>A]A','C[C>A]C','C[C>A]G','C[C>A]T','G[C>A]A','G[C>A]C','G[C>A]G','G[C>A]T','T[C>A]A','T[C>A]C','T[C>A]G','T[C>A]T','A[C>G]A','A[C>G]C','A[C>G]G','A[C>G]T','C[C>G]A','C[C>G]C','C[C>G]G','C[C>G]T','G[C>G]A','G[C>G]C','G[C>G]G','G[C>G]T','T[C>G]A','T[C>G]C','T[C>G]G','T[C>G]T','A[C>T]A','A[C>T]C','A[C>T]G','A[C>T]T','C[C>T]A','C[C>T]C','C[C>T]G','C[C>T]T','G[C>T]A','G[C>T]C','G[C>T]G','G[C>T]T','T[C>T]A','T[C>T]C','T[C>T]G','T[C>T]T','A[T>A]A','A[T>A]C','A[T>A]G','A[T>A]T','C[T>A]A','C[T>A]C','C[T>A]G','C[T>A]T','G[T>A]A','G[T>A]C','G[T>A]G','G[T>A]T','T[T>A]A','T[T>A]C','T[T>A]G','T[T>A]T','A[T>C]A','A[T>C]C','A[T>C]G','A[T>C]T','C[T>C]A','C[T>C]C','C[T>C]G','C[T>C]T','G[T>C]A','G[T>C]C','G[T>C]G','G[T>C]T','T[T>C]A','T[T>C]C','T[T>C]G','T[T>C]T','A[T>G]A','A[T>G]C','A[T>G]G','A[T>G]T','C[T>G]A','C[T>G]C','C[T>G]G','C[T>G]T','G[T>G]A','G[T>G]C','G[T>G]G','G[T>G]T','T[T>G]A','T[T>G]C','T[T>G]G','T[T>G]T')
    # if (cosmic_ver=='v2'){
    #     # get cosimic_cancer_signatrues v2 30 sigs
    #     sp_url <- system.file('extdata/signatures_probabilities.txt',package = 'GCrutils')
    #     stopifnot(file.exists(sp_url))
    #     cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
    #     # Match the order of the mutation types to MutationalPatterns standard
    #     new_order = match(tnn_order, cancer_signatures$Somatic.Mutation.Type)
    #     # Reorder cancer signatures dataframe
    #     cancer_signatures = cancer_signatures[as.vector(new_order),]
    #     # Add trinucletiode changes names as row.names
    #     row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
    #     # Keep only 96 contributions of the signatures in matrix
    #     cancer_signatures = as.matrix(cancer_signatures[,4:33])
    # }else if (cosmic_ver =='v3'){
    #     # get cosmic_caner_signatures v3 78 sigs
    #     sp_url <- system.file('extdata/COSMIC_v3.2_SBS_GRCh37.txt',package = 'GCrutils')
    #     stopifnot(file.exists(sp_url))
    #     cancer_signatures2 = read.table(sp_url, sep = "\t", header = TRUE)
    #
    #     new_order = match(tnn_order, cancer_signatures2$Type)
    #     # Reorder cancer signatures dataframe
    #     cs = cancer_signatures2[as.vector(new_order),]
    #     #         cs = cancer_signatures2
    #     # Keep only 96 contributions of the signatures in matrix
    #     row.names(cs) = tnn_order
    #     cancer_signatures = cs[,-1]
    # }else{
    #     stop("wrong COSMIC_VER value, must in ['v2','v3']")
    # }

    ## shorts
    if (cosmic_ver=='v2'){
        # get cosimic_cancer_signatrues v2 30 sigs, originated from signatures_probabilities.txt
        sp_url <- system.file('extdata/cosmic_v2.RDS',package = 'GCrutils')
        ref_signatures = readRDS(sp_url)
    }else if(cosmic_ver=='v3'){
        # get cosimic_cancer_signatrues v3 SBS 78 sigs, originated from COSMIC_v3.2_SBS_GRCh37.txt
        sp_url <- system.file('extdata/cosmic_v3_SBS.RDS',package = 'GCrutils')
        ref_signatures = readRDS(sp_url)
    }else{
        stop("wrong COSMIC_VER value, must in ['v2','v3']")
    }
    ref_signatures
}


#' extract tnm, signatures and signatures contributions using maf objects
#'
#' @details use maftools::trinucleotideMatrix to get tri-nucleotides mutations
#' use MutationalPatterns::extract_signatures to extract maf specific signatures
#' use MutationalPatterns::fit_to_signatures to get cosine similaries and contributions
#'
#' @param laml Maf type obj, by maftools::read.maf
#' @param cosmic_ver {v2:30 classical, v3: SBS signatures}
#' @param rank 3, numbers of signatures to extract
#' @param nrun 3, re sample times in extracting signatures
#' @param single_core
#' @param ...
#'
#' @return list, of list(mut_mat=mut_mat,extracted_n_signatures_suit=custom_sigs,contribution=fit_res$contribution,tnm=laml.tnm)
#' @export
#'
#' @examples
yget_signatures_related_all = function(laml, cosmic_ver = 'v2',rank=3, nrun=3,single_core=TRUE,...){
    "require MutationalPatterns::extract_signatures, yget_known_cancer_signatures, MutationalPatterns::fit_to_signatures"
    # get mut_mat of 96 tri-nucleotide
    laml.tnm = maftools::trinucleotideMatrix(maf = laml
                                   #, prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = laml.tnm$nmf_matrix %>% t %>% as.matrix
    # get custom signatures
    custom_sigs = MutationalPatterns::extract_signatures(mut_mat, rank = rank, nrun = nrun, single_core = single_core,...)
    # get contributions
    cancer_signatures = yget_known_ref_cancer_signatures(cosmic_ver=cosmic_ver) %>% data.matrix
    fit_res <- MutationalPatterns::fit_to_signatures(data.matrix(mut_mat), cancer_signatures)
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
        tnm = maftools::trinucleotideMatrix(maf = laml_or_tnm
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
        laml.sign = maftools::estimateSignatures(mat = tnm, nTry = nTry)
        plotCophenetic(res = laml.sign)
        return (laml.sign)
    }
    if (cosmic_ver=='v2'){
        sig_db = 'legacy'
    }else if (cosmic_ver=='v3'){
        sig_db = 'SBS'
    }
    tnm$nmf_matrix[tnm$nmf_matrix==0] = 0.01
    laml.sig = maftools::extractSignatures(mat = tnm, n = rank)
    laml.cosm = maftools::compareSignatures(nmfRes = laml.sig, sig_db = sig_db)
    if (plot ==TRUE){
        ComplexHeatmap::Heatmap(mat = laml.cosm$cosine_similarities,col=circlize::colorRamp2(breaks = c(0, 0.5,1), colors =  c("blue", "white", "red"), space = 'RGB'),
                                cluster_rows = FALSE, column_title = "cosine similarity against validated signatures")
        maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = sig_db)
    }
    # concat the two result lists
    c(laml.sig,laml.cosm)
}

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


#' Title
#'
#' @param laml.cos
#' @param laml.sig
#' @param titles
#' @param .with_ori_data
#'
#' @return
#' @export
#'
#' @examples
yplot_signatrues_extracted_gg = function (laml.cos, laml.sig, titles = "", .with_ori_data=FALSE) {

    #     stopifnot('comparedSIGs' %in% class(laml.cos))
    #     stopifnot('extractedSIGs' %in% class(laml.sig))

    COL_TITV = structure(c('#0169B0','#8BC3E0','#F7DAC6','#ECA382','#C92027','#CAB3D7')
                         , name =c('C>A','C>G','C>T','T>A','T>C','T>G'))

    res = list()

    res$heatmap = pheatmap(laml.cos$cosine_similarities,name = 'Similarity') %>% ggplotify::as.ggplot()

    res$byExtractedSig = list()
    n_sigs = laml.cos$cosine_similarities %>% rownames
    mat = laml.sig$signatures %>% data.frame %>% rownames_to_column("N96") %>%
        mutate(TN = factor(str_sub(N96, 3, 5)))

    for (i in seq_along(n_sigs)) {
        signame = n_sigs[[i]]
        index = which.max(laml.cos$cosine_similarities[i, ])
        best.match = names(index)
        cosine.similarity = laml.cos$cosine_similarities[i, index]
        tit = laml.cos$best_match[[signame]]$best_match #paste0("Best match: ", best.match, " [cosine-similarity:", cosine.similarity, "]")
        aet = laml.cos$best_match[[signame]]$aetiology #paste0("Aetiology: ", laml.cos$aetiology_db[names(index), "aetiology"])
        subtitle = paste0(tit, "\n", aet)
        gg = mat %>% ggplot(aes(x = N96, y = !!sym(signame), fill = TN)) +
            geom_col(width = 0.8) + theme_classic(base_size = 12) +
            facet_wrap(~TN, scales = "free_x", nrow = 1, strip.position = "bottom") +
            theme(axis.text.x = element_blank(), plot.title = element_text(size = 14,
                                                                           hjust = 0.5, vjust = -14), plot.subtitle = element_text(size = 14,
                                                                                                                                   hjust = 0.5, vjust = -14), strip.background.x = element_rect(fill = NA,
                                                                                                                                                                                                color = NA), strip.text = element_text(size = 16),
                  axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(), panel.spacing = unit(0,
                                                                       "pt"), legend.position = "none") + scale_fill_manual(breaks = names(COL_TITV),
                                                                                                                            values = COL_TITV) + scale_y_continuous(limits = c(0,
                                                                                                                                                                               0.3)) + ylab(signame) + ggtitle(label = "", subtitle = subtitle)
        if (titles != "") {
            gg = gg + ggtitle(label = titles)
        }
        res$byExtractedSig[[signame]] = list(gg = gg, title = titles, subtitle = subtitle,
                                             best.match = best.match, cosine.similarity = cosine.similarity)
    }
    res$barplot = ggpubr::ggarrange(plotlist = res$byExtractedSig %>% purrr::map(.f = function(i)i$gg), ncol = 1)
    if (.with_ori_data==TRUE){
        res = c(res,laml.sig,laml.cos)
    }
    res
}


#' Title
#'
#' @param res.maftools.sigs
#' @param sel_sig
#' @param sel_anno
#' @param title
#'
#' @return
#' @export
#'
#' @examples
yplot_tnm_cosmic = function(res.maftools.sigs, sel_sig='1',sel_anno='best',title=''){
    if (res.maftools.sigs$sig_db == 'legacy') PREFIX = 'COSMIC_' else
        if (res.maftools.sigs$sig_db == 'SBS') PREFIX = 'SBS_' else stop('wrong res.maftools.sigs obj, must contain sig_db value !')
        name_sigs = res.maftools.sigs$signatures %>% colnames
        name_cosmics = colnames(res.maftools.sigs$cosine_similarities)
        if (is.character(sel_sig)){
            if (check.numeric(sel_sig)) sig_name = paste0('Signature_',sel_sig) else sig_name = sel_sig
        }else if(is.integer(sel_sig)){
            sig_name = name_sigs[[sel_sig]]
        }else stop('wrong sel_sig value, should be: integer / character')
        stopifnot(sig_name %in% name_sigs)

        if (length(sel_anno) == 1 && sel_anno == 'best'){
            cos_name = which.max(res.maftools.sigs$cosine_similarities[sig_name,]) %>% names
            title_prefix = paste0("Best match: ", cos_name)
        }else{
            if (is.character(sel_anno)){
                if (check.numeric(sel_anno)){
                    cos_name = paste0(PREFIX,sel_anno)
                }else{
                    cos_name = sel_anno
                }
            }else if(is.integer(sel_anno)){
                cos_name = name_cosmics[sel_anno]
            }else stop('wrong sel_anno value, should be: integer / character')
            title_prefix = cos_name
        }
        stopifnot(cos_name %in% name_cosmics)

        cosine.similarity = res.maftools.sigs$cosine_similarities[sig_name, cos_name]
        tit = paste0(title_prefix, " [cosine-similarity:", cosine.similarity, "]")
        aet = paste0("Aetiology: ", res.maftools.sigs$aetiology_db[cos_name, "aetiology"])
        subtitle = paste0(tit, "\n", aet)

        sig_name_sym = ensym(sig_name)
        #names_drop = setdiff(name_sigs,sig_name)
        mat2 = res.maftools.sigs$signatures %>% data.frame %>% rownames_to_column("N96") %>%
            mutate(TN = factor(str_sub(N96, 3, 5))) #%>% select(!(!!names_drop))

        gg = mat2 %>% ggplot(aes(x = N96, y = !!sig_name_sym, fill = TN)) +
            geom_col(width = 0.8) + theme_classic(base_size = 12) +
            facet_wrap(~TN, scales = "free_x", nrow = 1, strip.position = "bottom") +
            theme(axis.text.x = element_blank(), plot.title = element_text(size = 14,
                                                                           hjust = 0.5, vjust = -14), plot.subtitle = element_text(size = 14,
                                                                                                                                   hjust = 0.5, vjust = -14), strip.background.x = element_rect(fill = NA,
                                                                                                                                                                                                color = NA), strip.text = element_text(size = 14),
                  axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(), panel.spacing = unit(0,
                                                                       "pt"), legend.position = "none") + scale_fill_manual(breaks = names(COL_TITV),
                                                                                                                            values = COL_TITV) + scale_y_continuous(limits = c(0,
                                                                                                                                                                               0.3)) + ylab(NULL) + ggtitle(label = "", subtitle = subtitle) +
            theme(plot.margin = unit(c(-10,2,-5,2),units = 'pt'))
        if (title != "") {
            gg = gg + ggtitle(label = title)
        }
        make.custom(6.5,3)
        gg
}

##@

#' Title
#'
#' @param cosmic_ver
#'
#' @return
#' @export
#'
#' @examples
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

#' Extract cosmic signatures from Maf object
#'
#' @description
#' using a maftools::readMaf returned MAF object
#' sequentially use maftools::trinucleotideMatrix ->
#' MutationalPatterns::extract_signatures ->
#' MutationalPatterns::fit_to_signatures to generate cosmic signatures
#' cosmic_ver, v2 means 30 v2 COSMIC signatures, v3 means the 79 v3 SBS signatures
#' rank, nrun and single_corepassed to MutationalPatterns::extract_signatures
#'
#'
#'
#' @param laml maftools::maf object
#' @param cosmic_ver c('v2','v3')
#' @param rank number of signatures extracted
#' @param nrun resample times
#' @param single_core [MutationalPatterns::extract_signatures()]
#' @param ... pass to MutationalPatterns::extract_signatures
#'
#' @return list contains:
#'  mut_mat: the tri-nucleotide-matrix:TNM
#'  extracted_n_signatures_suit: extracted signatures
#'  contribution: matrix of each samples' known COSMIC signatures' contribution components
#'  tnm: list of maftools::trinucleotideMatrix return values
#'
#' @export
#'
#' @examples
yget_signatures_related_all = function(laml,cosmic_ver = 'v2',rank=3,nrun=3,single_core=TRUE,...){
    "require MutationalPatterns::extract_signatures, yget_known_cancer_signatures, MutationalPatterns::fit_to_signatures"
    # get mut_mat of 96 tri-nucleotide
    laml.tnm = maftools::trinucleotideMatrix(maf = laml
                                   #, prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = laml.tnm$nmf_matrix %>% t %>% as.matrix
    # get custom signatures
    custom_sigs = MutationalPatterns::extract_signatures(mut_mat, rank = rank, nrun = nrun, single_core = single_core,...)
    # get contributions
    cancer_signatures = yget_known_cancer_signatures(cosmic_ver=cosmic_ver) %>% data.matrix
    fit_res <- MutationalPatterns::fit_to_signatures(data.matrix(mut_mat), cancer_signatures)
    list(mut_mat=mut_mat,extracted_n_signatures_suit=custom_sigs,contribution=fit_res$contribution,tnm=laml.tnm)
}


