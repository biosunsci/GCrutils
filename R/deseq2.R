



#' Title use DESeq2 to do normalization and/or diff expr analysis for count data
#'
#' @param count_file_ you either provide @count_file_ which refer to a file path which will be read by `yload_dfx` or @cnt_ a data.frame/matrix obj containing rna count data with rows are genes and columns are patients
#' @param group_file_ you either provide @group_file_ which refer to a file path which will be read by `yload_dfx` or @glx_ which is a data.frame/matrix obj contaning group design information
#' @param cnt_ see @count_file_
#' @param glx_ see @group_file_
#' @param frm will be pass to `yload_dfx` if you pass @count_file_ / @group_file_
#' @param levels which levels will be compared, if NULL or NA use the first 2 levels sorted alphabetically in Clin_classification column of @glx_
#' @param export output data results, default FALSE
#' @param outputdir set the output data directory
#' @param norm default TRUE, boolean, whether to do normalization
#' @param diff default TRUE, boolean, whether to do diff expression analysis
#'
#' @importFrom tibble column_to_rownames rownames_to_column
#' @return
#' A data.frame, the DESeq2 different expression analysis output data.frame
#' The DESeq2 different expression analysis output contains several columns.
#' The “baseMean” column contains the mean of normalized counts for each gene across all samples.
#' The “log2FoldChange” column contains the log2 fold change between two groups1.
#' The “lfcSE” column contains the standard error of the log2 fold change1.
#' The “stat” column contains the Wald test statistic. The Wald test is commonly used for hypothesis testing when comparing two groups.
#'    And “stat” column is determined by the log2FoldChange and pvalue column values which makes it a good ranking index for further order based analysis like GSEA pre-ranked procedure.
#' The “pvalue” column contains the p-value of the Wald test.
#' The “padj” column contains the adjusted p-value for multiple testing correction.
#' @export
#'
#' @examples
ynormalize_count_bydeseq2 = function(count_file_
                                     ,group_file_
                                     ,cnt_
                                     ,glx_
                                     ,frm='export/'
                                     ,levels = NULL
                                     ,export=F
                                     ,outputdir='export/'
                                     ,norm=T
                                     ,diff=T
){
    # "Tumor_Sample_Barcode and Clin_classification should in [gl]"
    # "IF output==TRUE, output sig_table_new.dfx and colData.dfx"
    print('-------- Please library("DESeq2") first ! -------')
    # ylog('read table ',count_file_)

    cnt = yload_dfx(count_file_,row.names=NULL,frm=frm,ext='dfx')
    # check whether the first col in cnt is gene name or not, if true convert it to rownames
    if (cnt[,1] %>% typeof != 'integer') cnt = cnt %>% column_to_rownames({{cnt}} %>% colnames %>% `[`(1))
    if ((cnt %>% apply(2, is.numeric) == FALSE) %>% sum != 0) stop('input count matrix must be numeric (only can except for the first column)')
    # conver cnt to matrix
    cnt = cnt %>% as.matrix
    patient_n = cnt %>% colnames %>% length
    # filter values out of range
    cnt <- cnt[apply(cnt, 1, sum) > 2*patient_n, ]

    # ylog('reading table ',group_file_)
    glx = yload_dfx(group_file_,row.names=NULL,frm=frm,ext='dfx')

    tsb_in_glx = glx$Tumor_Sample_Barcode
    tsb_in_count = cnt %>% colnames

    tsb_x = tsb_in_count %>%intersect(tsb_in_glx) %>% sort
    if(glx %>% has_rownames){rownames(group_table)=NULL}

    grps = glx$Clin_classification %>% unique %>% sort
    if (!is.null(levels)){
        lvls = levels
    }else{
        lvls = grps
    }
    colData2 = glx %>% dplyr::select(c('Tumor_Sample_Barcode','Clin_classification')) %>%
        filter(Tumor_Sample_Barcode %in% tsb_x) %>%
        as.data.frame %>%
        column_to_rownames("Tumor_Sample_Barcode") %>%
        mutate(Clin_classification= factor(Clin_classification, levels=lvls)) %>%
        relocate(Clin_classification=Clin_classification)

    tsb_x = colData2 %>% rownames
    count_table = cnt[,tsb_x]


    x = list(input_count_table=count_table
             ,col_data=colData2
             ,input_group_list=glx
    )

    if (norm==TRUE){
        # DO normalization
        ddsNormMatrix = DESeqDataSetFromMatrix(count_table, colData = colData2, design= ~ 1)
        tryrrrr <- try({
            vsd <- vst(ddsNormMatrix,nsub = min(1000,dim(ddsNormMatrix)[[1]]))
            'OK'
        },silent = TRUE)
        if (tryrrrr %>% strsplit('\\n') %>% unlist=='OK'){
            # do nothing
        } else if ((tryrrrr %>% strsplit('\\n') %>% unlist)[[3]] == '  it is recommended to use varianceStabilizingTransformation directly'){
            write("[warning] in func vsd: less than 'nsub' rows with mean normalized count > 5, use varianceStabilizingTransformation directly",2)
            vsd =  varianceStabilizingTransformation(ddsNormMatrix)
        }
        sig_table_new=assay(vsd)
        x$normd_count=sig_table_new
        # x$dds_norm_matrix=ddsNormMatrix
    }
    if (diff==TRUE){
        # DO diff expr analysis
        #
        # design, the formula design, is a named vector with values is the TSB, and names is the grouping like
        # c(余先梅='GBM',刘利新='Midline',刘杰='Midline',刘锡全='GBM',...)
        ddsDifExpMatrix = DESeqDataSetFromMatrix(count_table, colData = colData2, design= ~ Clin_classification)
        lvls = lvls[1:2]
        print(paste('Compare between',lvls %>% str_flatten(collapse = " vs "),'. NOTE: grps with more than 2 group will only compare the first 2 groups'))
        diff_expr = ddsDifExpMatrix %>% DESeq %>% results(contrast=c('Clin_classification',lvls))
        x$diff_expr=diff_expr
        # x$dds_dif_exp_matrix=ddsDifExpMatrix
        x$compare_order=lvls
        x$tag_to_add = lvls %>% str_flatten(collapse = '_vs_')
        x$diff_expr_details = x$diff_expr %>% data.frame %>%
            mutate(log2FC_abs = abs(log2FoldChange),.after = log2FoldChange) %>%
            mutate(FC_Ins = log2FoldChange >= 0) %>%
            arrange(desc(log2FC_abs))
    }
    print('Done analysis.')

    # return
    x
}


#' A data.frame
#' The “stat” column in DESeq2 different expression analysis output refers to the Wald test statistic. The Wald test is commonly used for hypothesis testing when comparing two groups. And it is determined by the log2FoldChange and pvalue column values which makes it a good ranking index for further order based analysis like GSEA pre-ranked procedure.
#' The DESeq2 different expression analysis output contains several columns. The “baseMean” column contains the mean of normalized counts for each gene across all samples1. The “log2FoldChange” column contains the log2 fold change between two groups1. The “lfcSE” column contains the standard error of the log2 fold change1. The “stat” column contains the Wald test statistic2. The “pvalue” column contains the p-value of the Wald test2. The “padj” column contains the adjusted p-value for multiple testing correction3.
#'
#' @param cnt.sorted
#' @param colData.sorted
ydo_count_normlization_deseq2 = function(cnt.sorted,colData.sorted){
    # DO normalization
    ddsNormMatrix = DESeq2::DESeqDataSetFromMatrix(cnt.sorted, colData = colData.sorted, design= ~ 1)
    tryrrrr <- try({
        vsd <- DESeq2::vst(ddsNormMatrix,nsub = min(1000,dim(ddsNormMatrix)[[1]]))
        'OK'
    },silent = TRUE)
    if (tryrrrr %>% strsplit('\\n') %>% unlist=='OK'){
        # do nothing
    } else if ((tryrrrr %>% strsplit('\\n') %>% unlist)[[3]] == '  it is recommended to use varianceStabilizingTransformation directly'){
        write("[warning] in func vsd: less than 'nsub' rows with mean normalized count > 5, use varianceStabilizingTransformation directly",2)
        vsd =  DESeq2::varianceStabilizingTransformation(ddsNormMatrix)
    }
    sig_table_new=assay(vsd)

    # x$dds_norm_matrix=ddsNormMatrix

    # normd_count=sig_table_new
    sig_table_new
}



#' Title
#'
#' @param cnt.sorted
#' @param colData.sorted
#'
#' @return
#' @export
#'
#' @examples
ydo_count_diffexpr_deseq2 = function(cnt.sorted,colData.sorted){
    x = alist()
    # DO diff expr analysis
    # design, the formula design, is a named vector with values is the TSB, and names is the grouping like
    # c(余先梅='GBM',刘利新='Midline',刘杰='Midline',刘锡全='GBM',...)
    ddsDifExpMatrix = DESeq2::DESeqDataSetFromMatrix(cnt.sorted, colData = colData.sorted, design= ~ Clin_classification)
    stopifnot(is.factor(colData.sorted$Clin_classification))
    lvls = levels(colData.sorted$Clin_classification)
    lvls = lvls[1:2]
    print(paste('Compare between',lvls %>% str_flatten(collapse = " vs "),'. NOTE: grps with more than 2 group will only compare the first 2 groups'))
    diff_expr = ddsDifExpMatrix %>% DESeq2::DESeq() %>% DESeq2::results(contrast=c('Clin_classification',lvls))
    x$diff_expr=diff_expr
    # x$dds_dif_exp_matrix=ddsDifExpMatrix
    x$compare_order=lvls
    x$tag_to_add = lvls %>% str_flatten(collapse = '_vs_')
    x$diff_expr_details = x$diff_expr %>% data.frame %>%
        mutate(log2FC_abs = abs(log2FoldChange),.after = log2FoldChange) %>%
        mutate(FC_Ins = log2FoldChange >= 0) %>%
        arrange(desc(log2FC_abs))
    x
}
#' Title
#'
#' @param cnt
#' @param colData
#' @param col.id
#' @param col.group
#' @param levels
#' @param auto.filter
#'
#' @return
#' @export
#'
#' @examples
yutils_wash_cnt_and_colData = function(cnt
                                       ,colData
                                       ,col.id = 'Tumor_Sample_Barcode'
                                       ,col.group = 'Clin_classification'
                                       ,levels = NULL
                                       ,auto.filter = TRUE){
    # check whether the first col in cnt is gene name or not, if true convert it to rownames
    #
    if (is.data.frame(cnt)){
        if (!yhas_rownames(cnt)) {cnt = cnt %>% column_to_rownames({{cnt}} %>% colnames %>% `[`(1)) %>% as.matrix} else {cnt = as.matrix(cnt)}
    }else if (is.matrix(cnt)){
        # pass
    }else{
        stop('wrong @cnt type, choose numeric data.frame or matrix')
    }
    # conver cnt to matrix
    # cnt = cnt %>% as.matrix

    if (auto.filter==TRUE){
        patient_n = cnt %>% colnames %>% length
        # filter values out of range
        print('filter out genes with sum(expression) < 2*patient_n, applying `cnt[apply(cnt, 1, sum) > 2*patient_n, ]`')
        cnt <- cnt[apply(cnt, 1, sum) > 2*patient_n, ]
    }
    named.colData = colData %>% has_rownames
    if(!named.colData){
        tsb_in_colData = colData[[col.id]]
    }else{
        tsb_in_colData = rownames(colData)
    }
    tsb_in_count = cnt %>% colnames
    tsb_x =  intersect(tsb_in_count,tsb_in_colData) %>% sort

    grps = colData[[col.group]] %>% unique %>% sort
    if (!is.null(levels)){
        lvls = levels
    }else{
        lvls = grps
    }
    col.id = sym(col.id)
    col.group = sym(col.group)
    if (!named.colData){
        colData2 = colData %>% dplyr::select(Tumor_Sample_Barcode=!!col.id,Clin_classification=!!col.group) %>%
            filter(Tumor_Sample_Barcode %in% tsb_x) %>%
            as.data.frame %>%
            column_to_rownames("Tumor_Sample_Barcode") %>%
            mutate(Clin_classification= factor(Clin_classification, levels=lvls)) %>%
            relocate(Clin_classification=Clin_classification)
        colData2 = colData2[tsb_x,,drop=FALSE]
    }else{
        colData2 = colData[tsb_x,,drop=FALSE] %>%
            dplyr::mutate(Clin_classification=factor(!!col.group,levels=lvls)) %>%
            dplyr::select(Clin_classification)
    }
    count_table = cnt[,tsb_x,drop=FALSE]
    return(list(cnt=count_table,colData=colData2))
}

#' Title
#'
#' @param cnt
#' @param colData
#' @param col.id
#' @param col.group
#' @param levels
#' @param auto.filter
#' @param mode
#'
#' @return
#' @export
#'
#' @examples
ydo_count_deseq2 = function(cnt
                            ,colData
                            ,col.id = 'Tumor_Sample_Barcode'
                            ,col.group = 'Clin_classification'
                            ,levels = NULL
                            ,auto.filter = TRUE
                            ,mode = "auto"
){

    r = yutils_wash_cnt_and_colData(cnt,colData
                                    ,col.id = col.id, col.group = col.group, levels = levels, auto.filter = auto.filter)
    count_table = r[[1]]
    colData2 = r[[2]]

    if (mode=='auto'){
        lvls = base::levels(colData2$Clin_classification)
        norm = TRUE
        if (length(lvls) == 2) {diff = TRUE} else {diff = FALSE}
    }else if (mode == 'norm'){
        norm = TRUE
        diff = FALSE
    }else if (mode == 'diff'){
        norm = FALSE
        diff = TRUE
    }else if(mode == 'all' || mode =='both'){
        norm = TRUE
        diff = TRUE
    }else{
        stop('wrong @mode values, should in c(NULL,"auto","norm","diff","all","both")')
    }
    print(paste('norm =',norm,", diff =",diff))

    x = alist(input_count_table=count_table
              ,col_data=colData2
              ,input_group_df=colData
    )

    if (norm==TRUE){
        x$normd_count = ydo_count_normlization_deseq2(count_table,colData2)
    }
    if (diff==TRUE){
        r = ydo_count_diffexpr_deseq2(count_table,colData2)
        x = utils::modifyList(x,r)
    }
    print('Done analysis.')

    # return
    x
}
