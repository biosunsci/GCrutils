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
        ddsNormMatrix = DESeq2::DESeqDataSetFromMatrix(count_table, colData = colData2, design= ~ 1)
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
        sig_table_new=SummarizedExperiment::assay(vsd)
        x$normd_count=sig_table_new
        # x$dds_norm_matrix=ddsNormMatrix
    }
    if (diff==TRUE){
        # DO diff expr analysis
        #
        # design, the formula design, is a named vector with values is the TSB, and names is the grouping like
        # c(余先梅='GBM',刘利新='Midline',刘杰='Midline',刘锡全='GBM',...)
        ddsDifExpMatrix = DESeq2::DESeqDataSetFromMatrix(count_table, colData = colData2, design= ~ Clin_classification)
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
    sig_table_new=SummarizedExperiment::assay(vsd)

    # x$dds_norm_matrix=ddsNormMatrix

    # normd_count=sig_table_new
    sig_table_new
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
        r = ydo_count_diffexpr_deseq2(count_table,colData2,)
        x = utils::modifyList(x,r)
    }
    print('Done analysis.')

    # return
    x
}


#' DO diff expr analysis of paired or unpaired data
#'
#' @description assert colnames(cnt.sorted)==rownames(colData.sorted)
#'
#' @param cnt.sorted rows are genes, columns are TSBs (Tumor Sample Barcode of samples)
#' @param colData.sorted rows are TSBs and columns are groups, 列名至少包含['Tumor_Sample_Barcode','Clin_classification'] for
#'   @paired=F 列名至少包含['ap','a','g'] for @paired=T, make sure all this columns are not factors!
#' @param paired Do paired test or un-paired test
#' @param levels set compare order, the first in @levels is used as the control group.
#' @param ...
#'
#' @return list named c(diff_expr,compare_order,tag_to_add,diff_expr_details)
#' @export
#'
#' @examples
ydo_count_diffexpr_deseq2 = function(cnt.sorted,colData.sorted, col.id='Tumor_Sample_Barcode'
                                     ,col.group='Clin_classification',paired=FALSE,levels=NULL,.retSimpleTable=TRUE,...){

    if (is.null(levels)){
        if (paired==TRUE){
            stopifnot('g' %in% colnames(colData.sorted))
            if (col.group=='Clin_classification'){
                col.group = 'g'
            }
        }
        if (is.factor(colData.sorted[[col.group]])){
            print(paste('warning',col.group,'is factor convert it into charactors'))
            colData.sorted[[col.group]] = as.character(colData.sorted[[col.group]])
        }
        levels = colData.sorted[[col.group]] %>% unique
    }

    if (colData.sorted %>% has_rownames) {
        stopifnot((colnames(cnt.sorted) == rownames(colData.sorted)) %>% all)
        colData.sorted = colData.sorted %>% rownames_to_column(col.id)
    }else{
        stopifnot((colnames(cnt.sorted) == colData.sorted[[col.id]]) %>% all)
    }

    if (is.data.frame(cnt.sorted)){
        cnt.sorted = as.matrix(cnt.sorted)
    }
    sym.col.id = sym(col.id)
    sym.col.group = sym(col.group)

    if (paired==TRUE){
        #         stopifnot(colData.sorted %>% has_colnames('a'))
        #         stopifnot(colData.sorted %>% has_colnames('g'))
        colData = colData.sorted %>%
            select(!!sym.col.id,a,g) %>%
            mutate(a=factor(a),g=factor(g,levels=levels)) %>%
            arrange(a,g) %>%
            column_to_rownames(col.id)
    }else{
        colData = colData.sorted %>%
            select(a=!!sym.col.id,!!sym.col.group) %>%
            mutate(g=factor(!!sym.col.group,levels=levels)) %>%
            arrange(a,g) %>%
            column_to_rownames("a")
    }

    cnt.sorted = cnt.sorted[,row.names(colData)]

    x = alist()
    x$compare_ref_group = levels[[1]]
    x$compare_order = rev(levels)
    x$tag_to_add = x$compare_order %>% str_flatten(collapse = '_vs_')

    if(paired==TRUE){
        dds_paired  = DESeq2::DESeqDataSetFromMatrix(cnt.sorted, colData = colData, design= ~ a + g)
        paired_prop = 'paired'
    }else{
        dds_paired  = DESeq2::DESeqDataSetFromMatrix(cnt.sorted, colData = colData, design= ~ g)
        paired_prop = 'unpaired'
    }
    print('-------------------------------------------------------------')
    print(paste0('Comparing > ',paired_prop,' < ',x$tag_to_add,', REF_GROUP = ',x$compare_ref_group))
    print('-------------------------------------------------------------')

    # set control group for comparation
    dds_paired$group <- relevel(dds_paired$g, ref = levels[[1]])
    dds_paired <- DESeq2::DESeq(dds_paired)
    # get result table
    diff_expr <- as.data.frame(DESeq2::results(dds_paired)) %>%
        mutate(abs_log2FC = abs(log2FoldChange),.before = lfcSE) %>%
        arrange(pvalue,log2FoldChange)
    # pack add info
    x$cnt_table = cnt.sorted
    x$colData = colData
    x$diff_expr = diff_expr
    x$paired = paired
    # x$dds_dif_exp_matrix=ddsDifExpMatrix
    class(x) <- c('DEA.results',class(x))
    print('Done')
    return (x)
}


# library(limma)

#' @description  like ydo_count_diffexpr_deseq2,  but used with numeric value but not count values
#'
#' @param numeric.sorted rows are genes, columns are TSBs (Tumor Sample Barcode of samples)
#' @param colData.sorted rows are TSBs and columns are groups, 列名至少包含[col.id = 'Tumor_Sample_Barcode',col.group = 'Clin_classification']
#' @param col.id id as rownames of colData.sorted, usually Tumor_Sample_Barcode
#' @param col.group group column in colData.sorted
#' @param levels compare order, The first level in levels will be control group(Denominator) in DEG
#'   results. the first one is used as the control CTL, for example, if levels = c(CTL, EXP) then
#'   the result log2FoldChange > 0 if EXP > CTL
#'
#' @return
#' @export
#'
#' @examples
ydo_numeric_DEG_limma = function(numeric.sorted, colData.sorted,col.id = 'Tumor_Sample_Barcode', col.group = 'Clin_classification', levels=NULL){

    if (!colData.sorted %>% has_rownames){
        colData.sorted = colData.sorted %>% column_to_rownames(col.id)
    }
    stopifnot((colnames(numeric.sorted) == rownames(colData.sorted)) %>% all)
    if (is.data.frame(numeric.sorted)){
        numeric.sorted = as.matrix(numeric.sorted)
    }
    stopifnot(col.group %in% colnames(colData.sorted))
    if (is.null(levels)){
        levels = colData.sorted[[col.group]] %>% unique()
        print(levels)
    }else{
        .i = levels[1]
        levels = c(levels[2:length(levels)], .i)
    }
    contrast_order = stringr::str_flatten(levels, collapse = ' - ')

    print("-------------------------------------------------------------")
    print(paste0("Comparing > ", 'unpaired', " < ", contrast_order,
                 ", REF_GROUP = ", levels[length(levels)]))
    print("-------------------------------------------------------------")

    g = colData.sorted[[col.group]] %>% factor(levels=levels)

    design <- stats::model.matrix(~ g - 1)
    colnames(design) = levels

    fit <- limma::lmFit(numeric.sorted, design)
    contrast.matrix <- limma::makeContrasts(contrasts = contrast_order, levels = design)

    fit2 = limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    results <- limma::topTable(fit2, coef=1, number=Inf)

    colnames(results) = c("log2FoldChange",'baseMean','t','pvalue','padj','B')
    # c('baseMean','log2FoldChange','log2FC_abs','lfcSE','stat','pvalue','padj','FC_Ins')
    #     print('assign stat = t')
    results$stat = results$t
    results$abs_log2FC = abs(results$log2FoldChange)
    results$FC_Ins = results$log2FoldChange > 0

    # cal mean of each group
    for (i in levels(g)){
        mat = numeric.sorted[,g == i]
        results[[paste0("median_",i)]] = matrixStats::rowMedians(mat)
    }
    #     mean.df = apply(numeric.sorted, 1, function(x) tapply(x, g, mean, na.rm = TRUE)) %>% t %>% as.data.frame
    #     results = cbind(results,mean.df)
    results %>% rownames_to_column('symbol')
}

