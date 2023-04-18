library(magrittr)

ydumpto = function(x,fname=NULL,export=NULL,outputdir=OUTPUTROOT,ext=NULL,
                   prefix=NULL,flag=NULL,
                   worker=NULL,verbose=T,
                   suffix=NULL,mkdir=F,
                   ...) {
    # OUTER_USAGE = c('fname',"fname",'export',
    #                 'outputdir','ext','prefix',
    #                 'flag','worker','verbose',
    #                 'suffix','mkdir')

    if(is.null(fname)&&is.null(export)) stop('fname export both NULL')
    if (!is.null(export)){
        if (export == TRUE){
            fname = 'plot'
        } else if (export == FALSE){
            return("export==FALSE, exporting aborted")
        } else if (typeof(export) == 'character'){
            fname = export %>% str_flatten()
        } else {
            stop('export must be boolean or string')
        }
    }

    # fname contains path '/'
    if (fname %>% str_detect('/')) {
        tmp = fname %>% ysplit_path()
        if (length(tmp) > 1){
            outputdir = tmp[[1]]
            fname = tmp[[2]]
        }else{
            fname = tmp[[2]]
        }
    }

    # fname contains ext i.e. '.'
    if (fname %>% str_detect('.')) {
        tmp = fname %>% ysplit_file_name(mode = 'auto', ext = ext)
        if (length(tmp)>1){
            fname = tmp[[1]]
            ext = tmp[[2]]
        }else{
            fname = tmp[[1]]
            ext = NULL
        }
    }

    # labels order : [worker,prefix,tag,export(fname),suffix]
    fname = c(worker,flag, prefix,fname,suffix)  %>%
        str_flatten(collapse = '_') %>%
        str_replace_all(pattern = fixed('__'),'_') %>%
        str_remove_all(pattern='^_+|_+$')

    if ((!is.null(ext)) && ext=='json') {
        attr(x,'class') = c('json',class(x))
    }
    # out is generated
    #return(list(fname=fname,ext=ext,outputdir=outputdir,verbose=verbose))
    # if (class(x)=='function') .ydumpto(x,fname,ext,outputdir,verbose=verbose,args=args,...)
    if (mkdir==TRUE && !file.exists(outputdir)) {
        dir.create(outputdir,recursive = TRUE,mode = '0755')
    } else if (mkdir==FALSE && !file.exists(outputdir)) {
        ylog(outputdir,'does not exist.',.addtime = FALSE)
        stop('!')
    }
    argv = yget_args(...)
    if (is.null(x)){
        return(argv)
    }else{
        .ydumpto %>%  do.call(argv)
    }
}
# generic function for export
.ydumpto = function(x,...){
    UseMethod('.ydumpto')
}

.ydumpto.call = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,row.names=NULL,...){
    if (row.names %>% is.null) row.names= FALSE
    if (is.null(ext)) ext='pdf'
    file = ypath_join(outputdir,paste(fname,ext,sep='.'))

    if(ext %in% c('png','svg')){
        dev = get(ext,.GlobalEnv)
    } else if(ext =='pdf') {
        dev = get('cairo_pdf',.GlobalEnv)
    }

    tryCatch({
        dev(file,...)
        eval(x)
        ylog('[plot args]',args %>% names,',',.addtime = F)
    }, finally = {
        dev.off()
    })

}
.ydumpto.data.frame = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,row.names=NULL,...){
    if (row.names %>% is.null) row.names= FALSE
    if (is.null(ext)) ext='dfx'
    file = ypath_join(outputdir,paste(fname,ext,sep='.'))
    argv =  yget_args(...,.f = write.table)
    if (argv$quote |> is.null()) argv$quote=FALSE
    if (argv$sep |> is.null()) argv$sep='\t'
    #     print(c("[argv used]",names(argv)))
    write.table %>% do.call(args=argv)
    ylog('write 1',ext,'at', file, verbose=verbose)
    argv$file
}

.ydumpto.ggplot = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,...){
    if (is.null(ext)) ext='pdf'
    filename = ypath_join(outputdir,paste(fname,ext,sep='.'))
    if(ext %in% c('png','svg')){
        dev = get(ext,.GlobalEnv)
    } else if(ext =='pdf') {
        dev = get('cairo_pdf',.GlobalEnv)
    }
    argv = yget_args(...,.filter = dev)
    tryCatch({
        dev |> do.call(argv)
        print(x)
    }, finally = {
        dev.off()
    })
    ylog('write 1',ext,'at', filename, verbose=verbose)
}

.ydumpto.pheatmap = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,...){
    if (is.null(ext)) ext='pdf'
    #     filename = ypath_join(outputdir,paste(fname,ext,sep='.'))
    filename = ypath_join(outputdir,paste(fname,ext,sep='.'))
    if(ext %in% c('png','svg')) {
        dev = get(ext,.GlobalEnv)
    } else if (ext =='pdf') {
        dev = get('cairo_pdf',.GlobalEnv)
    }
    argv = yget_args(...,.filter = dev)
    tryCatch({
        dev |> do.call(argv)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
    }, finally = {
        dev.off()
    })
    ylog('write 1',ext,'at', filename, verbose=verbose)
}

.ydumpto.matrix = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,row.names=NULL,...){
    if (row.names %>% is.null) row.names= FALSE
    if (is.null(ext)) ext='dfx'
    file  = ypath_join(outputdir,paste(fname,ext,sep='.'))
    argv = yget_args(..., .filter = write.table)
    path = ypath_join(outputdir,paste(fname,ext,sep='.'))
    if (argv$quote |> is.null()) argv$quote=FALSE
    if (argv$sep |> is.null()) argv$sep='\t'
    ylog('[argv used]',names(argv),',',.addtime = F)
    write.table |> do.call(args=argv)
    ylog('write 1',ext,'at', path, verbose=verbose)
}

.ydumpto.json = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,...){
    if (is.null(ext)) ext='json'
    path = ypath_join(outputdir,paste(fname,ext,sep='.'))
    argv = yget_args(...,.f = jsonlite::write_json)
    if (argv$auto_unbox |> is.null()) argv$auto_unbox = TRUE
    ylog('[argv used]',names(argv),',',.addtime = F)
    if (typeof(x)=='character'){
        fileConn<-file(path)
        writeLines(x, fileConn)
        close(fileConn)
    }else{
        jsonlite::write_json |> do.call(args = argv)
    }
    ylog('write 1 json at',path,verbose=verbose)
}

# args is the list of args used by the plotting function
.ydumpto.function = function(x,args,verbose,ext,outputdir,fname,...){
    " return plot_func and its args used for plotting named 'plot_func' & 'plot_args' "
    argv = yget_args(...)
    if ('args' %in% (argv %>% names)){
        print('⚠ run .ydumpto.function with no @args argument supplied')
        args = argv %>% ygetlast('args',.squeeze = T)
        argv = argv %>% yrmlast(c('args','x'))
    }else{
        args = NULL
        argv = argv %>% yrmlast(c('x'))
    }

    if (verbose==TRUE) ylog('you are dumping a function, assuming its a plotting one',.addtime = F)
    if (is.null(ext)) {
        ext='pdf'
    } else if (nchar(as.character(ext)) > 20) {
        print(paste0('There maybe somthing wrong with <ext>, it is to long: "',str_sub(as.character(ext),1,20),'..." total ',nchar(as.character(ext)),' chars'))
        stop('!')
    }
    filename = ypath_join(outputdir,paste(fname,ext,sep='.'))
    argv$filename = filename
    if(ext %in% c('png','svg')){
        dev = get(ext,.GlobalEnv)
    } else if(ext =='pdf') {
        dev = get('cairo_pdf',.GlobalEnv)
    }
    tryCatch({
        argv1 = argv %>% ygetlast(formals(dev) %>% names)
        ylog('[dev args]', argv1 %>% names,',',.addtime = F)
        dev %>% do.call(args = argv1)
        ylog('[plot args]',args %>% names,',',.addtime = F)
        res = x %>% do.call(args)
    }, finally = {
        dev.off()
    })
    ylog('write 1', ext, 'at', filename, verbose=verbose)
    list(plot_func=x,plot_args=args)
}

library(tidyverse)
library(scales)

setwd('~/Downloads/')
getwd()
npg = pal_npg()(10)

# 纵图
for (i in c('104','100')){
    x = read.csv(paste0('./lines.txt'),sep='\t')
    x[3,3] = as.numeric(i) / 100

    gg = x %>% ggplot(aes(x=Group,y=Proliferation)) +
        geom_bar(stat = 'identity',alpha=1,aes(fill=Group),show.legend = FALSE) +
        facet_wrap(~ Cell.lines , nrow = 1,scales = 'free_x', strip.position = 'bottom') +
        theme_bw(base_size = 10) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 2))+
        theme(
            panel.spacing.y = unit(0.5,'inch')
            ,axis.title.y = element_blank()
            ,panel.border = element_blank()
            ,strip.background = element_blank()
            ,strip.placement = "outside"
            ,strip.text = element_text(face='bold',size = 8)
        ) +
        scale_fill_manual(values = npg[c(6,3)])
    # ylab('Clone Formation Rate')

    print(gg)
    # pdf(paste0('lines_',i,'.pdf'),width = 7,height = 4)
    # print(gg)
    # dev.off()
}
# 横图
for (i in c('104','100')){
    x = read.csv(paste0('./lines.txt'),sep='\t')
    x[3,3] = as.numeric(i) / 100

    gg = x %>% ggplot(aes(x=Proliferation,y=Group)) +
            geom_bar(stat = 'identity',alpha=1,aes(fill=Group),show.legend = FALSE) +
            facet_wrap(Cell.lines ~ ., ncol = 1,scales = 'free_y', strip.position = 'left') +
            theme_bw(base_size = 10) +
            scale_x_continuous(labels = scales::percent_format(accuracy = 2))+
            theme(
                panel.spacing.y = unit(0.5,'inch')
                  ,axis.title.y = element_blank()
                  ,panel.border = element_blank()
                  ,strip.background = element_blank()
                  ,strip.placement = "outside"
                  ,strip.text = element_text(face='bold',size = 8)
                  ) +
            scale_fill_manual(values = npg[c(6,3)])
            # ylab('Clone Formation Rate')

    print(gg)
    pdf(paste0('lines_',i,'.pdf'),width = 7,height = 4)
    print(gg)
    dev.off()
}

mat = x %>%
        pivot_wider(id_cols = Cell.lines,values_from = Clone.formation.rate,names_from = Group) %>%
        column_to_rownames('Cell.lines')%>%
        as.matrix
mat = mat * 100

st = chisq.test(t(mat))
st


#' Title
#'
#' @param cnt.sorted
#' @param colData.sorted
#'
#' @return
#' @export
#'
#' @examples
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

    cnt = yload_dfx(count_file_,frm=frm,ext='dfx')
    # check whether the first col in cnt is gene name or not, if true convert it to rownames
    if (cnt[,1] %>% typeof != 'integer') cnt = cnt %>% column_to_rownames({{cnt}} %>% colnames %>% `[`(1))
    if ((cnt %>% apply(2, is.numeric) == FALSE) %>% sum != 0) stop('input count matrix must be numeric (only can except for the first column)')
    # conver cnt to matrix
    cnt = cnt %>% as.matrix
    patient_n = cnt %>% colnames %>% length
    # filter values out of range
    cnt <- cnt[apply(cnt, 1, sum) > 2*patient_n, ]

    # ylog('reading table ',group_file_)
    glx = yload_dfx(group_file_,frm=frm,ext='dfx')

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

