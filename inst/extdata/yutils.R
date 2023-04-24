SSGSEA_BASE_REF_DIR = '/share_storage/pipeline/REF_FILES/'

library(GenomicRanges)
#' Title
#'
#' @param w
#'
#' @return
#' @export
#'
#' @examples
yget_chr_window = function(w = 1e6){
    chr_df = circlize::read.chromInfo()$df
    chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
    chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))

    chr_window = EnrichedHeatmap::makeWindows(chr_gr, w = w)
    chr_window
}

#' Title
#'
#' @param window
#' @param gr
#' @param v
#' @param method
#' @param empty_v
#'
#' @return
#' @export
#'
#' @examples
yaverage_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {

    if(missing(v)) v = rep(1, length(gr))
    if(is.null(v)) v = rep(1, length(gr))
    if(is.atomic(v) && is.vector(v)) v = cbind(v)

    v = as.matrix(v)
    if(is.character(v) && ncol(v) > 1) {
        stop("`v` can only be a character vector.")
    }

    if(length(empty_v) == 1) {
        empty_v = rep(empty_v, ncol(v))
    }

    u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))

    mtch = as.matrix(findOverlaps(window, gr))
    intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
    w = width(intersect)
    v = v[mtch[,2], , drop = FALSE]
    n = nrow(v)

    ind_list = split(seq_len(n), mtch[, 1])
    window_index = as.numeric(names(ind_list))
    window_w = width(window)

    if(is.character(v)) {
        for(i in seq_along(ind_list)) {
            ind = ind_list[[i]]
            if(is.function(method)) {
                u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
            } else {
                tb = tapply(w[ind], v[ind], sum)
                u[window_index[i], ] = names(tb[which.max(tb)])
            }
        }
    } else {
        if(method == "w0") {
            gr2 = reduce(gr, min.gapwidth = 0)
            mtch2 = as.matrix(findOverlaps(window, gr2))
            intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])

            width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
            ind = unique(mtch2[, 1])
            width_setdiff = width(window[ind]) - width_intersect

            w2 = width(window[ind])

            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
                u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
            }

        } else if(method == "absolute") {
            for(i in seq_along(ind_list)) {
                u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
            }

        } else if(method == "weighted") {
            for(i in seq_along(ind_list)) {
                ind = ind_list[[i]]
                u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
            }
        } else {
            if(is.function(method)) {
                for(i in seq_along(ind_list)) {
                    ind = ind_list[[i]]
                    u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
                }
            } else {
                stop("wrong method.")
            }
        }
    }

    return(u)
}




# yhas_rownames = function(df) {
#   !all(row.names(df)==seq(1, nrow(df)))
# }



# CODE IN CREATE PROJECT
library(tidyverse)
# library(ggfortify)
library(ggsci)

if (!exists("WORKER",.GlobalEnv)){
    WORKER = NULL
}







##@





#' Title
#'
#' @param f
#'
#' @return
#' @export
#'
#' @examples
ylocate_funs <- function(f) {
    # Returns dataframe with two columns:
    # `package_name`: packages(s) which the function is part of (chr)
    # `builtin_package`:  whether the package comes with standard R (a 'builtin'  package)

    # Arguments:
    # f: name of function for which the package(s) are to be identified.


    if ("tidyverse" %in% rownames(installed.packages()) == FALSE) {
        cat("tidyverse is needed for this fuction. Please install. Stopping")
        stop()}

    suppressMessages(library(tidyverse))


    # search for help in list of installed packages
    help_installed <- help.search(paste0("^",f,"$"), agrep = FALSE)

    # extract package name from help file
    pckg_hits <- help_installed$matches[,"Package"]

    if (length(pckg_hits) == 0) pckg_hits <- "No_results_found"


    # get list of built-in packages

    pckgs <- installed.packages()  %>% as_tibble
    pckgs %>%
        dplyr::filter(Priority %in% c("base","recommended")) %>%
        dplyr::select(Package) %>%
        distinct -> builtin_pckgs_df

    # check for each element of 'pckg hit' whether its built-in and loaded (via match). Then print results.

    results <- data_frame(
        package_name = pckg_hits,
        builtin_pckage = match(pckg_hits, builtin_pckgs_df$Package, nomatch = 0) > 0,
        loaded = match(paste("package:",pckg_hits, sep = ""), search(), nomatch = 0) > 0
    )

    return(results)
}




#' Title
#'
#' @param x
#' @param func
#'
#' @return
#' @export
#'
#' @examples
yapply = function(x,func){
    l = length(formals(func))
    if (l==3) for (i in 1:length(x)){
        res = func(x[[i]],names(x[i]),i)
    }else if (l==2)  for (i in 1:length(x)){
        res = func(x[[i]],names(x[i]))
    }else if (l==1)  for (i in 1:length(x)){
        res = func(x[[i]])
    }
    res
}
#' Title
#'
#' @param .x
#' @param .f
#' @param .init
#' @param .right
#' @param .accumulate
#'
#' @return
#' @export
#'
#' @examples
yreduce = function(.x,.f,.init,.right=FALSE,.accumulate=FALSE){
    Reduce(.f,.x,init = .init,right=.right, accumulate=.accumulate)
}

#' Title
#'
#' @param x
#' @param func
#'
#' @return
#' @export
#'
#' @examples
ymap = function(x, func) {
    " func (item, name, index) "
    l = length(formals(func))
    res = list()
    if (l == 3){
        for (i in 1:length(x)) {
            n = names(x[i])
            res[[n]] = func(x[[i]], names(x[i]), i)
        }
    }else if (l == 2){
        for (i in 1:length(x)) {
            n = names(x[i])
            res[[n]] = func(x[[i]], names(x[i]))
        }
    } else if (l == 1)
        for (i in 1:length(x)) {
            res[[i]] = func(x[[i]])
        }
    res
}


yconfig_global_load = function(keys=NULL
                               ,.json='global_config.json'
                               ,.frm='/GCI/jup/A_TSF/GLOBAL'
                               ,retmode = 'local'
                               ,verbose = TRUE
){
    # use package jsonlite to load frm/json file, assign its containing dict keys to .GlobalEnv with values
    # at current, 2022年2月19日, only support un-nested json dict
    # @params retmode ['local','global']
    library(jsonlite)
    log = c()
    i = 0
    json_file = yfile_path(.frm,.json)
    if (file.exists(json_file)){
        parameters = jsonlite::fromJSON(json_file)
        if ('list' %in% class(parameters)){
            row = parameters
        }
        if ('data.frame' %in% class(parameters)){
            row = parameters[1,] %>% list
        }
        if (keys %>% is.null == FALSE) params = ygetlast(row,keys)
        else { params = row }

        if (retmode=='global'){
            for (var in names(params)){
                log = c(log,paste0("  ",var," <- '",params[var],"' "))
                assign(var,params[[var]],.GlobalEnv)
                i = i + 1
            }
        }else if (retmode=='local'){
            i = params %>% length
        }
        else stop('wrong <retmode> value')
        ylog("INFO ",json_file," loaded. [",i,'] vars Parsed ',log %>% str_flatten,verbose=verbose)
        if (retmode=='local') return(params)
        NULL
    }else{
        # nothing to do
        ylog('Config_file <',json_file,'> does not exist, Nothing parsed')
    }
}

##@


##@

##@

ycurrent_worker = function (set = NULL){
    if (set %>% is.null) set = getwd()

    if (set %>% str_detect('/GCI/jup/(.+)/{0,2}.*')) key = (set %>% str_split(fixed('/')))[[1]][[4]]
    else if (set %>% str_detect('^[A-Z][0-9]$')) return(set)
    else key = set

    conf = yload_global_config('GET_WORKER_KEY',verbose = F)
    set = conf %>% ygetlast(key)
    if (!is.null(set)) return(set)
    else stop(paste('invalid key [',key,']: worker search key invalid, check whether this key is in db'))
    # stop('wrong set value, NULL or matches ^[A-Z][0-9]\\_$')

}














#' 暂时没有用到@2023-03-09
yutils_df_add_titleline = function(dat, table_title_line) {
    if(is.matrix(dat)){
        dat = data.frame(dat)
    }
    if (yhas_rownames(dat)){
        dat = dat %>% rownames_to_column('index')
    }
    rbind(c(paste0('## ',table_title_line), rep('', ncol(dat)-1)), # title
          # rep('', ncol(dat)), # blank spacer row
          names(dat), # column names
          unname(sapply(dat, as.character))
    )
}





#' subsetMaf using group data.frame glx
#'
#' @param maf `read.maf()` generated maf class object
#' @param glx data.frame contains the Group information
#' @param col the column in `glx` used as Group vector to generate the submaf
#'  default NULL, if NULL check either one of c('Group','Clin_classification','g') is in the colnames(`glx`)
#'  only use the most prior one
#' @param .container default NULL, if not NULL, must be list class obj,
#'  the subgroups will be assign to this list then returned
#' @param .assignGlobal default TRUE, controls the return of the function, see return section
#'
#' @return NULL if `.assignGlobal` == TRUE else list: names is subGroup names and values is submaf
#' @export
#'
#' @examples
ysubset_mafs_ = function(maf,glx=NULL,col = NULL,.container=NULL,.assignGlobal=TRUE){
    warning('USE ygen_subMafs_ with raw data.frame instead, Usage of this function is not recommonded, cause it will delete TSB with no mutations in subgroup')
    if (is.null(glx)){
        glx = maf %>% getClinicalData
    }
    if (!is.null(.container)){
        .assignGlobal = FALSE
        stopifnot(is.list(.container))
    }
    if (is.null(col)){
        col.names = colnames(glx)
        if ('Group' %in% col.names) col = 'Group'
        else if ('Clin_classification' %in% col.names) col = 'Clin_classification'
        else if ('g' %in% col.names) col = 'g'
        else stop('can not auto determine group col arg @col, please set it manually')
    }
    grps = glx[[col]] %>% unique
    if (length(grps)>2){
        submafs = list()
        for (g in grps){
            query = paste0(col," == '",g,"' ")
            print(paste0('Applying filter: ',query))
            submafs[[g]] = maf %>% subsetMaf(clinQuery = query)
        }
        return (submafs)
    }else{
        g1name = grps[[1]]
        g1 = maf %>% subsetMaf(clinQuery = paste0(col," == '",g1name,"' "))
        g2name = grps[[2]]
        g2 = maf %>% subsetMaf(clinQuery = paste0(col," == '",g2name,"' "))
        print(paste0('Applying filter: ',col," == '",g1name,"' or '",g2name,"' "))
        if (.assignGlobal){
            assign('g1',g1,envir = .GlobalEnv)
            assign('g2',g2,envir = .GlobalEnv)
            assign('g1name',g1name,envir = .GlobalEnv)
            assign('g2name',g2name,envir = .GlobalEnv)
            print('[!] 4 Global variables set: g1, g2, g1name, g2name')
        }else if (!is.null(.container)){
            .container$g1 = g1
            .container$g2 = g2
            .container$g1name = g1name
            .container$g2name = g2name
            print('[!]returnning a list .container with 4 variables being set: g1, g2, g1name, g2name')
            return(.container)
        }else{
            submafs = list(g1,g2)
            names(submafs) = c(g1name,g2name)
            return(submafs)
        }
    }
}



#' gen_subMafs_ using data.frame snv/cnv, glx
#'
#' @param mol `read.maf()` generated maf class object
#' @param glx data.frame contains the Group information
#' @param col group_column, the column in `glx` used as Group vector to generate the submaf
#'  default NULL, if NULL check either one of c('Group','Clin_classification','g') is in the colnames(`glx`)
#'  only use the most prior one
#' @param id_col='Tumor_Sample_Barcode', no need to change
#' @param .container default NULL, if not NULL, must be list class obj,
#'  the subgroups will be assign to this list then returned
#' @param .assignGlobal default TRUE, controls the return of the function, see return section
#'
#' @return NULL if `.assignGlobal` == TRUE else list: names is subGroup names and values is submaf
#' @export
#'
#' @examples
ygen_subMafs_ = function(mol,glx,col = NULL,id_col='Tumor_Sample_Barcode',.container=NULL,.assignGlobal=TRUE,...){
    if (!is.null(.container)){
        .assignGlobal = FALSE
        stopifnot(is.list(.container))
    }
    if (class(mol)=='MAF'){
        stop('Please not input MAF read but raw snv/cnv data.frame and raw glx data.frame')
    }
    if (is.null(col)){
        col.names = colnames(glx)
        if ('Group' %in% col.names) col = 'Group'
        else if ('Clin_classification' %in% col.names) col = 'Clin_classification'
        else if ('g' %in% col.names) col = 'g'
        else stop('can not auto determine group col arg @col, please set it manually')
    }
    grps = glx[[col]] %>% unique

    if (length(grps)>2){
        submafs = list()
        for (g in grps){
            mask = glx[[col]] == g
            sub.glx = glx[mask,]
            mask = mol[[id_col]] %in% sub.glx[[id_col]]
            sub.mol = mol[mask,]
            submafs[[g]] = read.maf(sub.mol,sub.glx,...)
        }
        return (submafs)
    }else{
        g1name = grps[[1]]
        sub.glx1 = glx[glx[[col]] == g1name,]
        sub.mol1 = mol[mol[[id_col]] %in% sub.glx1[[id_col]],]
        g1 = read.maf(sub.mol1,sub.glx1,...)
        g2name = grps[[2]]
        sub.glx2 = glx[glx[[col]] == g2name,]
        sub.mol2 = mol[mol[[id_col]] %in% sub.glx2[[id_col]],]
        g2 = read.maf(sub.mol2,sub.glx2,...)
        print(paste0('Applying filter: ',col," == '",g1name,"' or '",g2name,"' "))
        if (.assignGlobal){
            assign('g1',g1,envir = .GlobalEnv)
            assign('g2',g2,envir = .GlobalEnv)
            assign('g1name',g1name,envir = .GlobalEnv)
            assign('g2name',g2name,envir = .GlobalEnv)
            print('[!] 4 Global variables set: g1, g2, g1name, g2name')
        }else if (!is.null(.container)){
            .container$g1 = g1
            .container$g2 = g2
            .container$g1name = g1name
            .container$g2name = g2name
            print('[!]returnning a list .container with 4 variables being set: g1, g2, g1name, g2name')
            return(.container)
        }else{
            submafs = list(g1,g2)
            names(submafs) = c(g1name,g2name)
            return(submafs)
        }
    }
}



library(DESeq2)
library(GSVA)
library(clusterProfiler)

.y_generate_gmt = function(){

}



#' load gmt files using clusterProfile::read.gmt into data.frame format obj
#'
#' @description see docs of yload_list_gmt, yload_gmt
#'
#' @param key one of names(`ref`), default one of c('hm','kegg','immune','immune.jnj'),
#' @param ref key=filename pairs list, the `key` is search in. No need to change
#' @param base_ref_dir =NULL, if NULL use global SSGSEA_BASE_REF_DIR value
#' search gmt files in the base_ref_dir folder
#'
#' @return data.frame with first column is term, and second columns is gene, like
#'                               term    gene
#' 1 HALLMARK_TNFA_SIGNALING_VIA_NFKB    JUNB
#' 2 HALLMARK_TNFA_SIGNALING_VIA_NFKB   CXCL2
#' 3 HALLMARK_TNFA_SIGNALING_VIA_NFKB    ATF3
#' 4 HALLMARK_TNFA_SIGNALING_VIA_NFKB  NFKBIA
#' 5 HALLMARK_TNFA_SIGNALING_VIA_NFKB TNFAIP3
#' ...
#' @export
#'
#' @examples
#' gmt = yload_gmt('hm')
yload_gmt = function(key=NULL,ref = list(
    hm='h.all.v7.5.symbols.gmt'
    ,kegg='c2.cp.kegg.v7.5.symbols.gmt'
    ,immune='immune.cell.rep.symbols.gmt'
    ,immune.jnj='immune.jnj.immunity.symbols.gmt')
    ,base_ref_dir = NULL
){
    if (is.null(base_ref_dir) || is.null(base_ref_dir)){
        base_ref_dir = SSGSEA_BASE_REF_DIR
    }
    if(!is.null(key)){
        stopifnot(key %in% names(ref))
        if (key %in% c('hm','kegg','immune','immune.jnj')){
            v = ref[[key]]
            path = file.path(base_ref_dir,v)
            gmt = clusterProfiler::read.gmt(path)
        }else{
            stop('key not supported')
        }
    }
    gmt
}

#' read gmt file as a list
#'
#' read Gene Matrix Transposed (gmt) file and output a list with the the first
#' column as the names of items in the list. see
#' \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{Gene Matrix Transposed file format}
#' for more details.
#'
#' @param annofile a gmt file. Examples are from MSigDB Collections.
#' A list of gene set could be find in the vignette of cogena
#' @return a gmt list
#' @export
#'
#' @seealso gmtlist2file
#' @examples
#' anno <- "c2.cp.kegg.v7.01.symbols.gmt.xz"
#' annofile <- system.file("extdata", anno, package="cogena")
#' gl <- gmt2list(annofile)
#'
yload_gmt2list <- function(annofile){

    if (!file.exists(annofile)) {
        stop("There is no such a gmt file!")
    }

    if (tools::file_ext(annofile) == "xz") {
        annofile <- xzfile(annofile)
        x <- scan(annofile, what="", sep="\n", quiet=TRUE)
        close(annofile)
    } else if (tools::file_ext(annofile) == "gmt") {
        x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    } else {
        stop ("Only gmt and gmt.xz are accepted for gmt2list")
    }

    y <- strsplit(x, "\t")
    names(y) <- sapply(y, `[[`, 1)

    annoList <- lapply(y, `[`, c(-1,-2))
    return(annoList)
}

#' load gmt files using read.table into list format obj
#'
#' @description see docs of yload_list_gmt, yload_gmt
#'
#' @param key one of names(`ref`), default one of c('hm','kegg','immune','immune.jnj')
#' @param ref key=filename pairs list, the `key` is search in. No need to change
#' @param base_ref_dir =NULL, if NULL use global SSGSEA_BASE_REF_DIR value
#' @param gmt_txt_file_path =NULL, a file path string, file must the same as the files listed in the `ref`
#'  if provided, read the file of the gmt_txt_file_path, ignore all other parameters
#' @return list with names is term, and value is corresponding genes, like
#'  list(HALLMARK_ADIPOGENESIS=c('FABP4','ADIPOQ','PPARG',...),
#'      HALLMARK_ALLOGRAFT_REJECTION=c('PTPRC','IL12B','TGFB1','IL12A','CD3E',...),
#'      HALLMARK_ANDROGEN_RESPONSE=c('KLK3','KLK2','ACSL3','PIAS1','CAMKK2',...),
#'      HALLMARK_ANGIOGENESIS=c('VCAN','POSTN','FSTL1',...),
#'      ...)
#' @export
#'
#' @examples
yload_list_gmt = function(key=NULL,ref = list(
    hm='h_all_v7_5_symbols.txt'
    ,kegg='c2_cp_kegg_v7_5_symbols.txt'
    ,immune='immune_list.txt'
    ,immune.jnj='immune_list_jnj.txt')
    ,base_ref_dir = NULL
    ,gmt_txt_file_path = NULL
){
    if (!is.null(gmt_txt_file_path)){
        stopifnot(file.exists(gmt_txt_file_path))
        path = gmt_txt_file_path
    }else{
        if (is.null(base_ref_dir) || is.null(base_ref_dir)){
            base_ref_dir = SSGSEA_BASE_REF_DIR
        }
        if(!is.null(key)){
            stopifnot(key %in% names(ref))
            # (key %in% c('hm','kegg','immune','immune.jnj'))
            path = file.path(base_ref_dir,ref[[key]])
        }
    }
    tmp = read.table(path
                     ,quote='"'
                     ,header = T
                     ,sep=','
                     ,row.names=1)
    gmt_list = split(tmp[,1],tmp[,2])
    gmt_list
}


#' do_ssGSEA analysis
#'
#' @description depends on yload_list_gmt
#'  use the GSVA::gsea method which need 2 critical parameters
#'  a <expr_matrix> is the numeric matrix of expression, normed count of DSeq2 or tpm
#'  a <ref_gmt_list> which is the returns of yload_list_gmt, see details of yload_list_gmt documentation
#'
#' @param expr_matrix a numeric matrix witch rows are genes and columns are patients
#' @param todo list, combinations of elements in list('hm','kegg','immune','immune.jnj')
#'  the collections todo ssGSEA, default run all 4 process, all the keys is the keys of parameter
#'  <ref> in function `yload_list_gmt`
#' @param ref.gmt =NULL, if provided, must be a list like yload_list_gmt returns
#'  or a file_path to a gmt.txt file like in yload_list_gmt <ref>, use it as <ref_gmt_list> todo ssGSEA
#' @param ref.gmt.name =NULL, if `ref.gmt` is provided, this will be the key-name of the returned results list
#'  if is NULL, the returned list is unamed
#' @param mx.diff,verbose,method,kcdf,parallel.sz,... pass to GSVA::gsva, details see GSVA::gsva docs
#'  please keep `method` 'ssgsea' to do ssGSEA analysis
#' @return list
#' @export
#'
#' @examples
ydo_ssGSEA = function(expr_matrix
                      ,todo = list('hm','kegg','immune','immune.jnj')
                      ,ref.gmt = NULL
                      ,ref.gmt.name = NULL
                      ,mx.diff=FALSE
                      ,verbose=FALSE
                      ,method='ssgsea'
                      ,kcdf='Gaussian'
                      ,parallel.sz=1
                      ,...){
    res = list()
    if (!is.null(ref.gmt)){
        if (is.character(ref.gmt)){
            ref_df = yload_list_gmt(gmt_txt_file_path = ref.gmt)
        }
        if (is.null(ref.gmt.name)){
            ref.gmt.name = 1
        }
        es = GSVA::gsva(expr_matrix
                        , ref_df
                        , mx.diff=mx.diff
                        , verbose=verbose
                        , method=method
                        , kcdf=kcdf
                        , parallel.sz=1)
        res[[ref.gmt.name]] = es
    }else{
        for (i in 1:length(todo)){
            k = todo[[i]]
            if (!is.null(ref.gmt)){
                ref_df = ref.gmt
            }else{
                ref_df = yload_list_gmt(k)
            }
            es = GSVA::gsva(expr_matrix
                            , ref_df
                            , mx.diff=mx.diff
                            , verbose=verbose
                            , method=method
                            , kcdf=kcdf
                            , parallel.sz=1
                            ,...
            )
            # es %>% ydumpto(flag='ssGSEA_result_',export=k,outputdir=OUTPUTROOT,row.names=T)
            res[[k]] = es
        }
    }
    res
}


yutils_filter_stat_res = function(df,condition.str = NULL, t.p.val = NA, t.log.fc = 1,  t.p.adj = 0.05, col.p.val = 'pvalue', col.p.adj = 'padj',col.log.fc = 'log2FC_abs',...){
    if (is.null(condition.str)){
        cols = colnames(df)
        stopifnot(c(col.pval,col.p.adj,col.log.fc) %in% cols)
        str = c()
        for (i in c("p.val","p.adj")){
            t = get(paste0('t.',i))
            if (!(is.na(t) || is.null(t))){
                str = c(str,get(paste0('col.',i)),'<',t,';')
            }
        }
        if (!(is.na(t.log.fc) || is.null(t.log.fc))){
            str = c(str,col.log.fc,'>=',t.log.fc)
        }
        e.str = str_flatten(str)
    }else{
        e.str = condition.str
    }
    e1 = rlang::parse_exprs(e.str)
    e2 = enexprs(...)
    e = c(e1,e2)
    r = df %>% filter(!!! e)
    print(paste0('[!] nrow ',nrow(df),' => ',nrow(r)))
    print(paste0('Filtered by ', e %>% as.character %>% str_flatten_comma()))
    r
}

# GSEA
library(clusterProfiler)
library('org.Hs.eg.db')
library(ggsci)

#



# GO analysis using different expression results
#' Title
#'
#' @param diff a data.frame with columns: symbol, log2FC, pvalue, padj, stat
#' @param ont Ontology, in c('BP', 'MF', 'CC', 'ALL')
#' @param by.expr.direction =TRUE, if TRUE genes in diff is group_by change direction before doing GO analysis
#' @param p.t pvalue threshold, default 0.01
#' @param p.adjust.t padj threshold, default 0.05
#' @param log2fc.abs.t log2FC threshold, default 1
#' @param ...
#'
#' @return a list obj of class 'y.GO.res', used for input for yplot_GO_res
#' @export
#'
#' @examples
ydo_GO_bydiff = function(diff, ont = "ALL", by.expr.direction=FALSE
                         , go.tb.filter = 'top10'
                         , p.t = 0.01, p.adjust.t = 0.05
                         , t.p.val = NA, t.log.fc = 1, t.p.adj = 0.05
                         ,  ...) {
    "ont is Ontology, in c('BP', 'MF', 'CC', 'ALL')"
    ""
    "如果genes is.null, use p.t or p.adjust.t && log2fc.abs.t (as FILTERS) to filter diff to get genes"
    "diff必须包含全部基因(as ALL)而不是基因子集,GO分析才准确(https://mp.weixin.qq.com/s/zxiMXqRriGHcGAlCNk0e7A)\n    ,若genes(格式为symbols)不为空,使用genes/ALL进行GO分析,若genes为空, 使用FILTERS从diff中筛选出genes"

    ont = toupper(ont)
    sig_genes = yutils_filter_stat_res(diff,t.p.val = t.p.val, t.log.fc = t.log.fc, t.p.adj = t.p.adj)

    if (!yhas_rownames(diff)) {
        diff = diff %>% rownames_to_column("symbol")
    }
    all_genes = rownames(diff)

    Name = "GO_enrich"
    tryCatch(expr = {
        GENE.REF = get("GENE.REF", envir = .GlobalEnv)
        print("Got Global var GENE.REF")
    }, error = function(err) {
        assign("GENE.REF", value = bitr(all_genes, fromType = "SYMBOL",
                                        toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
                   arrange(SYMBOL) %>% distinct(SYMBOL, .keep_all = TRUE) %>%
                   column_to_rownames("SYMBOL"), envir = .GlobalEnv)
        print("Global var GENE.REF assigned")
    })

    if (by.expr.direction==TRUE){
        res.tb = list()
        res.ego = list()
        for (direction in c('Up','Down')){
            genes = sig_genes %>% filter(FC_Ins==(direction=="Up")) %>% rownames
            ego = enrichGO(gene = GENE.REF[genes, "ENTREZID"], universe = GENE.REF$ENTREZID,
                           ont = ont, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                           pAdjustMethod = "BH", pvalueCutoff = p.t, qvalueCutoff = p.adjust.t,
                           readable = TRUE, ...)
            res.tb[[direction]] = ego@result %>% mutate(Change = factor(!!direction))
            res.ego[[direction]] = ego
        }
        res.tb = rbind(res.tb[[1]],res.tb[[2]]) %>%
            mutate(Description = fct_reorder(Description, p.adjust, .desc = TRUE)
                   , ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
                   , EnrichedGeneNum = Count
                   ,setSize = BgRatio %>%
                       str_split("/", simplify = TRUE) %>% .[,1] %>% as.integer
                   , EnrichedPer = round(EnrichedGeneNum/setSize * 100))
    }else{
        genes = sig_genes %>% rownames
        res.ego = enrichGO(gene = GENE.REF[genes, "ENTREZID"], universe = GENE.REF$ENTREZID,
                           ont = ont, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                           pAdjustMethod = "BH", pvalueCutoff = p.t, qvalueCutoff = p.adjust.t,
                           readable = TRUE, ...)
        res.tb = res.ego@result %>%
            mutate(Description = fct_reorder(Description, p.adjust, .desc = TRUE)
                   , ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
                   , EnrichedGeneNum = Count
                   ,setSize = BgRatio %>%
                       str_split("/", simplify = TRUE) %>% .[,1] %>% as.integer
                   , EnrichedPer = round(EnrichedGeneNum/setSize * 100)
                   ,Change = 'DiffExpr'
            )%>% dplyr::select(!geneID)
    }
    attr(res.tb,'class') = c(class(res.tb),'y.GO.res.dtable')
    print('GO analysis Done')
    res = list(gsea.set = Name
               , gene_list = sig_genes
               , res.tb = res.tb
               , res.go = res.ego
               , args = list(ont=ont
                             ,by.expr.direction=by.expr.direction
                             ,p.t = p.t
                             ,p.adjust.t = p.adjust.t
                             ,log2fc.abs.t = t.log.fc
               )
    )
    attr(res,'class') = c(class(res),'y.GO.res')
    message('GO analysis done, now you can use yplot_GO_res to get the figure')

    res
}

#' do GO analysis of a vector of gene symbols
#'
#' @param sig_genes sig_genes是all_genes的子集,sig_genes与all_genes的占比决定了GO富集的结果
#' @param all_genes all_genes必须包含全部基因(as ALL)而不是基因子集,GO分析才准确(https://mp.weixin.qq.com/s/zxiMXqRriGHcGAlCNk0e7A),
#'  genes(格式为基因正式的名称symbols,不是ID)不为空,使用genes/ALL进行GO分析, 如果NULL,使用
#' @param ont ont is Ontology, one of c('BP', 'MF', 'CC', 'ALL')
#'  ALL会进行全部3个富集分析:'BP', 'MF', 'CC'
#' @param ggsci_col_scheme 默认使用pal_npg()(10)[c(3,6,9)]来标识BP,MF,CC的富集结果
#' @param ... 其他传入enrichGO函数的参数
#'
#' @return a list obj of class 'y.GO.res', used for input for yplot_GO_res
#' @export
#'
#' @examples
ydo_GO = function(sig_genes
                  ,all_genes =NULL
                  ,ont = 'ALL'
                  ,p.t = 0.01
                  ,p.adjust.t = 0.05
                  ,ggsci_col_scheme = NULL
                  ,...
){
    .T_ylab_wrap = 50
    ont = toupper(ont)
    if (is.null(ggsci_col_scheme)){
        npg = pal_npg()(10)[c(3,6,9)]
    }else{
        npg = ggsci_col_scheme
    }
    if (is.null(all_genes)){
        all_genes = readRDS('/GCI/jup/A_TSF/dev/ref/all_genes_symbol.RDS')
        print('all_genes is NULL, Read 19039 gene symbols from /GCI/jup/A_TSF/dev/ref/all_genes_symbol.RDS -> all_genes')
    }
    if (!is.character(sig_genes)){
        print('Try to Convert sig_genes into charactor')
        sig_genes = as.character(sig_genes)
    }

    figureTitle = paste("DiffExpr GO enrichment",ont)
    Name = 'GO_enrich'

    # 生成symbol-> entrezID的map：GENE.REF
    tryCatch(
        expr = {
            GENE.REF = get("GENE.REF", envir = .GlobalEnv)
            stopifnot(is.data.frame(GENE.REF))
            stopifnot(nrow(GENE.REF)>0)
            print('Got Global var GENE.REF')
        },
        error = function(err){
            assign(
                "GENE.REF",
                value = bitr(
                    all_genes,
                    fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db
                )  %>%
                    arrange(SYMBOL) %>%
                    distinct(SYMBOL, .keep_all = TRUE) %>%
                    column_to_rownames('SYMBOL')
                ,
                envir = .GlobalEnv
            )
            GENE.REF = get("GENE.REF", envir = .GlobalEnv)
            print('Global var GENE.REF assigned')
        }
    )

    rm.genes = setdiff(sig_genes,rownames(GENE.REF))

    if (length(rm.genes)>0){
        sig_genes = intersect(sig_genes, rownames(GENE.REF))
        print(paste0('remove ',length(rm.genes),' genes symbol not in GENE.REFß:'))
        print(rm.genes)
    }
    # sig_genes/diff为GO分析的主体数据, 使用ENTREZID做分析, 使用全局变量GENE.REF作为map转换symbols->ENTREZID
    ego = enrichGO(gene       = GENE.REF[sig_genes,'ENTREZID'],
                   universe      = GENE.REF$ENTREZID,
                   ont           = ont,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = p.t,
                   qvalueCutoff  = p.adjust.t,
                   readable      = TRUE,
                   ...
    )
    if (nrow(ego@result)==0){
        print('GO analysis result in 0 terms enriched')
        res.tb = data.frame()
    }else{
        res.tb = ego@result %>% mutate(
            Description = fct_reorder(Description,p.adjust,.desc = TRUE),
            ONTOLOGY = factor(ONTOLOGY,levels=c('BP','MF','CC')),
            EnrichedGeneNum = GeneRatio %>% str_split('/',simplify = TRUE) %>% `[`(,1) %>% as.integer,
            setSize = BgRatio %>% str_split('/',simplify = TRUE) %>% `[`(,1) %>% as.integer,
            EnrichedPer = round(EnrichedGeneNum / setSize * 100)
        )
    }

    res = list(gsea.set = Name
               , gene_list = sig_genes
               , res.tb = res.tb
               , res.go = ego
               , args = list(ont=ont
                             ,by.expr.direction=FALSE
                             ,p.t = p.t
                             ,p.adjust.t = p.adjust.t
                             ,log2fc.abs.t = NA
               )
    )
    attr(res,'class') = c(class(res),'y.GO.res')
    message('GO analysis done, now you can use yplot_GO_res to get the figure')

    return(res)
}



yplot_GO_res = function(go.res, go.tb.filter='top10',.T_ylab_wrap = 50,ggsci_col_scheme = NULL,...){
    "@param go.tb.filter ='top10' value in c(NULL,'topN',function), if 'topN' select the top N small p.adj value of each Ontology, else a filter function which will be applied to go.res.tb before plotting, if NULL do not apply filter on this level"
    if ('y.GO.res' %in% class(go.res)){
        res.tb = go.res$res.tb
        ont = go.res$args$ont
        figureTitle = paste("DiffExpr GO:", ont)
        by.expr.direction = go.res$args$by.expr.direction

    }else if ('y.GO.res.dtable' %in% class(go.res)){
        figureTitle = 'GO'
        ont =
            res.tb = go.res
        ont = res.tb$ONTOLOGY %>% unique
        ont = ifelse( ont %>% length > 1,'ALL', ont)
        by.expr.direction = ifelse(res.tb$Change[[1]] == 'DiffExpr', FALSE, TRUE)
    }else{
        stop('wrong go.res type, must be y.GO.res or y.GO.res.dtable')
    }

    if (is.null(ggsci_col_scheme)) {
        npg = pal_npg()(10)
    }
    else {
        npg = ggsci_col_scheme
    }

    if (!is.null(go.tb.filter)){
        if (is.character(go.tb.filter)){
            if (str_sub(go.tb.filter,0,3)=='top'){
                print(paste('using filter',go.tb.filter))
                top = go.tb.filter %>% str_sub(4) %>% as.integer
                data.tb = res.tb %>% group_by(ONTOLOGY) %>% top_n(-top, wt = pvalue) %>% ungroup
            }else{
                warning("go.tb.filter<str> value is wrong, No filtered applied")
                data.tb = res.tb
            }
        }else if(is.function(go.tb.filter)){
            print('using custom filter funcion')
            data.tb = go.tb.filter(res.tb)
        }else{
            data.tb = res.tb
        }
    }
    n.res = data.tb$ID %>% unique %>% length
    print(paste0("[ ", n.res, " ] GO:", ont, " results left after filtering"))

    tab = data.tb %>% dplyr::select(Description, EnrichedGeneNum,
                                    EnrichedPer)
    maxlen_ylabel = data.tb$Description %>% sapply(str_length) %>%
        max
    if (maxlen_ylabel > .T_ylab_wrap) {
        ylabel_size = 10
        height_ratio = 1 + 0.15 * (maxlen_ylabel%/%.T_ylab_wrap)
        maxlen_ylabel = .T_ylab_wrap
    }
    else {
        ylabel_size = 14
        height_ratio = 1
    }
    w = ceiling((maxlen_ylabel/10 + 7) * 2)/2
    h = ((round(min(30, nrow(data.tb))/3.5 + 3) * height_ratio - 2.5) *
             2)/2
    gg.go = data.tb %>% ggplot(aes(x = -log10(pvalue),y = Description)) +
        facet_grid(ONTOLOGY~.,scales = 'free', space = 'free') +
        ylab(NULL) +
        ggtitle(figureTitle) +
        theme_bw(base_size = 16) +
        theme(axis.text.y = element_text(size = ylabel_size, lineheight = 0.5)) +
        scale_y_discrete(breaks=tab$Description,
                         labels = function(x)str_wrap(x,width = maxlen_ylabel)) +

        #
        guides(
            y.sec = ggh4x::guide_axis_manual(
                breaks = tab$Description,
                labels = paste0(tab$EnrichedGeneNum, "(", tab$EnrichedPer,
                                "%)")
            ),
            size = guide_legend(title = "Enriched%"),

        )
    if (by.expr.direction==TRUE){
        gg.go = gg.go +
            geom_point(aes(fill = Change,  shape = ONTOLOGY,
                           size = EnrichedPer),alpha = 0.65) +
            scale_shape_manual(breaks=c('BP',"MF",'CC'),values = c(21,22,24)) +
            scale_fill_manual(breaks=c('Up','Down'), values=npg[1:2]) +
            guides(shape = guide_legend(title='ONTOLOGY',override.aes = list(size=4))
                   ,fill = guide_legend(title = "Direction",override.aes = list(size=4,shape=24))
            )
    }else{
        gg.go = gg.go +
            geom_point(aes(fill=ONTOLOGY, size = EnrichedPer),alpha = 0.65,shape=21) +
            #scale_shape_manual(breaks=c('BP',"MF",'CC'),values = c(21,22,24)) +
            scale_fill_manual(breaks=c('BP',"MF",'CC'),values = npg[3:5])+
            guides(fill = guide_legend(title = "ONTOLOGY",override.aes = list(size=4,shape=21)))
    }
    make.custom(w, h)
    figure_size = c(w, h)
    res = list(plot.tb = data.tb
               , gg = gg.go
               , figsize = figure_size
               , w = w
               , h = h
               , height_ratio = height_ratio
               , maxlen_ylabel = maxlen_ylabel
    )
    attr(res,'class') = c(class(res),'y.GO.res.plot')
    res

}

# ydo_GO = function(sig_genes
#                   ,all_genes =NULL
#                   ,ont = 'ALL'
#                   ,p.t = 0.01
#                   ,p.adjust.t = 0.05
#                   ,ggsci_col_scheme = NULL
#                   ,...
# ){

#     .T_ylab_wrap = 50
#     ont = toupper(ont)
#     if (is.null(ggsci_col_scheme)){
#         npg = pal_npg()(10)[c(3,6,9)]
#     }else{
#         npg = ggsci_col_scheme
#     }
#     if (is.null(all_genes)){
#         all_genes = readRDS('/GCI/jup/A_TSF/dev/ref/all_genes_symbol.RDS')
#         print('all_genes is NULL, Read 19039 gene symbols from /GCI/jup/A_TSF/dev/ref/all_genes_symbol.RDS -> all_genes')
#     }
#     if (!is.character(sig_genes)){
#         print('Try to Convert sig_genes into charactor')
#         sig_genes = as.character(sig_genes)
#     }

#     figureTitle = paste("DiffExpr GO enrichment",ont)
#     Name = 'GO_enrich'

#     # 生成symbol-> entrezID的map：GENE.REF
#     tryCatch(
#         expr = {
#             GENE.REF = get("GENE.REF", envir = .GlobalEnv)
#             print('Got Global var GENE.REF')
#         },
#         error = function(err){
#             assign(
#                 "GENE.REF",
#                 value = bitr(
#                     all_genes,
#                     fromType = "SYMBOL",
#                     toType = c("ENTREZID"),
#                     OrgDb = org.Hs.eg.db
#                 )  %>%
#                     arrange(SYMBOL) %>%
#                     distinct(SYMBOL, .keep_all = TRUE) %>%
#                     column_to_rownames('SYMBOL')
#                 ,
#                 envir = .GlobalEnv
#             )
#             print('Global var GENE.REF assigned')
#         }
#     )

#     rm.genes = setdiff(sig_genes,rownames(GENE.REF))
#     if (length(rm.genes)>0){
#         sig_genes = intersect(sig_genes, rownames(GENE.REF))
#         print(paste0('remove ',length(rm.genes),' symbols in sig_genes:'))
#         print(rm.genes)
#     }
#     # sig_genes/diff为GO分析的主体数据, 使用ENTREZID做分析, 使用全局变量GENE.REF作为map转换symbols->ENTREZID
#     ego = enrichGO(gene       = GENE.REF[sig_genes,'ENTREZID'],
#                    universe      = GENE.REF$ENTREZID,
#                    ont           = ont,
#                    OrgDb         = org.Hs.eg.db,
#                    keyType       = "ENTREZID",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = p.t,
#                    qvalueCutoff  = p.adjust.t,
#                    readable      = TRUE,
#                    ...
#     )

#     res.tb = ego@result %>% mutate(
#         Description = fct_reorder(Description,p.adjust,.desc = TRUE),
#         ONTOLOGY = factor(ONTOLOGY,levels=c('BP','MF','CC')),
#         EnrichedGeneNum = GeneRatio %>% str_split('/',simplify = TRUE) %>% `[`(,1) %>% as.integer,
#         setSize = BgRatio %>% str_split('/',simplify = TRUE) %>% `[`(,1) %>% as.integer,
#         EnrichedPer = round(EnrichedGeneNum / setSize * 100)
#     )


#     n.res = res.tb %>% nrow
#     print(paste0('[ ',n.res,' ] GO:',ont,' results left after filtering'))

#     tab = res.tb %>% dplyr::select(Description,EnrichedGeneNum,EnrichedPer)

#     maxlen_ylabel = res.tb$Description %>% sapply(str_length) %>% max
#     if (maxlen_ylabel > .T_ylab_wrap){
#         ylabel_size = 10
#         height_ratio = 1 + 0.15 * (maxlen_ylabel %/% .T_ylab_wrap)
#         maxlen_ylabel = .T_ylab_wrap
#     }else{
#         ylabel_size = 14
#         height_ratio = 1
#     }
#     w = ceiling((maxlen_ylabel / 10 + 7)*2) / 2 # to 0.5的倍数
#     h = ((round(min(30,nrow(res.tb)) / 3.5 + 3) * height_ratio)*2) / 2
#     gg.go = res.tb %>% ggplot(aes(x = -log10(pvalue), y = Description, color=ONTOLOGY)) +
#         geom_point(aes(size = EnrichedPer),alpha=0.8) +
#         #         scale_colour_gradient(limits=c(1,50), low="red") +
#         # scale_shape_manual(values = c(15,16,17,18)) +
#         scale_color_manual(values = npg) +
#         ylab(NULL) +
#         ggtitle(figureTitle) +
#         theme_bw(base_size = 16) +
#         theme(axis.text.y = element_text(size = ylabel_size, lineheight = 0.6)) +
#         scale_y_discrete(labels = function(x) str_wrap(x, width = maxlen_ylabel)) +
#         guides(y.sec = ggh4x::guide_axis_manual(
#             breaks= tab$Description
#             ,labels = paste0(tab$EnrichedGeneNum,'(',tab$EnrichedPer,'%)')
#         ),
#         size = guide_legend(title = 'Enriched%'),
#         color = guide_legend(title = 'Ontology')
#         )

#     make.custom(w,h)
#     figure_size = c(w,h)
#     list(gsea.set=Name
#          ,gene_list=sig_genes
#          ,res.go=ego
#          ,tb.go=res.tb
#          ,gg=gg.go
#          ,figsize=figure_size
#          ,w=w
#          ,h=h
#          ,height_ratio=height_ratio
#          ,maxlen_ylabel=maxlen_ylabel)
# }
# else if (key == ){
#           v = ref[[key]]
#             path = file.path(base_ref_dir,v)
#             original_gmt_GSVA = readLines(path)
#             strsplit_no_name = function(gmt.list_layer){ as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))[-2] }
#             ref_df = lapply(original_gmt_GSVA, strsplit_no_name)
#             for (layers in 1:length(ref_df)) {
#             names(ref_df)[layers] = ref_df[layers][[1]][1]
#                 ref_df[layers][[1]] = ref_df[layers][[1]][-1]
#             }
#         }




yread_table = function(table,frm=INPUTROOT,...){
    # table is a dataframe,matrix or a string represent the df,matrix file name
    if ((table %>% class == 'character')&&(length(table)==1)) {
        table2 = yload_dfx(li=table,frm=frm,...)
        return(table2)
    }else if (length(intersect(c('data.frame','matrix'),table %>% class) ) > 0){
        # table is data.frame or tibble or matrix
        return(table)
    }else{
        stop('TYPE error of table')
    }
}

#' like make.names but use '.' instead of '.'
#'
#' @description same as make.names, parameters doc see make.names details
#'
#' @param names character, strings to be format into legal symbols
#' @param unique
#'
#' @return character vector
#' @export
#'
#' @examples
ymake_names = function(names, unique = FALSE){rr
    names %>% make.names(unique = unique, allow_ = TRUE) %>% str_replace(fixed('.'),fixed('_'))
}


library(RColorBrewer)
yroll_colors_n = function(n,N=7,random=FALSE){
    library(RColorBrewer)
    # N = 7
    if (n <= N){
        # sample(brewer.pal(7,"Set2"),n)
        if (random==TRUE) c = sample(brewer.pal(8,"Set1"),n)
        else c =brewer.pal(7,"Set1")[1:n]
    }else if (n >N){
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

        if (random==TRUE) c = sample(col_vector,n)
        else c =col_vector[1:n]
    }
    c
}


yloadz_params = function(keys=NULL,.json='.ex.parameters.json',.frm='./',retmode='local',verbose=TRUE){
    # use package jsonlite to load frm/json file, assign its containing dict keys to .GlobalEnv with values
    # at current, 2022年2月19日, only support un-nested json dict
    # @params retmode ['local','global']
    library(jsonlite)
    log = c()
    i = 0
    json_file = yfile_path(.frm,.json)
    if (file.exists(json_file)){
        parameters = jsonlite::fromJSON(json_file,simplifyVector = T)
        if ('list' %in% class(parameters)){
            row = parameters
        }
        if ('data.frame' %in% class(parameters)){
            row = parameters[1,] %>% list
        }

        if (keys %>% is.null == FALSE) pars = ypop(row,keys)
        else { pars = row; row = NULL }

        if (retmode=='global'){
            for (var in names(pars)){
                log = c(log,paste0("  ",var," <- '",pars[var],"' "))
                assign(var,pars[[var]],.GlobalEnv)
                i = i + 1
            }
        }else if (retmode=='local'){
            i = pars %>% length
        }
        else stop('wrong <retmode> value')

        row %>% yexport_params(verbose=F)
        ylog("INFO ",json_file," loaded. [",i,'] vars Parsed ',log %>% str_flatten,verbose=verbose)
        if (retmode=='local') return(pars)
        NULL
    }else{
        # nothing to do
        ylog('Config_file <',json_file,'> does not exist, Nothing parsed')
    }
}

##@
yexport_params = function(li,.json = '.ex.parameters.json',.outputdir = './',verbose=TRUE){
    out = yfile_path(.outputdir,.json)
    jsonlite::write_json(li,out,auto_unbox=TRUE)
    ylog('params exported at ',out,verbose=verbose)
}


yloadz_config = function(keys=NULL,.json='config.json',.frm='./',retmode='local',verbose=TRUE){
    # use package jsonlite to load frm/json file, assign its containing dict keys to .GlobalEnv with values
    # at current, 2022年2月19日, only support un-nested json dict
    # @params retmode ['local','global']
    library(jsonlite)
    log = c()
    i = 0
    json_file = yfile_path(.frm,.json)
    if (file.exists(json_file)){
        parameters = jsonlite::fromJSON(json_file)
        if ('list' %in% class(parameters)){
            row = parameters
        }
        if ('data.frame' %in% class(parameters)){
            row = parameters[1,] %>% list
        }

        if (keys %>% is.null == FALSE) pars = ygetlast(row,keys)
        else { pars = row }

        if (retmode=='global'){
            for (var in names(pars)){
                log = c(log,paste0("  ",var," <- '",pars[var],"' "))
                assign(var,pars[[var]],.GlobalEnv)
                i = i + 1
            }
        }else if (retmode=='local'){
            i = pars %>% length
        }
        else stop('wrong <retmode> value')
        ylog("INFO ",json_file," loaded. [",i,'] vars Parsed ',log %>% str_flatten,verbose=verbose)
        if (retmode=='local') return(pars)
        NULL
    }else{
        # nothing to do
        ylog('Config_file <',json_file,'> does not exist, Nothing parsed')
    }
}







library("ggsci")


library("RColorBrewer")


# `+` = function(a,b){
#     if (a %>% is.character && b %>% is.character){
#         return (paste0(a,b))
#     }else{
#         return(.Primitive("+")(a,b))
#     }
# }



#' ypsig_mark
#'
#' @description convert a p value into sig marks like '***'
#'  see also ymark_psig
#'
#' @param p numeric vector of length 1, represent p values of statistics, should 0 < p < 1
#'  if p value is not numeric, 'N/A' will be returned at that position
#'
#' @return character vector
#'  * -> p < 0.05
#'  ** -> p < 0.01
#'  *** -> p < 0.001
#'  **** -> p < 0.0001
#'  n.s. -> p > 0.05 or is.na(p)
#'  N/A -> not a number, can't apply mark sig
#' @export
#'
#' @examples
#' list(0.2,0.015,0.0005,NA,'a') %>% apply(ypsig_mark) # c('n.s.','*','***','n.s.','N/A')
ymark_psig = function(p){
    if (p %>% is.na){
        return('n.s.')
    }else if (p %>% is.numeric){
        if (p < 0.0001){
            return('****')
        }else if (p < 0.001){
            return('***')
        }else if (p < 0.01){
            return("**")
        }else if (p < 0.05){
            return("*")
        }else{
            return('n.s.')
        }
    }else{
        return('N/A')
    }
}

# ymark_psig = ypsig_mark

#' format p value into sig digits
#'
#' @param p numeric charactor, represents a p value
#' @param digits how many digits to keep
#' @param .ret.prefix ='p = ' control the returned string prefix of formatted p values
#'  if .ret.prefix is not a character, return numeric fortmatted p
#'
#' @return numeric / charactor depending on .ret.prefix value
#'  if p is not numeric, return 'N/A'
#' @export
#'
#' @examples
ypsig_format = function(p, digits = 3, .ret.prefix = 'p = ') {
    if (p %>% is.numeric) {
        if (is.character(.ret.prefix)) {
            if (p %>% is.na) {
                return('n.s.')
            } else if (p < 0.05) {
                return(paste(.ret.prefix, signif(p, digits)))
            } else{
                return('n.s.')
            }
        } else{
            if (p %>% is.na) {
                return(NA)
            } else {
                return(signif(p, digits))
            }
        }
    } else{
        return('N/A')
    }
}

ystat_chi2.every.row = function(data,y='cnt',x = 'feature2',group='InDel',...,.retdetails=FALSE,.hint=TRUE){
    "y ~ x , groupby group"
    "x is levels, y is values(numeric) to be calculated, group %>% unique will be 2x2 table column names"
    "returns 2x2 cross table chi-square test result of every level"
    res = list()
    ps = c()
    methods = c()
    grp1 = c()
    grp2 = c()
    levels = data[[x]]%>% unique
    mat = data[,c(y,x,group)] %>% pivot_wider(names_from = group,values_from = y) %>% column_to_rownames(x)
    rn = rownames(mat)
    for (i in levels){
        mat2 = rbind(
            mat[rn==i,],
            mat[rn!=i,] %>% summarise(across(.fns = sum)) %>% `rownames<-`('Others')
        )
        stat = chisq.test(mat2,...)
        ps %<>% c(stat$p.value)
        methods %<>% c(stat$method)
        grp1 %<>% c(stat$observed[1,1])
        grp2 %<>% c(stat$observed[1,2])
        res[[i]] = stat
    }
    res[['FDR']] = p.adjust(ps)
    res[['pvalues']] = ps
    #     print('FDR is')
    #     print(paste(levels,'=',res[['FDR']]))
    p.table =data.frame(levels=as.character(levels),grp1=grp1,grp2=grp2,pval=ps,FDR=res[['FDR']],method=methods) %>% dplyr::arrange(FDR)
    .sum = mat %>% summarise(across(.fns = sum)) %>% c
    grps = colnames(mat)
    p.table %<>%
        rename(!!grps[[1]]:=grp1,!!grps[[2]]:=grp2) %>%
        mutate(sig.pval = p.table$pval %>% sapply(ymark_psig),.after = pval) %>%
        mutate(sig.FDR = p.table$FDR %>% sapply(ymark_psig),.after = FDR) %>%
        rbind(list('Total',.sum[[1]],.sum[[2]],NA,'-',NA,'-','-'))
    if (.hint==TRUE){
        print('each row calculate the p value of a 2x2 cross table of [level_grp1,level_grp2,the_other_grp1,the_other_grp2] colnames meaning: levels, grp1_oberved, grp2_oberved, pval, sig.pval, FDR, sig.FDR, static method used')
    }
    if (.retdetails==TRUE){
        return(list(p=p.table,stat.info=res))
    }else{
        return(p.table)
    }
}



