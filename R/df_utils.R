#' Check whether a DataFrame has true rownames
#'
#' @details check whether a data.frame has rownames or not TRUE if there is a rownames vector not like 1:N (The natural one) else FALSE 一般用来确定diff expr等data.frame是否有symbol作为rownames，还是没有rownames，symbol为data.frame的一列
#'
#' @param df must be an obj whose class is data.frame
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
yhas_rownames = function(df){
    if(is.data.frame(df)){
        rn = .row_names_info(df)
        if (rn <= 0) return(FALSE)
        else return(TRUE)
    }else{
        stop('df must be data.frame type')
    }
}



#' Check DataFrame rownames and column names
#'
#' @description like [yhas_rownames()] check whether a data.frame has rownames or not TRUE if there is a rownames vector else FALSE
#'
#' @param df must be a class of data.frame
#'
#' @return TRUE or FASLE
#' @export
#'
#' @examples
yhas_colnames = function(df){
    stopifnot(is.data.frame(df))
    cn = Negate(is.null)(names(df))
    cn
}


#' Convert column named data.frame to named list
#'
#' Convert a data.frame into list in length of ncol(df), with names(list) is the df's colnames
#' and each
#'
#' @param df a data.frame object
#'
#' @return list
#' @export
#'
#' @examples
ydeframe = function(df){
    rn = rownames(df)
    res = list()
    for (col in colnames(df)){
        res[[col]] = df[,col]
        names(res[[col]]) <- rn
    }
    res
}




#' display a dataframe/matrix table with more rows and columns displayed
#'
#' @description used with Jupyter notebook display of data table, Need View function available at Global.Env
#'
#' @param mat the dataframe/matrix to be displayed
#' @param nrowmax the max number of rows to display
#' @param ncolmax the max number of columns to display
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' \dontrun{
#'  data(mtcars)
#'  echo(mtcars)
#' }
#'
#'
echo = function(mat,nrowmax = 100,ncolmax=100){
    cur_nrow = getOption("repr.matrix.max.rows")
    cur_ncol = getOption("repr.matrix.max.cols")
    col = min(ncolmax,ncol(mat))
    row = min(nrowmax,nrow(mat))
    options(repr.matrix.max.cols=max(30,col))
    options(repr.matrix.max.rows=max(10,row))
    utils::View(mat)
    options(repr.matrix.max.cols=cur_ncol)
    options(repr.matrix.max.rows=cur_nrow)
}


#' calculate Z-score of each row
#'
#' Z-score = (x-mean) / sd, calculate Z-score of each row
#'
#' @param x matrix of number
#'
#' @return row Z-scored matrix
#' @export
#'
#' @examples
#'
yscale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, stats::sd, na.rm = T)
    return((x - m) / s)
}



#' subset named list by keys
#'
#' @description subset the list x by names(keys) if keys is named list or vector
#' subset the list x by keys if keys is not named list or vector
#'
#' @param x is.list(x) must be TRUE
#' @param keys vector or list, return the element in x which match the keys, others in the keys are ignored.
#'  if names(keys) is null, use values of keys
#'  else use names(keys) to subset
#'
#' @return list with only names in intersect(names(x),~keys) is returned
#' @export
#'
#' @examples
#' a = list(a=1,b=2,c=3)
#' k1 = c('a','c','d')
#' k2 = c(a="x",b=9,d='qq')
#' k3 = list(a='x',x='x2',c='ff')
#' ysubset_list_named(a,k1) # list(a=1,c=3)
#' ysubset_list_named(a,k2) # list(a=1,b=2)
#' ysubset_list_named(a,k3) # list(a=1,c=3)
#'
ysubset_list_named = function(x,keys){
    stopifnot(is.list(x))
    xnames = names(x)
    knames = names(keys)
    if (is.null(knames)){
        knames = keys
    }
    nm = intersect(xnames,knames)
    return(x[nm])
}

#' SubsetMaf maftools::maf object
#'
#' @param maf maftools::maf object
#' @param glx group data used to subset the maf obj
#' @param col group column name used to subset the maf
#' @param id_col id_col used to link glx and maf Tumor_Sample_Barcode field
#' @param .container
#' @param .assignGlobal
#' @param ...
#'
#' @return if .assignGlobal==FALSE, return list else return NULL
#' @export
#'
#' @examples
yget_subset_maf = function(maf,glx,col = NULL,id_col='Tumor_Sample_Barcode',.container=NULL,.assignGlobal=TRUE,...){
    if (!is.null(.container)){
        .assignGlobal = FALSE
        stopifnot(is.list(.container))
    }
    if (class(maf)=='MAF'){
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
            mask = maf[[id_col]] %in% sub.glx[[id_col]]
            sub.maf = maf[mask,]
            submafs[[g]] = read.maf(sub.maf,sub.glx,...)
        }
        return (submafs)
    }else{
        g1name = grps[[1]]
        sub.glx1 = glx[glx[[col]] == g1name,]
        sub.maf1 = maf[maf[[id_col]] %in% sub.glx1[[id_col]],]
        g1 = read.maf(sub.maf1,sub.glx1,...)
        g2name = grps[[2]]
        sub.glx2 = glx[glx[[col]] == g2name,]
        sub.maf2 = maf[maf[[id_col]] %in% sub.glx2[[id_col]],]
        g2 = read.maf(sub.maf2,sub.glx2,...)
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
