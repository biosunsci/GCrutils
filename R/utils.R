

#'
#' Log messages
#'
#' Log messages with timestamps to the stdout or stderr
#'
#' @param ... parameters of string, which will be concat to print to
#' @param verbose if verbose is FALSE, stop logging anything
#' @param .addtime add timestamp at the front when print
#' @param .dev in [1,2] default 2(stderr), .dev==2 mean stderr, .dev==1 means stdout
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' file = 'abc.txt'
#' dir = '/path/to/save'
#' ylog(file,'is written at',dir)
#' }
#'
ylog = function(...,verbose=TRUE, .addtime=TRUE, .dev=2){
    if (verbose==TRUE){
        if (.addtime){
            msg = c(list('[',Sys.time() |> as.character(),'] '),list(...)) %>%
                unlist %>%
                as.character %>%
                str_flatten(collapse = ' ')
        }else{
            msg = list(...) %>% unlist %>% as.character %>% str_flatten(collapse = ' ')
        }
        print(msg)
        # message(msg)
    }
}

#' check string is a valid R var/obj name
#'
#' check string is a valid R var/obj name, return TRUE if string is valid False if it is not
#'
#' @param string
#'
#' @return TRUE, if string is a valid else FALSE
#' @export
#'
#' @examples
#' yis_valid_unreserved('.jjj') # returns TRUE
#' yis_valid_unreserved('_jjj') # returns FALSE
yis_valid_unreserved <- function(string) {
    make.names(string) == string
}

# -------------- list utils ----------------

#' ypush up-sert (upsert) items to list
#'
#' Unlike append, ypush is an update-append(like SQL upsert) version of append
#'
#' @param ...
#' @param .expand
#' @param .as.whole
#' @param iter an iterable obj like list
#'
#' @return updated `iter`
#' @export
#'
#' @examples
#' \dontrun{
#' list(a=1,b=2) %>% append(list(a=3,c=4)) # list(a=1,b=2,a=3,c=4)
#' list(a=1,b=2) %>% ypush(list(a=3,c=4)) # list(a=3,b=2,c=4)
#' list(a=1,b=2) %>% ypush(list(a=3,c=4),.expand=F) # list(a=3,b=2,c=4)
#' }
ypush = function(iter, ..., .expand=TRUE,.as.whole=FALSE) {
    # update on keys for named list, push for unnamed list
    # you must let var = ypush(var,...) to update var, for the convenience of do.call
    # if partial named list, all no name items is ignored
    # if (iter %>% class =='character') iter = get(iter ,envir = parent.env(environment()))
    items = list(...)
    if (.expand==TRUE && length(items) == 1 && names(items) %>% is.null && 'list'%in%class(items[[1]])){
        items = items[[1]]
    }
    un_named = names(items) %>% is.null()
    if (un_named){
        if (.as.whole == FALSE){
            return(c(iter,items))
        }else{
            # .as.whole == TRUE
            iter[[length(iter)+1]] = items
            return(iter)
        }
    }else{
        # items is named list
        if (.as.whole == FALSE){
            ns = names(items)
            for (i in 1:length(items)){
                n = ns[[i]]
                v = items[[i]]
                if (n == ""){
                    if (is.null(v)){
                        iter[length(iter)+1] = list(NULL)
                    }else{
                        iter[[length(iter)+1]] = v
                    }
                }else{
                    if (is.null(v)){
                        iter[n] = list(NULL)
                    }else{
                        iter[[n]] = v
                    }
                }
            }
        }else{
            # .as.whole == TRUE
            iter[[length(iter)+1]] = items
            return(iter)
        }
        return(iter)
    }
}


#' apply func over a list x
#'
#' @param x the list to be itered
#' @param func apply the func over x items, func could be 1,2,3 parameters
#' 1: func(value)
#' 2: func(value,key)
#' 3: func(value,key,index)
#'
#' @return None
#' @export
#'
#' @examples
yloop = function(x,func){
    l = length(formals(func))
    if (l==3) for (i in 1:length(x)){
        func(x[[i]],names(x[i]),i)
    }else if (l==2)  for (i in 1:length(x)){
        func(x[[i]],names(x[i]))
    }else if (l==1)  for (i in 1:length(x)){
        func(x[[i]])
    }
}



#
#' Get last element(s)
#'
#' Get the last element which match the keys, of an iterable, like ypop, but do not change .obj,
#' *** Usually used for screening out unused arguments passed by ... in a function ***
#' Used by .ydumpto.function
#'
#' Caution: list(1,3,4) |> names() -> NULL
#' list(1,3,a=4) |> names() -> c('','','a')
#'
#' @param .obj the iterable, list
#' @param keys keys must be a vector of keys used for searching in names(.obj<list>)
#' @param .squeeze boolean,  .squeeze = T equals .retmode == 'squeeze'
#' @param .retmode in c('list','squeeze'), 'list' is default
#'
#' @return
#' @export
#'
#' @examples
ygetlast = function(.obj, keys=NULL, .squeeze = FALSE, .retmode='list') {
    # keys is NULL,c(),list(), but not NA for length(NA) == 1
    # key not provided, return the last element of the .obj
    if (length(keys)==0) {
        out = .obj[[length(.obj)]]
        if (.retmode == 'list'){
            out = list(out)
        }
        return(out)
    }
    # else
    if (.squeeze == TRUE) .retmode = 'squeeze'
    stopifnot(keys |> is.vector())
    out = list()
    EmptystringKEY_FLAG = FALSE
    for (name in keys){
        # if name == '',
        if (name=="" && EmptystringKEY_FLAG==FALSE){
            # get all unamed values of a partial named list
            # like list(1,2,a=4,6) by key='' will result list(1,2,6)
            # stop(' get partial named list  by key = "" is not supported ')
            tmp = .obj[names(.obj)==""]
            out = append(out,tmp)
            EmptystringKEY_FLAG = TRUE
        }else{
            r = .obj[[name]]
            if (!is.null(r)){
                if (is.double(name) || is.integer(name)){
                    out[[length(out)+1]] = r
                }else{
                    out[[name]] = r
                }
            }
        }
    }
    if (.retmode=='squeeze' && length(out) == 1){
        out = out[[1]]
    }
    return(out)
}

#' Remove the last matched by keys
#'
#' yrmlast is the counterpart of the ygetlast, details See examples
#'
#' @param .obj is a named or unamed list or vector
#' @param keys is a unnamed vector
#' @param .byindex boolean
#'
#' @return the .obj with keys removed
#'
#' @examples
#' \dontrun{ x = list(a=4,b=3,c=2,d=1) # if .byindex == FALSE (default) # if .obj is named list or
#' vector, removal depents on keys, and `keys` is treated as names in .obj # if .obj is unnamed list
#' or vector, removal depents on values, and `keys` is treated as values in .obj x %>%
#' ygetlast(c('a','b','e')) # list(a=4,b=3) x %>% yrmlast(c('a','b','e')) # list(c=2,d=1) # if
#' `keys` is NULL,c(),list(), but not NA, remove nothing, `keys` must be an unnamed vector or list #
#' if .byindex == TRUE, keys must be positive numbers # representing the index of the .obj ranging
#' at 1~length(.obj) x %>% yrmlast(c(1,2),.byindex = T) # list(c=2,d=1) }
yrmlast = function(.obj, keys = NULL, .byindex=FALSE){
    stopifnot((keys |> is.vector())&&(names(keys) |> is.null()))
    # keys is NULL,c(),list(), but not NA for length(NA) == 1
    if (length(keys)==0){
        return(.obj)
    }
    judge_unamed = .obj |> names() |> is.null()
    if (.byindex==FALSE){
        if (judge_unamed==TRUE){
            return (.obj[!(.obj %in% keys)])
        }else{
            mask = names(.obj) %in% keys
            return(.obj[!mask])
        }
    }else{
        # remove by index
        stopifnot( (keys |> is.numeric()) && all(keys > 0) )
        return(.obj[keys*-1])
    }
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





# ----------- df utils ----------------
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
ygen_subMafs_ = function(maf,glx,col = NULL,id_col='Tumor_Sample_Barcode',.container=NULL,.assignGlobal=TRUE,...){
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

