SSGSEA_BASE_REF_DIR = '/share_storage/pipeline/REF_FILES/'

library(GenomicRanges)
yget_chr_window = function(w = 1e6){
    chr_df = circlize::read.chromInfo()$df
    chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
    chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))

    chr_window = EnrichedHeatmap::makeWindows(chr_gr, w = w)
    chr_window
}

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

#-----------------

yhas_rownames = function(df){
    rn = .row_names_info(df)
    if (rn < 0) return(FALSE)
    else return(TRUE)
}
yhas_colnames = function(df){
    cn = Negate(is.null)(names(df))
    cn
}

ydeframe = function(df){
    rn = rownames(df)
    res = list()
    for (col in colnames(df)){
        res[[col]] = df[,col]
        names(res[[col]]) = rn
    }
    res
}

# yhas_rownames = function(df) {
#   !all(row.names(df)==seq(1, nrow(df)))
# }

#-----------------

# CODE IN CREATE PROJECT
library(tidyverse)
# library(ggfortify)
library(ggsci)

if (!exists("WORKER",.GlobalEnv)){
    WORKER = NULL
}


#' display a dataframe/matrix table with more rows and columns displayed
#'
#' @description used with Jupyter notebook display of data table
#'
#' @param mat the dataframe/matrix to be displayed
#' @param nrowmax the max number of rows to display
#' @param ncolmax the max number of columns to display
#'
#' @return NULL
#' @export
#'
#' @examples
#' data(mtcars)
#' echo(mtcars)
#'
echo = function(mat,nrowmax = 100,ncolmax=100){
    cur_nrow = getOption("repr.matrix.max.rows")
    cur_ncol = getOption("repr.matrix.max.cols")
    col = min(ncolmax,ncol(mat))
    row = min(nrowmax,nrow(mat))
    options(repr.matrix.max.cols=max(30,col))
    options(repr.matrix.max.rows=max(10,row))
    View(mat)
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
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
#-----------------

#' ypath_join
#' join given list of folders, likely the R version of os.path.join
#' may cause some unpredictable return values at certain circumstances
#' @param ... parameters of path to join
#'
#' @return character represent the full path
#' @export
#'
#' @examples
library(magrittr)
library(dplyr)

ypath_join = function(...){
    args = list(...)
    args = as.character(args) %>% str_replace('\\\\','/')
    res = paste(args, sep = "/", collapse = "/")
    res = gsub("/[.]{0,1}/", "/", res)
    #     res = gsub("/[\\\\/]{2}/", "/", res)
    res
}

#' see documentation of ypath_join
#' @export
yfile_path = ypath_join

#' check string is a valid R var/obj name, TRUE if string is valid False if it is not
#'
#' @param string
#'
#' @return TRUE, if string is a valid
#' @export
#'
#' @examples
#' is_valid_unreserved('.jjj') # returns TRUE
#' is_valid_unreserved('_jjj') # returns FALSE
yis_valid_unreserved <- function(string) {
    make.names(string) == string
}

#' yslice
#' slice a iterable like `[:]`, but use a python-like style, negative-n means the last-n-th(python style) but not except the n-th(R style)
#'
#' @param iterable vector or list
#' @param start start from 1, can be negative, if negative n then return the last-n-th one element
#' @param end if set, will slice from start to end, negative mean the last start to last end elements
#' @param .style determine the behavior of the function, if pyhton/py like python slice(start,end), default
#'
#' @return sliced iterable
#' @export
#'
#' @examples
#' a = c(0,1,2,3,'b')
#' a %>% yslice(-1) # 'b'
#' a %>% yslice(-3:-1) # c('2','3','b')
#' a %>% yslice(1:3) # c('0','1','2')
#' a %>% yslice(3) #  c('2','3','b')
yslice = function(iterable, seqs=NULL, .style_negative_index='py'){
    " if .style_negative_index == 'py' (default)
     negative values will be treated as to get the last N
     like python does except index is from 1 to length(@iterable)
     and the end of @seqs is included, too
     if .style_negative_index == 'r' (optional)
     negative values will be treated as to remove those at the positive position
     just like in r code []
    "
    l = length(iterable)
    if ( .style_negative_index=='py' || .style_negative_index=='python'){
        start = seqs[1]
        end = seqs[length(seqs)]
        if (start < 0){
            start = start + l + 1
        }
        if (end < 0){
            end = end + l + 1
        }
        seqs = start:end
        return(iterable[seqs])
    }else if (.style_negative_index=='r' || .style_negative_index=='R'){
        # pass
        return(iterable[seqs])
    }else{
        stop('wrong .style value, must in ["python","py","R","r"]')
    }
}

#' ypush up-sert (upsert) items to list
#'
#' Unlike append, ypush is an update-append version of append
#'
#' @param iter an iterable obj like list
#' @param ...,.expand an unamed parameter of type list or a series of parameters,
#'  items to be upsert into `iter`, items = list(...) when `.expand` == FALSE
#'  items = list(...)[[1]] when `.expand` == TRUE only when ... contains only an unamed parameter
#'  .expand default is TRUE, this is useful when updating argv list for a do.call
#'  set to FALSE to stop the behavior and treat items = list(...) for upsert
#' @param .as.whole, control how to treat items generated from ..., default FALSE
#'  if TRUE items will be treat as a whole part to upsert into the last position of the iter
#'  else FALSE items will be iterated and each element item will be tried to upsert the iter
#'
#' @return updated `iter`
#' @export
#'
#' @examples
#' list(a=1,b=2) %>% append(list(a=3,c=4)) -> list(a=1,b=2,a=3,c=4)
#' list(a=1,b=2) %>% ypush(list(a=3,c=4)) -> list(a=3,b=2,c=4)
#' list(a=1,b=2) %>% ypush(list(a=3,c=4),.expand=F) -> list(a=3,b=2,c=4)

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
#-----------------

##@

#-----------------

#' yrmlast, counterpart of the ygetlast
#'
#' @param .obj is a named or unamed list or vector
#' @param keys is a unnamed vector
#' @param .byindex boolean
#'
#' @return the .obj with keys removed
#' @export nothing
#'
#' @examples
#' counterpart of the ygetlast,
#' if .byindex == FALSE (default)
#' if .obj is named list or vector, removal depents on keys, and @keys is treated as names in .obj
#' if .obj is unnamed list or vector, removal depents on values, and @keys is treated as values in .obj
#' e.g. you have x = list(a=4,b=3,c=2,d=1)
#' x %>% ygetlast(c('a','b','e')) -> list(a=4,b=3)
#' x %>% yrmlast(c('a','b','e')) -> list(c=2,d=1)
#' if @keys is NULL,c(),list(), but not NA, remove nothing, @keys must be an unnamed vector or list
#' if .byindex == TRUE, keys must be positive numbers
#' representing the index of the .obj ranging at 1~length(.obj)
#' x %>% yrmlast(c(1,2),.byindex = T) -> list(c=2,d=1)
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
#-----------------

yloop = function(x,func){
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
yreduce = function(.x,.f,.init,.right=FALSE,.accumulate=FALSE){
    Reduce(.f,.x,init = .init,right=.right, accumulate=.accumulate)
}

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
#-----------------

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
#' @param str
#'
#' @param ext set extention name of
#' @param mode in ['auto','keep','add','replace'] default 'auto', behave, see example
#' @param ... diff depend on str and ext values
#'
#' @example
#' ysplit_file_name('a.b.c.txt',ext=NULL,mode='auto') # c('a.b.c', 'txt')
#' ysplit_file_name('abc',ext='dfx',mode='auto') # c('abc', 'dfx')
#' ysplit_file_name('abc.txt',ext='dfx',mode='auto') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',ext='dfx',mode='auto') # c('a.b.c', 'dfx')
#'
#' ysplit_file_name('a.b.c.txt',ext=NULL,mode='replace') # c('a.b.c')
#' ysplit_file_name('abc',ext='dfx',mode='replace') # c('abc', 'dfx')
#' ysplit_file_name('abc.txt',ext='dfx',mode='replace') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',ext='dfx',mode='replace') # c('a.b.c', 'dfx')
#'
#' ysplit_file_name('a.b.c.txt',ext=NULL,mode='keep') # c('a.b.c', 'txt')
#' ysplit_file_name('abc',ext='dfx',mode='keep') # c('abc')
#' ysplit_file_name('abc.txt',ext='dfx',mode='keep') # c('abc', 'txt')
#' ysplit_file_name('a.b.c.txt',ext='dfx',mode='keep') # c('a.b.c', 'txt')
#'
#' ysplit_file_name('abc',ext='dfx',mode='add') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',ext=NULL,mode='add') # c('a.b.c', 'txt')
#' ysplit_file_name('abc.txt',ext='dfx',mode='add') # c('abc', 'txt', 'dfx')
#' ysplit_file_name('a.b.c.txt',ext='dfx',mode='add') # c('a.b.c', 'txt', 'dfx')

ysplit_file_name = function(str,ext=NULL,mode='auto',...){
    if (!is.null(ext)) {
        # 判断ext首字符是不是点“.” 有点则去掉
        if (ext %>% str_sub(1,1)=='.') ext = ext %>% str_sub(2,)
        .count = str_count(str, fixed('.'))
        if (.count == 0) res = c(str)
        else if(.count >= 1) res = str_split(str, fixed('.'))[[1]]
    }else{
        res = str_split(str, fixed('.'))[[1]]
    }
    l = length(res)
    if (l==1){
        # 'abc'
        file_name = res[[1]]
        ext_name = NULL
    }else if (l==2){
        # 'abc.txt'
        file_name = res[[1]]
        ext_name = res[[2]]
    }else if(l>2){
        # 'a.b.c.txt'
        file_name = res %>% yslice(1:-2) %>% str_flatten(collapse='.')
        ext_name = res %>% yslice(-1)
    }
    if (mode=='auto') {
        if (ext |> is.null()){
            res = c(file_name, ext_name)
        }else{
            res = c(file_name, ext)
        }
    } else if (mode=='keep') res = c(file_name,ext_name)
    else if (mode=='add') res = c(file_name, ext_name, ext)
    else if (mode=='replace') res = c(file_name, ext)
    res
}



##@
#' ylog for logging information to the stdout/stderr
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
        write(msg,.dev)
    }
}


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
#-----------------

# like ypop, but do not change .obj
ygetlast = function(.obj, keys=NULL, .squeeze = FALSE, .retmode='list') {
    " *** Usually used for screening out unused arguments passed by ... in a function ***
     @.retmode in c('list','squeeze'), 'list' is default
     @.squeeze = T equals .retmode == 'squeeze'
     @keys must be a vector of keys
     Caution
     list(1,3,4) |> names() -> NULL
     list(1,3,a=4) |> names() -> c('','','a')
    "
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
#' input args ..., filter all the args by the formals of .f, return all matched,
#' .f is.function() and is the filter function, the returned args list will match the formals of the .f
#' .filter is the alias of .f
#' returns like list(...)[intersect(list(...)%>%names,formals(.f)%>%names)]
yget_args = function (..., .f = NULL, .filter = NULL) {

    if ((.filter %>% is.null %>% not)&&(.f %>% is.null)){
        .f = .filter
    }
    dots = list(...)
    pf <- parent.frame()
    pf = pf %>% as.list %>% modifyList(dots)
    if (.f %>% is.null %>% not) {
        formal_args = formals(.f) %>% names
        pf = pf %>% ysubset_list_named(formal_args)
    }
    pf
}

#' split a string into a list of dirnames
#'
#' if path start with /(root), the 1st element of returned list
#' is '/', multiple '/' like '//' will be treated as single '/'
#'
#' @param path, path in string format of length 1
#'
ysplit_path = function(path){
    path.list = str_split(path,pattern = '/+')[[1]]
    if (str_sub(path,1,1)=='/'){
        path.list[[1]] = '/'
    }
    path.list
}

yfunc_args <- function(func=NULL,orig_values = FALSE) {
    ""
    if (is.null(func)){
        # get formaMls for parent function
        parent_formals <- formals(sys.function(sys.parent(n = 1)))
    }else{
        stopifnot(is.function(func))
        # get formals for func
        parent_formals = formals(func)
    }
    # Get names of implied arguments
    fnames <- names(parent_formals)
    # Remove '...' from list of parameter names if it exists
    fnames <- fnames[fnames != '...']

    if (!is.null(func)) return(fnames)

    # Get currently set values for named variables in the parent frame
    args <- evalq(as.list(environment()), envir = parent.frame())
    # Get the list of variables defined in '...'
    args <- c(args[fnames], evalq(list(...), envir = parent.frame()))

    if(orig_values) {
        # get default values
        defargs <- as.list(parent_formals)
        defargs <- defargs[unlist(lapply(defargs, FUN = function(x) class(x) != "name"))]
        args[names(defargs)] <- defargs
        setargs <- evalq(as.list(match.call())[-1], envir = parent.frame())
        args[names(setargs)] <- setargs
    }
    return(args)
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


##@
#' yload_dfx
#'
#' @description 依赖于 yslice, ypush, ylog, ysplit_path, ysplit_file_name, ypath_join, <stringr>
#'  if `li` && `ext` both NULL and `frm` is not NULL, read all files from `frm`
#'
#' @param li overwritting <pattern> settings
#'  a character vector, if length(<li>)>1 generate a named list
#'  else generate the dataframe read, return is determined further by <retmode>
#'  a character vector whose element <i> refers to dfx names,
#'  if <i> is not started with c('/','~','./'), find <i>.<ext> in folder <frm>
#'  else use <i> as the path, ignore <frm> for the current <i>
#'  if input li is a data.frame it will be returned immediately
#' @param frm
#'  a folder.path, all the data will be read from this path
#' @param pattern
#'  a RegExpr to glob, return matched files, only acts when `li` is NULL
#' @param ext
#'  basicly, <li> members do not include an extention name, ext is used to set the <li> extention name. set ext=NULL to only load exactly the name <li> provided.
#'  ext不会修改已经显式声明的扩展名，如li=c('abx.txt'),ext='dfx' 则读入仍是'abx.txt'而不是'abx.dfx'或'
#' @param retmode , in ['local','global']
#'  local: the dataframe read will be returned, you need to declare a local var to accecpt it.
#'  global: the dataframe read will be directly assign to .GlobalEnv
#'
#' @return df or list of dfs if retmode=='local', assign every df to .GlobalEnv if retmode=='global'
#' @export
#'
#' @examples
#'
yload_dfx = function(li = NULL,
                     frm = INPUTROOT,
                     pattern = NULL,
                     ext = NULL,
                     worker= NULL,
                     retmode = 'local',
                     verbose = TRUE,
                     row.names=NULL,
                     ...
) {
    # retmode in ['local','global']
    # ext in [str,NULL]
    # default parameters with custom one
    std_args = yfunc_args(read.table)
    argv = list(sep = "\t",
                header = T,
                quote = "",
                stringsAsFactors = F,
                comment.char = "",
                na.strings = "",
                row.names=row.names
    ) %>% modifyList(yfunc_args(),keep.null = TRUE) %>% ysubset_list_named(std_args)

    if ('data.frame' %in% (li %>% class)) return(li)
    if (length(li) >= 1 || !is.null(li)) {
        # files [(frm1,f_name1,f_ext1),(frm1,f_name2,f_ext2)...]
        files = list()
        if (is.null(li) && !is.null(pattern)) {
            li = list.files(frm, pattern = pattern) %>% as.list
        } else if (is.null(li)) {
            li = list.files(frm) %>% as.list
        }
        for (i in li) {
            tmp = ysplit_path(i)
            if(length(tmp)==1){
                ffname = tmp[[1]]
                myfrm = frm
            }else if (length(tmp) >= 2){
                ffname = tmp %>% yslice(-1)
                myfrm = tmp %>% yslice(1:-2)
                if (myfrm[[1]] %in% c('/','~','.')){
                    myfrm = myfrm %>% str_flatten(collapse='/')
                }else{
                    myfrm = c(frm,myfrm) %>% str_flatten(collapse='/')
                }
            }else{
                stop("item in @li is wrong")
            }
            # 检查i是否能作为变量名
            if (retmode=='global' && make.names(i) != i) {
                # 不能
                ylog('[WARN] ignore ', i, ' for its not a legal var name',.addtime = F)
            } else{
                tmp = ysplit_file_name(ffname,ext=ext,mode='auto')
                fname = tmp[[1]]
                if (length(tmp) == 1){
                    ext = 'dfx'
                }else if(length(tmp) == 2){
                    ext = tmp[[2]]
                }else{
                    stop('wrong tmp length')
                }
                if (fname %>% str_detect(fixed('.'))){
                    write(paste('[INFO] skip ', i,' for its not right format: multi-dot in filename'),2)
                }else{
                    files = files %>% ypush(c(myfrm,fname,ext))
                }
            }
        }
    }
    res = list()
    for (file in files) {
        # files is list of ('path','file_name','ext_name')
        myfrm = file[[1]]
        vname = file[[2]]
        ext = file[[3]]
        if (vname %>% str_sub(1,3) %>% str_detect('[A-Z][0-9]\\_')) worker=NULL
        if (!is.null(worker)) {
            ffname = str_flatten(c(worker,'_',vname,'.',ext))
            input = ypath_join(myfrm,ffname)
            if (!file.exists(input)){
                ylog('[WARN] prefix with <worker> file not found, retry without <worker> one', .addtime = FALSE)
                ffname = str_flatten(c(vname,'.',ext))
            }
        }else{
            ffname = str_flatten(c(vname,'.',ext))
        }
        input = ypath_join(myfrm,ffname)

        if (verbose == TRUE) ylog("read [", ffname , "] as '", vname, "' from > ",myfrm)

        if (ext == 'gmt'){
            library(GSVA)
            original_gmt_GSVA <- readLines(input)
            strsplit_no_name <- function(gmt.list_layer){ as.character(unlist(strsplit(gmt.list_layer, split = '\t', fixed = T)))[-2] }
            ref_df <- lapply(original_gmt_GSVA, strsplit_no_name)
            for (layers in 1:length(ref_df)) {
                names(ref_df)[layers] <- ref_df[layers][[1]][1]
                ref_df[layers][[1]] <- ref_df[layers][[1]][-1]
            }
            res[[vname]] = ref_df
        }else{
            "ext == {txt,dfx,...} "
            argv$file = input
            res[[vname]] = read.table %>% do.call(args=argv)
        }
        if (retmode == "global") {
            assign(vname, res[[vname]], envir = .GlobalEnv)
        }
    }
    if (retmode == 'local') {
        if (length(res) == 1) {
            return(res[[1]])
        } else{
            return(res)
        }
    } else if(retmode=='global'){

    } else{
        # invalid retmode value
        stop('invalid retmode value')
    }
}
#-----------------
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
#-----------------

# library(maftools)
library(scales)
library(openxlsx)
# color_scheme = scale_fill_aaas
#' save pdf/xlsx
#'
#' @description auto save Robject into corresponding file format
#' output file type is auto determined by the class(fc), currently support classes are
#' 'ggplot', 'pheatmap','Heatmap' -> pdf files
#' 'call' -> eval(fc) -> pdf files
#' 'data.frame', 'matrix' -> xlsx files
#' 'list' of 'data.frame'/'matrix's -> xlsx files with each element in the list -> a sheet in the xlsx
#'
#' @param fc the obj to export
#' @param ...,sep all these parameters will be concat to be the output filename, with sep=`sep`
#' @param table_title_lines a list of title lines, element i will be the title line in the corresponding table fc[[i]], make sure the length of fc and table_title_lines are the same
#' @param outputdir output file into these dir
#' @param family used in cairo_pdf(family=family)
#'
#' @return NULL
#' @export
#'
#' @examples
ysave = function (fc=NULL, ..., table_title_lines=list(), outputdir = OUTPUTROOT, sep = "_",family = "sans" )
{
    t = class(fc)
    l = length(t)
    name = str_flatten(c(...), collapse = sep)
    if (l==1 && t=='NULL'){
        return(c(outname=name,outputdir=outputdir,width=WIDTH,height=HEIGHT,sep=sep))
    }
    if (name == "") {
        name = "default"
    }
    ext = 'pdf'
    if (l==1 && t=='list'){
        ext = 'xlsx'
        fpath = file.path(outputdir, paste0(name , '.' , ext))

        # save list of data.frames into excel
        #         file <- file.path('report',"data_titv.xlsx")
        wb <- createWorkbook()
        #if (names(fc) %>% is.null)
        seqnames <- paste0("Sheet", seq_along(fc)) # or names(datas) if provided
        # fill '' names
        sheetnames = c()
        N = names(fc)
        for (i in seq_along(fc)){
            n = N[[i]]
            if (n==""){
                sheetnames = sheetnames %>% append(seqnames[[i]])
            }else{
                sheetnames = sheetnames %>% append(n) %>% str_sub(1,25)
            }
        }
        sheetnames = sheetnames %>% make.unique

        wb = createWorkbook()
        for (i in seq_along(fc)){
            sheetname = sheetnames[[i]]
            addWorksheet(wb,sheetname)
            tt = class(fc[[i]])
            title_line = table_title_lines[i]
            rowNames = yhas_rownames(fc[[i]])
            if (tt %>% intersect(c('tibble','data.frame','matrix')) %>% length > 0){
                if (!is.null(title_line) && title_line!=""){
                    writeData(wb, sheet =  sheetname,x =  paste0('## ',title_line), startRow = 1, startCol = 1)
                    writeData(wb, sheet =  sheetname,x =  fc[[i]], startRow = 2, startCol = 1,rowNames=rowNames)
                }else{
                    writeData(wb, sheet =  sheetname,x =  fc[[i]], startRow = 2, startCol = 1,rowNames=rowNames)
                }
            }else{
                ylog('Skip [',tt,']',sheetname,', cause it is in an unsupport format')
            }
        }
        saveWorkbook(wb, file = fpath, overwrite = TRUE)

    } else if (l == 1 && t=="call") {
        ext = 'pdf'
        fpath = file.path(outputdir, paste0(name ,'.', ext))

        tryCatch(expr = {
            cairo_pdf(fpath, width = WIDTH, height = HEIGHT)
            eval(fc, envir = .GlobalEnv)
        }, finally = {
            dev.off()
        })
    } else if (c('ggplot','pheatmap','Heatmap') %>% intersect(t) %>% length > 0) {
        ext = 'pdf'
        fpath = file.path(outputdir, paste0(name ,'.', ext))

        tryCatch(expr = {
            cairo_pdf(fpath, width = WIDTH, height = HEIGHT, family=family)
            print(fc)
        }, finally = {
            dev.off()
        })
    } else if (c('tibble','data.frame','matrix') %>% intersect(t) %>% length > 0){
        ext = 'xlsx'
        fpath = file.path(outputdir, paste0(name ,'.', ext))
        title_line = table_title_lines[[1]]
        rowNames = yhas_rownames(fc)
        wb = createWorkbook()
        addWorksheet(wb,'Sheet1')
        if (!is.null(title_line) && title_line!=""){
            writeData(wb, sheet =  'Sheet1',x =  paste0('## ',title_line), startRow = 1, startCol = 1)
            writeData(wb, sheet =  'Sheet1',x =  fc, startRow = 2, startCol = 1,rowNames=rowNames)
        }else{
            writeData(wb, sheet =  'Sheet1',x =  fc, startRow = 2, startCol = 1,rowNames=rowNames)
        }
        saveWorkbook(wb, file = fpath, overwrite = TRUE)
    }
    ylog("write 1", ext, "at", fpath, verbose = TRUE)
}
#-----------------

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

#-----------------

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

#-----------------

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

# GO analysis using different expression results
#' Title
#'
#' @param diff a data.frame which is a output df of DESeq2 different expr analysis output, with columns: symbol, log2FC, pvalue, padj, stat
#' @param gene_list list, gene symbols as names and `stat` or `p-value` as values, decreasing arranged. will overrite the diff as input
#' @param set set in c('all','all3','gsea','hm','kegg','immune','immune.jnj')
#' 'all' - all the gmt sets in the  refs, currently: Hallmark, KEGG, Immune, Immune2(immune.jnj)
#' 'all3' - hm, kegg, immune do the 3 gmt sets
#' 'gsea' - hm, kegg do the 2 gmt sets
#' 'hm','kegg','immune','immune.jnj' use the corresponding single gmt set
#' @param p.t pvalue threshold, default 0.01
#' @param p.adjust.t padj threshold, default 0.05
#' @param log2fc.abs.t log2FC threshold, default 1
#' @param ...
#'
#' @return a list obj of class 'y.GO.res', used for input for yplot_GO_res
#' @export
#'
#' @examples
ydo_gsea = function (diff = NULL, gene_list = NULL, set = "gsea", p.adjust.t = 0.05,
                     p.t = NA, minGSSize = 10, maxGSSize = 500, .lineheight = 0.75, pvalueCutoff = 0.2,
                     ...)
{
    "set in c('all','all3','gsea','hm','kegg','immune','immune.jnj') "
    npg = pal_npg()(10)[2:1]
    .T_ylab_wrap = 40
    refs = list(hm = c(figureTitle = "DiffExpr Hallmark enricment",
                       Name = "Hallmark"), kegg = c(figureTitle = "DiffExpr KEGG enrichment",
                                                    Name = "KEGG"), immune = c(figureTitle = "DiffExpr Immune Infiltration enrichment",
                                                                               Name = "Immune"), immune.jnj = c(figureTitle = "DiffExpr Immune Infiltration2 enrichment",
                                                                                                                Name = "Immune2"))
    if (set == "all") {
        run = refs
    }
    else if (set == "gsea") {
        run = refs[c("hm", "kegg")]
    }
    else if (set == "all3") {
        run = refs[c("hm", "kegg", "immune")]
    }
    else if (set %in% c("hm", "kegg", "immune", "immune.jnj")) {
        run = refs[set]
    }
    else {
        stop("set must in in c('all','all3','gsea','hm','kegg','immune','immune.jnj')")
    }
    res = list()
    if (!is.null(diff)) {
        if (!("symbol" %in% colnames(diff))) {
            diff = diff %>% rownames_to_column("symbol")
        }
        gene_list = diff %>% distinct(symbol, .keep_all = TRUE) %>%
            pull(stat, name = symbol) %>% sort(decreasing = TRUE)
    }
    else if (!is.null(gene_list)) {
    }
    else {
        stop("input diff / gene_list can not be all NA")
    }
    for (key in names(run)) {
        v = run[[key]]
        figureTitle = v["figureTitle"]
        Name = v["Name"]
        gmt = yload_gmt(key)
        res.gsea = GSEA(geneList = gene_list, TERM2GENE = gmt,
                        minGSSize = minGSSize, maxGSSize = maxGSSize, pAdjustMethod = "BH", pvalueCutoff=pvalueCutoff,
                        ...)
        res.tb.gsea = res.gsea@result %>% dplyr::filter(p.adjust <
                                                            p.adjust.t)
        n.res = res.tb.gsea %>% nrow
        print(paste("after filter,", n.res, "results left"))
        if (n.res > 0) {
            res.tb.gsea %<>% mutate(Description = fct_reorder(Description,
                                                              p.adjust, .desc = TRUE), EnrichedGeneNum = core_enrichment %>%
                                        str_split("/") %>% sapply(length), direction = factor(sign(NES),
                                                                                              levels = c(-1, 1), labels = c("DEC Genes", "INS Genes")),
                                    EnrichedPer = round(EnrichedGeneNum/setSize *
                                                            100)) %>% dplyr::select(!core_enrichment)
            tab = res.tb.gsea %>% dplyr::select(Description,
                                                EnrichedGeneNum, EnrichedPer)
            maxlen_ylabel = res.tb.gsea$ID %>% sapply(str_length) %>%
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
            h = ((round(min(30, n.res)/4 + 3) * height_ratio) *
                     2)/2
            gg.gsea = res.tb.gsea %>% head(30) %>% ggplot(aes(x = -log10(pvalue),
                                                              y = Description, color = direction)) + geom_point(aes(size = EnrichedPer),
                                                                                                                alpha = 0.75) + scale_shape_manual(values = c(15,
                                                                                                                                                              16, 17, 18)) + scale_color_manual(breaks = c("DEC Genes",
                                                                                                                                                                                                           "INS Genes"), values = npg) + ylab(NULL) + ggtitle(figureTitle) +
                theme_bw(base_size = 16) + theme(axis.text.y.left = element_text(size = ylabel_size,
                                                                                 lineheight = .lineheight)) + scale_y_discrete(position = "left",
                                                                                                                               labels = function(x) str_wrap(str_replace_all(x,
                                                                                                                                                                             fixed("_"), " "), width = maxlen_ylabel)) +
                guides(y.sec = ggh4x::guide_axis_manual(breaks = tab$Description,
                                                        labels = paste0(tab$EnrichedGeneNum, "(", tab$EnrichedPer,
                                                                        "%)")), size = guide_legend(title = "Enriched%"),
                       color = guide_legend(title = "Direction"))
            make.custom(w, h)
            figure_size = c(w, h)
            res[[key]] = list(gsea.set = Name, gene_list = gene_list,
                              res.gsea = res.gsea, tb.gsea = res.tb.gsea, gg = gg.gsea,
                              figsize = figure_size, w = w, h = h, height_ratio = height_ratio,
                              maxlen_ylabel = maxlen_ylabel)
        }
        else {
            print(paste0("[!] Handling ", key, ": too less res got after filtering, consider to increase @p.adjust.t"))
            res[[key]] = list(gsea.set = Name, gene_list = gene_list,
                              res.gsea = res.gsea, res.tb.gsea = res.tb.gsea,
                              gg = NA, figsize = NA)
        }
    }
    message("GSEA analysis done, now you can use yplot_")
    res
}


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
#-----------------

#' Title
#'
#' @param x
#' @param fname
#' @param export
#' @param outputdir
#' @param ext
#' @param prefix
#' @param flag
#' @param worker
#' @param verbose
#' @param suffix
#' @param mkdir
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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

.ydumpto.call = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,row.names=FALSE,...){
    # if (row.names %>% is.null) row.names= FALSE
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
.ydumpto.data.frame = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,...){
    # if (row.names %>% is.null) row.names= TRUE
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

.ydumpto.matrix = function(x,fname,ext=NULL,outputdir=OUTPUTROOT,verbose=TRUE,...){
    # if (row.names %>% is.null) row.names= TRUE
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
#-----------------

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
#-----------------

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
#-----------------

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
#-----------------

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
#-----------------

library("survival")
library("survminer")

ydo_single_factor_coxph = function(data, covariates, key=NULL, join='_', tags=c('month','status')){
    "MAKE SURE that data values consist of int, which 0 mean no mut, others mean mut
 except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)
 mut_cases = column$values != 0
 wt_cases = total - mut_cases
"

    if (! is.null(key)){
        templ = stringr::str_flatten(c('Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'Surv(survivalMonth,survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(x) {
                                res = list()
                                formu=as.formula(paste(templ, x))
                                res[['formula']] = formu
                                res[['feature']] = x
                                res[['model']] = coxph(formu, data = data,)
                                res
                            },simplify = FALSE)

    univ_results <- lapply(univ_formulas,
                           function(li){
                               feature = li$feature
                               col = data %>% dplyr::pull(!!feature)
                               mut_cases = sum(col!=0,na.rm = TRUE)
                               wt_cases = sum(col==0,na.rm = TRUE)
                               all_cases = length(col)
                               if (all_cases > mut_cases + wt_cases){
                                   # col contains NA
                                   note = paste('contains',all_cases - mut_cases - wt_cases,'NA values')
                               }else{
                                   note = ""
                               }
                               ratio = signif(mut_cases/(mut_cases+wt_cases),digits = 3) * 100
                               x = li$model
                               x <- summary(x)
                               #获取p值
                               p.value<-signif(x$wald["pvalue"], digits=4)
                               p.sig = ymark_psig(p.value)
                               #获取HR
                               HR <-signif(x$coef[2], digits=2);
                               #获取95%置信区间
                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                               HR_95_CI_lower = HR.confint.lower
                               HR_95_CI_upper = HR.confint.upper
                               #feature = r$conf.int %>% row.names

                               res<-c(mut_cases, wt_cases, ratio, p.value
                                      , p.sig
                                      , HR, HR_95_CI_lower, HR_95_CI_upper,note)
                               names(res)<-c('mut_cases','wt_cases','mut_ratio%', "p_value"
                                             , 'p_sig'
                                             , "HR", 'HR_95_CI_lower', 'HR_95_CI_upper','note')
                               return(res)
                           })
    #转换成数据框，并转置
    res = t(as.data.frame(univ_results, check.names = FALSE))
    as.data.frame(res) %>% arrange(p_value)
}


ydo_single_factor_logrank = function(data, covariates, key=NULL, join='_', tags=c('month','status'), dev=FALSE){
    "MAKE SURE that data values consist of int, which 0 mean no mut, others mean mut
 except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)
 mut_cases = column$values != 0
 wt_cases = total - mut_cases
"
    if (! is.null(key)){
        templ = stringr::str_flatten(c('Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'Surv(survivalMonth, survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(feature) {
                                res = list()
                                formu=as.formula(paste(templ, feature))
                                res[['formula']] = formu
                                res[['feature']] = feature
                                res
                            },simplify = FALSE)
    if (dev ==TRUE){
        return(univ_formulas)
    }
    univ_results <- lapply(univ_formulas,
                           function(li){
                               feature = li$feature
                               col = data %>% dplyr::pull(!!feature)
                               mut_cases = sum(col!=0,na.rm = TRUE)
                               wt_cases = sum(col==0,na.rm = TRUE)
                               all_cases = length(col)
                               if (all_cases > mut_cases + wt_cases){
                                   # col contains NA
                                   note = paste('contains',all_cases - mut_cases - wt_cases,'NA values')
                               }else{
                                   note = ""
                               }
                               if (mut_cases==0 || wt_cases==0){
                                   res<-c(mut_cases, wt_cases, NA, NA
                                          , NA
                                          ,'A(1)/B(0)', NA,NA,NA)
                               }else{
                                   ratio = signif(mut_cases/(mut_cases+wt_cases),digits = 3) * 100
                                   model = survdiff(li[['formula']], data = data)
                                   # model is like
                                   #     survdiff(formula = formu, data = data)
                                   #
                                   #                 N Observed Expected (O-E)^2/E (O-E)^2/V
                                   # Radiotherapy=0 11        4    0.798      12.8      16.1
                                   # Radiotherapy=1 36        0    3.202       3.2      16.1
                                   #
                                   # Chisq= 16.1  on 1 degrees of freedom, p= 6e-05

                                   #获取p值
                                   p.value = 1 - pchisq(model$chisq, length(model$n) - 1)
                                   p.value = signif(p.value,4)
                                   p.sig = ymark_psig(p.value)
                                   #获取HR = A/B = (O(A)/E(B))/(O(B)/E(B))
                                   # here A is feature=1, B is feature=0
                                   HR = (model$obs[2]/model$exp[2])/(model$obs[1]/model$exp[1])
                                   #获取95%置信区间
                                   low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/model$exp[2]+1/model$exp[1]))
                                   up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/model$exp[2]+1/model$exp[1]))
                                   HR_95_CI_lower = signif(low95,2)
                                   HR_95_CI_upper = signif(up95,2)
                                   #feature = r$conf.int %>% row.names

                                   res<-c(mut_cases, wt_cases, ratio, p.value
                                          , p.sig
                                          ,'A(1)/B(0)', HR, HR_95_CI_lower, HR_95_CI_upper,note)
                               }

                               names(res)<-c('mut_cases','wt_cases','mut_ratio%', "p_value"
                                             , 'p_sig'
                                             ,'HR_formular', "HR", 'HR_95_CI_lower', 'HR_95_CI_upper',"note")
                               return(res)
                           })
    #转换成数据框，并转置
    res = t(as.data.frame(univ_results, check.names = FALSE))
    as.data.frame(res) %>% arrange(p_value)
}
#-----------------

ydo_single_factor_coxph2_helper = function(GRPS = NULL){
    " if GRPS == NULL: generate helper message
if GRPS is character: assign a global var ITER and returns a vector COLNAME
"
    if (is.null(GRPS)){
        print(" -- demo -- ")
        print("global var ITER is needed to run ysingle_factor_coxph2")
        print("GRPS = c('Good','Poor','Others');    ITER = combn(GRPS %>% sort,2,simplify = FALSE)")
        print('COLNAME is like:')
        print("COLNAME = c(paste0('cases_',GRPS,'_d'),'N_d',paste0('ratio%_',GRPS,'_d'),'p_value_d', 'p_sig', paste0('pval_',NAMES,'_d'),'Compares','HR.mean_d', paste0('HR_',NAMES,'_d'),paste0('CI_',NAMES),'note')")
    }else{
        ITER = combn(GRPS %>% sort,2,simplify = FALSE)
        NAMES = c()
        for (i in ITER){
            NAMES = c(NAMES,i %>% str_flatten(collapse = '_vs_'))
        }
        COLNAME = c(
            paste0('cases_',GRPS,'_d'),
            'N_d',
            paste0('ratio%_',GRPS,'_d'),
            "p_value_d", "p_sig",
            paste0('pval_',NAMES,'_d'),
            'Compares',
            "HR.mean_d",
            paste0('HR_',NAMES,'_d'),
            paste0('CI_',NAMES),
            "note"
        )
        assign('ITER',ITER)
        COLNAME
    }
}
ydo_single_factor_coxph2 = function (data, covariates, key = NULL, test.method='coxph', join = "_", tags = c("month", "status")) {
    "MAKE SURE that data values consist of int,
    which 0 mean no mut, others mean mut\n
    except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)\n
    mut_cases = column$values != 0\n wt_cases = total - mut_cases\n
    test.method in c('coxph','logrank','likelihood')
        coxph: wald-test
        logrank: logrank-test
"
    `+` = .Primitive("+")
    if (is.null(ITER)){
        ITER = data[,covariates]
    }
    if (!is.null(key)) {
        templ = stringr::str_flatten(c("Surv(", key, join, tags[[1]],
                                       ",", key, join, tags[[2]], ")~"))
    }
    else {
        templ = "Surv(survivalMonth,survivalEvent)~"
    }

    univ_formulas <- sapply(covariates, function(x) {
        res = list()
        formu = as.formula(paste(templ, x))
        res[["formula"]] = formu
        res[["feature"]] = x
        res
    }, simplify = FALSE)
    univ_results <- lapply(univ_formulas, function(li) {
        feature = li$feature
        col = data %>% dplyr::pull(!!feature)
        formu = li$formu

        grps = col %>% unique %>% .[!is.na(.)]

        ps = c()
        cmps = c()
        CIs = c()
        HRs = c()

        N = length(grps)

        nums = c(
            Good=sum(col == 'Good', na.rm = TRUE),
            Poor=sum(col == 'Poor', na.rm = TRUE),
            Others=sum(col=='Others', na.rm = TRUE)
        )
        all_cases = length(col)

        if (all_cases > sum(nums)) {
            note = paste("contains", all_cases - sum(nums), "NA values")
        } else {
            note = ""
        }
        ratios = signif(nums/all_cases,digits = 3) * 100

        for (x in ITER){
            a = x[[1]]
            b = x[[2]]
            if (nums[[a]]==0 | nums[[b]]==0){
                p.value = NA
                HR = NA
                HR_95_CI_lower = NA
                HR_95_CI_upper = NA
            }else{
                d = data[col %in% c(a,b),]
                model = coxph(formu, data = d )
                x <- summary(model)
                if (test.method == "logrank") {
                    p.value = x$sctest[["pvalue"]]
                }
                else if (test.method == "coxph" | test.method=='wald') {
                    p.value = x$waldtest[["pvalue"]]
                }
                else if (test.method == "likelihood") {
                    p.value = x$logtest[["pvalue"]]
                }
                p.value <- signif(p.value, digits = 4)
                HR <- signif(x$coef[2], digits = 2)
                HR_95_CI_lower <- signif(x$conf.int[, "lower .95"], 2)
                HR_95_CI_upper <- signif(x$conf.int[, "upper .95"], 2)
            }
            ps = ps %>% append(p.value)
            cmps = cmps %>% append(paste0(a,'_vs_',b))
            CIs = CIs %>% append(paste0('[',HR_95_CI_lower,' - ',HR_95_CI_upper,']'))
            HRs = HRs %>% append(HR)
        }

        p.value = mean(ps,na.rm = TRUE)
        # all_pval = ps %>% str_flatten(collapse = ', ')
        p.sig = ymark_psig(p.value)
        cmp.str = cmps %>% str_flatten(collapse = ', ')
        # CI = CIs %>% str_flatten(collapse = ', ')
        HR.mean = mean(HRs,na.rm = TRUE)
        # HR = HRs %>% str_flatten(collapse = ', ')


        res <- c(nums, N,ratios,
                 p.value, p.sig , ps, cmp.str,
                 HR.mean, HRs, CIs, note)
        names(res) = NULL
        #         names(res) = COLNAME
        return(res)
    })
    #     as.data.frame(univ_results, check.names = FALSE)
    return(univ_results)
    res = t(as.data.frame(univ_results, check.names = FALSE))  %>% as.data.frame
    colnames(res) = COLNAME
    res
}



ydo_single_factor_coxph3 = function (data, covariates, key = NULL, test.method = "coxph",
                                     join = "_", tags = c("month", "status"))
{
    "MAKE SURE that data values consist of int, \n    which 0 mean no mut, others mean mut\n \n    except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)\n \n    mut_cases = column$values != 0\n wt_cases = total - mut_cases\n\n    test.method in c('coxph','logrank','likelihood')\n        coxph: wald-test\n        logrank: logrank-test\n"
    `+` = .Primitive("+")


    if (!is.null(key)) {
        templ = stringr::str_flatten(c("Surv(", key, join, tags[[1]],
                                       ",", key, join, tags[[2]], ")~"))
    }
    else {
        templ = "Surv(survivalMonth,survivalEvent)~"
    }
    univ_formulas <- sapply(covariates, function(x) {
        res = list()
        formu = as.formula(paste(templ, x))
        res[["formula"]] = formu
        res[["feature"]] = x
        res
    }, simplify = FALSE)
    univ_results <- lapply(univ_formulas, function(li) {
        feature = li$feature
        col = data %>% dplyr::pull(!!feature)
        formu = li$formu
        grps = col %>% unique %>% .[!is.na(.)]
        ps = c()
        cmps = c()
        CIs = c()
        HRs = c()
        N = length(grps)
        nums = c()
        for (g in grps){
            nums = c(nums,sum(col == g, na.rm = TRUE))
        }
        # 使用grps_names,为了绕开下述 R语言设计缺陷:
        # 用0作为names的时候,不能读取
        # x = list('a')
        # names(x) = 0 # OK
        # x[[0]] # raise Error: attempt to select less than one element in get1index <real>
        grps_names = grps %>% make.names
        names(nums) = grps_names

        all_cases = length(col)
        if (all_cases > sum(nums)) {
            note = paste("contains", all_cases - sum(nums), "NA values")
        }
        else {
            note = ""
        }
        # R语言设计缺陷:
        # 用0作为names的时候,不能读取
        # x = list('a')
        # names(x) = 0 # OK
        # x[[0]] # raise Error: attempt to select less than one element in get1index <real>
        ratios = signif(nums/all_cases, digits = 3) * 100
        ITER = combn(grps_names,2,simplify = FALSE)
        # print(paste(feature,'-',nums))
        for (x in ITER) {
            x = x %>% sort
            a = x[[1]]
            b = x[[2]]
            if (nums[[a]] == 0 | nums[[b]] == 0) {
                p.value = NA
                HR = NA
                HR_95_CI_lower = NA
                HR_95_CI_upper = NA
            }
            else {
                a_ori = grps[match(a, grps_names)]
                b_ori = grps[match(b, grps_names)]
                d = data[col %in% c(a_ori,b_ori), ]
                model = coxph(formu, data = d)
                x <- summary(model)
                if (test.method == "logrank") {
                    p.value = x$sctest[["pvalue"]]
                }
                else if (test.method == "coxph" | test.method ==
                         "wald") {
                    p.value = x$waldtest[["pvalue"]]
                }
                else if (test.method == "likelihood") {
                    p.value = x$logtest[["pvalue"]]
                }
                p.value <- signif(p.value, digits = 4)
                HR <- signif(x$coef[2], digits = 2)
                HR_95_CI_lower <- signif(x$conf.int[, "lower .95"],
                                         2)
                HR_95_CI_upper <- signif(x$conf.int[, "upper .95"],
                                         2)
            }
            ps = ps %>% append(p.value)
            cmps = cmps %>% append(paste0(b_ori, "_vs_", a_ori))
            CIs = CIs %>% append(paste0("[", HR_95_CI_lower,
                                        " - ", HR_95_CI_upper, "]"))
            HRs = HRs %>% append(HR)
        }
        p.value = mean(ps, na.rm = TRUE)
        p.sig = ymark_psig(p.value)
        cmp.str = cmps %>% str_flatten(collapse = ", ")
        HR.mean = mean(HRs, na.rm = TRUE)
        res <- c(cmp.str, nums, N, ratios,test.method, p.value, p.sig, ps,
                 HR.mean, HRs, CIs, note)
        names(res) = NULL
        return(res)
    })
    # return(univ_results)
    res = t(as.data.frame(univ_results, check.names = FALSE)) %>%
        as.data.frame
    colnames(res) = c('Cmp_Groups','LeftGroup_cases','Right_Group_cases','SubGroupNumbers'
                      ,'LeftGroup_ratio%', 'Right_Group_ratio%','stat_method',"p_value"
                      , 'p_sig','p_val_list'
                      , "HR.mean", 'HR_list', 'HR_95_CI_list',"note")
    res
}
#-----------------

library("ggfortify")
yplot_survival = function(data, covariates
                          ,key=NULL
                          ,test.method='wald'
                          ,tags=c('month','status')
                          ,auto_anno=TRUE
                          ,join='_'
){
    " test.method in c('logrank','wald','likelihood') "
    if (! is.null(key)){
        templ = stringr::str_flatten(c('Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'Surv(survivalMonth,survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(x) {
                                res = list()
                                formu=as.formula(paste(templ, x))
                                res[['formula']] = formu
                                res[['feature']] = x
                                res[['model']] = survfit(formu, data = data,)
                                diff = coxph(formu, data = data) %>% summary
                                if (test.method=='logrank'){
                                    p.value = diff$sctest[['pvalue']]
                                }else if (test.method=='wald'){
                                    p.value = diff$waldtest[['pvalue']]
                                }else if (test.method=='likelihood'){
                                    p.value = diff$logtest[['pvalue']]
                                }
                                res[['pval']] = signif(p.value,4)
                                res[['test.method']] = test.method
                                res
                            },simplify = FALSE)
    plots = lapply(univ_formulas, function(li){
        model = li$model
        # fn = feature name
        fn = li$feature
        pval = li$pval
        if (pval < 0.05){
            pval.color = 'red'
        }else{
            pval.color = 'black'
        }
        # cal annotation y depend on group number of feature fn
        grps = data[[fn]] %>% unique
        n = length(grps)

        gg = autoplot(model) +
            labs(color=fn,fill=fn) +
            scale_y_continuous(limits = c(0,1),labels = scales::percent)
        gg = gg +  annotate(
            geom = "text", x=0,y = 0.042*n,
            color=pval.color,
            label = paste('p =',pval,'\n'),
            hjust = "inward"
        )
        if (auto_anno==TRUE){
            df = univ_formulas[[fn]]$model %>% summary %>% .$table
            text_label = paste('[',df[,'records'],']',df %>% row.names,' median survival: ',df[,'median'])
            gg = gg + annotate(
                geom = "text", x=0,y= 0.02*n,
                label = text_label %>% stringr::str_flatten('\n'),
                hjust = "inward"
                #                     ,vjust = 'inward'
            )
        }

        gg
    })
    names(plots) = names(univ_formulas)
    list(data=univ_formulas,plots=plots)
}
#-----------------

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
#-----------------


make.tall = function(.set_global=TRUE){
    options(repr.plot.width=6
            ,repr.plot.height=8
            #         ,repr.matrix.max.cols=30
            #         ,repr.matrix.max.rows=10
    )
    if (.set_global==TRUE){
        assign('WIDTH',value = 6, envir = .GlobalEnv)
        assign('HEIGHT',value = 8, envir = .GlobalEnv)
    }
}
make.wide = function(.set_global=TRUE){
    options(repr.plot.width=12
            ,repr.plot.height=8
            #         ,repr.matrix.max.cols=30
            #         ,repr.matrix.max.rows=10
    )
    if (.set_global==TRUE){
        assign('WIDTH',value = 12, envir = .GlobalEnv)
        assign('HEIGHT',value = 8, envir = .GlobalEnv)
    }
}
make.custom = function(w,h,.set_global=TRUE){
    options(repr.plot.width=w
            ,repr.plot.height=h
            #         ,repr.matrix.max.cols=30
            #         ,repr.matrix.max.rows=10
    )
    if (.set_global==TRUE){
        assign('WIDTH',value = w, envir = .GlobalEnv)
        assign('HEIGHT',value = h, envir = .GlobalEnv)
    }
}
library("ggsci")
#-----------------

library("RColorBrewer")
yset_color_scheme = function(name='wjw'){
    gg_color_hue <- function(n,start=15) {
        hues = seq(start, 360 + start, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }
    if (name=='sh'){
        # 胜寒1
        COL_SNV = c("#464E2B", "#A8352A", "#B2DF8A",
                    "#33A02C", "#FB9A99", "#559FCD",
                    "#FDBF6F", "#F1AE3C", "#CAB2D6",
                    "#8dd3c7", "#6B53BA")
        names(COL_SNV) = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins",
                           "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins",
                           "5'Flank","Multi_Hit")

        COL_CNV = c("#464E2B", "#559FCD", "#B2DF8A",
                    "#33A02C", "#FB9A99", "#A8352A",
                    "#FDBF6F", "#559FCD", "#A8352A",
                    "#8dd3c7", "#6B4B98")
        names(COL_CNV) = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins",
                           "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins",
                           "5'Flank","Multi_Hit")

        # # 蓝灰 藏蓝 深蓝 浅蓝
        # COL_ANNO = c('#b1b8d4','#3c2e43','#073763','#0B5394')[1:length(grps)]
        # names(COL_ANNO) = grps

        COL_TITV = c("#36312F",'#8A7D7D','#A8352A','#559FCD','#DAE2ED','#F2AC3C')
        names(COL_TITV) = c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")

    }else if (name == 'wjw'){
        # 金旺1
        COL_SNV = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                    "#33A02C", "#FB9A99", "#E31A1C",
                    "#FDBF6F", "#FF7F00", "#CAB2D6",
                    "#8dd3c7", "#6A3D9A")
        names(COL_SNV) = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins",
                           "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins",
                           "5'Flank","Multi_Hit")

        COL_CNV = c("#3653A5","#BFC2D7","#E21A21","#E5B3B8",
                    "#464E2B", "#559FCD", "#B2DF8A", "#33A02C",
                    "#A8352A", "#8dd3c7", "#6B4B98")

        names(COL_CNV) = c("Frame_Shift_Del",'In_Frame_Del',"Frame_Shift_Ins","In_Frame_Ins",
                           "Missense_Mutation", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation","5'Flank","Multi_Hit")
        # COL_ANNO = c('#C2171D','#3A7FB7','#964E9B','#F67C1E','#E7D712',
        #             "#36312F",'#8A7D7D','#A8352A','#559FCD','#DAE2ED','#F2AC3C')[1:length(grps)]
        # names(COL_ANNO) = grps

        COL_TITV = gg_color_hue(6,start=0) #c("#36312F",'#8A7D7D','#A8352A','#559FCD','#DAE2ED','#F2AC3C')
        names(COL_TITV) = c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")
    } else{
        COL_SNV = c("#464E2B", "#A8352A", "#B2DF8A",
                    "#33A02C", "#FB9A99", "#559FCD",
                    "#FDBF6F", "#F1AE3C", "#CAB2D6",
                    "#8dd3c7", "#6B53BA")
        names(COL_SNV) = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins",
                           "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins",
                           "5'Flank","Multi_Hit")

        COL_CNV = c("#464E2B", "#559FCD", "#B2DF8A",
                    "#33A02C", "#FB9A99", "#A8352A",
                    "#FDBF6F", "#559FCD", "#A8352A",
                    "#8dd3c7", "#6B4B98")
        names(COL_CNV) = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins",
                           "Splice_Site", "Translation_Start_Site","Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins",
                           "5'Flank","Multi_Hit")
        # COL_ANNO = gg_color_hue(length(grps))
        # names(COL_ANNO) = grps
        COL_TITV = gg_color_hue(6,start=50)#c("#36312F",'#8A7D7D','#A8352A','#559FCD','#DAE2ED','#F2AC3C')
        names(COL_TITV) = c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G")
    }
    assign('COL_CNV', COL_CNV, envir = .GlobalEnv)
    assign('COL_SNV', COL_SNV, envir = .GlobalEnv)
    # assign('COL_ANNO', COL_ANNO, envir = .GlobalEnv)
    assign('COL_TITV', COL_TITV, envir = .GlobalEnv)
    print("global var COL_SNV, COL_CNV, COL_TITV have been set.")
}

yset_random_anno_colors = function(df, palette=NULL, cols = c('Clin_classification'), .assign_global=FALSE){
    if (is.null(palette)){
        palette = pal_npg()(10)
    }
    col.names = df %>% colnames
    COL_ANNO = list()
    for (col in cols){
        if (col %in% col.names){
            grps = df[[col]] %>% unique %>% sort(na.last = TRUE)
            if (col == 'Clin_classification' || col =='Group'){
                COL_ANNO[[col]] = palette[1:length(grps)]
            }else{
                COL_ANNO[[col]] = sample(palette,length(grps))
            }
            names(COL_ANNO[[col]]) = grps
            COL_ANNO[[col]][grps %>% is.na] = '#BBBBBB'
        }else{
            print(paste('‼️',col,'not in df %>% colnames'))
        }
    }
    if (.assign_global){
        assign('COL_ANNO', COL_ANNO, envir = .GlobalEnv)
        print("global var COL_ANNO has been set.")
    }else{
        return (COL_ANNO)
    }
}

#' @title yset_anno_color_names
#' @description random colors for annotation columns (set @param colors=NULL)or set with a list of index of colors in palette
#'  used for generate COL_ANNO for maftools::oncoplot(annotationColor = COL_ANNO)
#'
#' @param df data.frame
#' @param palette color palette, provide HEX color string vectors or use pal_*()(MAX_COLOR_NUM) to get a palette
#'  if palette is NULL, use palette = pal_png()(10)
#' @param cols generate colors for these columns in @param df
#' @param colors a list of index of colors in palette @param colors is a list with the same length of cols,
#' each element will act on the corresponding element of cols, if colors is NULL, random colors will be used, the funcion will invoke yset_random_anno_colors
#' for each element of colors, if  NA or NULL, use palette(1:length()) to generate colors,
#'  if 'random', use sample(palette, length()) to generate colors,
#'  if a integer vector, use palette[color] to generate colors,
#' @param .autoNAcolor if TRUE, set color of NA value to '#BBBBBB'
#' @param .assign_global assign to global env or not, if colors is NULL, .assign_global will be set to FALSE
#'
#' @return list of colors
#' @export
#' @examples
#' df = data.frame(a = c('a','a','b'), b = c('a','c','c'), c = c('a','b','c'))
#' yset_anno_color_names(df, cols = c('a','b','c'),)
ygen_anno_color = function(df, palette=NULL, colors=NULL, cols = c('Clin_classification'), .autoNAcolor = TRUE, .assign_global=TRUE){
    if (is.null(palette)){
        palette = pal_npg()(10)
    }
    if(is.null(colors)){
        return(yset_random_anno_colors(df=df, palette=palette, cols=cols, .assign_global=.assign_global))
    }
    col.names = df %>% colnames
    COL_ANNO = list()
    L = length(cols)
    for (i in 1:L){
        col = cols[[i]]
        color = colors[[i]]
        if (col %in% col.names){
            grps = df[[col]] %>% unique %>% sort(na.last = TRUE)
            if (is.null(color) || is.na(color)){
                COL_ANNO[[col]] = palette[1:length(grps)]
            }else if (length(color)==1 && color=='random'){
                COL_ANNO[[col]] = palette %>% sample(length(grps))
            }else{
                # stopifnot(is.integer(color))
                COL_ANNO[[col]] = palette[color]
            }
            names(COL_ANNO[[col]]) = grps
            if (.autoNAcolor){
                COL_ANNO[[col]][grps %>% is.na] = '#BBBBBB'
            }
        }else{
            print(paste('‼️',col,'not in df %>% colnames'))
        }
    }
    if (.assign_global){
        assign('COL_ANNO', COL_ANNO, envir = .GlobalEnv)
        print("global var COL_ANNO has been set.")
    }else{
        return (COL_ANNO)
    }
}

#' @title show_anno_col
#' @description show colors for COL_ANNO
#' @param ... allowed value are c(labels = FALSE,borders = FALSE,cex_label = ?)
show_anno_col = function(col_anno,...){
    res = c()
    l = col_anno %>% sapply(length) %>% max
    for (i in col_anno){
        res = c(res,i,rep('#FFFFFF',l-length(i)))
    }
    print(col_anno)
    show_col(res,ncol = l,...)
}
#-----------------

# `+` = function(a,b){
#     if (a %>% is.character && b %>% is.character){
#         return (paste0(a,b))
#     }else{
#         return(.Primitive("+")(a,b))
#     }
# }

#-----------------

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
ypsig_mark = function(p){
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

ymark_psig = ypsig_mark

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

#-----------------

