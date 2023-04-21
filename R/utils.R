# -------------- IMPORT ----------------

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
        write(msg,.dev)
    }
}



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
