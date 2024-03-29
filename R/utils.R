#' same as Reduce
#'
#' @description same as Reduce except for x,f position exchanged for compilable with pipe op %>%
#'
#' @param x
#' @param f
#' @param init
#' @param ...
#'
#' @return reduce results
#'
#' @examples
reduce = function(x,f,init,...){
    Reduce(f,x,init,...)
}


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
#' @return Nothing is returned
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
#' @return Nothing is returned
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
#' @return last object in list .obj
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
yhas_colnames = function(df, col=NULL){
    stopifnot(is.data.frame(df))

    cn = Negate(is.null)(names(df))
    if (is.null(col)){
        reuturn(cn)
    }else{

    }
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
    View(mat)
    options(repr.matrix.max.cols=cur_ncol)
    options(repr.matrix.max.rows=cur_nrow)
}

#' fill NA values in a mat, prefer to be used with pipe \%>\% operators
#'
#' @description  can be easily placed in a pipe line \%>\%
#'
#' @param mat a matrix or data.frame
#' @param value value to replace NA values in the mat
#'
#' @return mat with NA filled with value
#' @export
#'
#' @examples
fill_na = function(mat,value){
    mat[is.na(mat)] <- value
    mat
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
    if ('MAF' %in% class(maf)){
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


#' Title
#'
#' @param maf_ a MAF style table file name
#' @param gl_ a group list file name
#' @param laml a maftools::MAF object
#' @param only_12_cols only use first 12 columns of read maf_ table (maftools standard input columns numbers)
#' @param ...
#'
#' @return maftools::MAF object
#' @export
#'
#' @examples
yload_laml_maf = function(maf_=NULL,gl_=NULL,laml=NULL,only_12_cols=TRUE,...){
    if (is.null(laml) && !is.null(maf_) && !is.null(gl_)) {
        maf_ = yload_dfx(maf_,...)
        if (only_12_cols==TRUE) maf_ = maf_[,1:12]
        gl_ = yload_dfx(gl_,...)
        laml = maftools::read.maf(maf_,clinicalData = gl_)
    }else if(!is.null(laml)){
        #pass
    }else stop('laml,maf_,gl_ all NULL')
    laml
}

#' Load all GC used gene symbols
#'
#' @return a character vector
#' @export
#'
#' @examples
yload_symbols_all_genes_18669 = function(){
    f = system.file('extdata/all_genes_symbol.RDS',package = 'GCrutils')
    readRDS(f)
}
#' Load GC curated cancer-related gene symbols
#'
#' @return a character vector
#' @export
#'
#' @examples
#' cr = yload_symbols_cancer_related_genes_1190()
yload_symbols_cancer_related_genes_1190 = function(){
    f = system.file('extdata/cancer_related_symbols.RDS',package = 'GCrutils')
    readRDS(f)
}


## Function Description:
##

#' Check whether a character could be safely convert to numeric
#'
#' @description   A function to assess if a vector can be interpreted as numbers. Code copied from
#'   varhandle::check.numeric.
#'
#' @details  This function checks if it is safe to convert the vector to numeric and this conversion
#'   will not end up in producing NA. In nutshell this function tries to mak sure provided vector
#'   contains numbers but in a non-numeric class. see example for better understanding. This
#'   function can be configured to only accept integer numbers (by setting the argument only.integer
#'   to TRUE).
#'
#'   It can also ignore NA values (na.rm argument) and ignore heading/tailing whitespaces
#'   (ignore.whitespace argument).
#'
#'   There is also room to manually define exceptions to be concidered
#'   as numbers (exceptions argument).
#'
#' @param v The character vector or factor vector. (Mandatory)
#' @param na.rm logical. Should the function ignore NA? Default value is FLASE since NA can be
#'   converted to numeric. (Optional)
#' @param only.integer logical. Only check for integers and do not accept floating point. Default
#'   value is FALSE. (Optional)
#' @param exceptions A character vector containing the strings that should be considered as valid to
#'   be converted to numeric. (Optional)
#' @param ignore.whitespace logical. Ignore leading and tailing whitespace characters before
#'   assessing if the vector can be converted to numeric. Default value is TRUE. (Optional)
#'
#' @return vector which could be safely as.numeric \?
#' @export
#'
#' @examples
check.numeric = function(v = NULL, na.rm = FALSE, only.integer = FALSE,
                          exceptions=c(""), ignore.whitespace = TRUE){
    #----[ checking the input ]----#
    {
        # if the only.integer is NOT a single TRUE or FALSE
        if (!is.logical(only.integer) | length(only.integer) != 1) {
            # complain
            stop("The parameter \"only.integer\" should be either TRUE or FALSE.")
        }

        # if user has not defined the vector v
        if (is.null(v)) {
            # complain
            stop("The parameter \"v\" is not defined. It can be character vector, numeric vector, factor vector or logical vector.")
            # if user has defined but the class is NOT character or factor
        }else if (!inherits(v, c("character", "factor"))) {
            # if the class is NOT numeric or integer either
            if (!inherits(v, c("numeric", "integer", "logical"))) {
                # complain
                stop("The parameter \"v\" can only be a character vector, numeric vector, factor vector or logical vector.")
                # if the class is numeric or integer
            }else{
                # if user wants to specifically filter out non-integers, there
                # is a chance that the vector contains some non-integer numbers
                # so we should turn the vector to character and run the function
                if(only.integer){
                    # convert the vector to character
                    v <- as.character(v)
                }else{
                    # since it is already a number
                    return(rep(x = TRUE, length(v)))
                }
            }
        }

        # if the na.rm is NOT a single TRUE or FALSE
        if (!is.logical(na.rm) | length(na.rm) != 1) {
            # complain
            stop("The parameter \"na.rm\" should be either TRUE or FALSE.")
        }



        # if the ignore.whitespace is NOT a single TRUE or FALSE
        if (!is.logical(ignore.whitespace) | length(ignore.whitespace) != 1) {
            # complain
            stop("The parameter \"ignore.whitespace\" should be either TRUE or FALSE.")
        }
    }


    #----[ pre-processing ]----#
    {
        # convert to character if it is vector
        if (inherits(v, "factor")) {
            # convert to character
            v <- as.character(v)
        }

        # if user wants to ignore NAs
        if (na.rm) {
            # if it has some NAs
            if (any(is.na(v))) {
                # remove NAs
                v <- v[-pin.na(v)]
            }
        }

        # if user wants to ignore leading or tailing white space
        if (ignore.whitespace) {
            # substitute whitespaces in the begining and at the ending of each item in v
            v <- gsub("^\\s+|\\s+$", "", v)
        }
    }


    #----[ processing ]----#
    {
        # if user wants to only detect integers
        if (only.integer) {
            regexp_pattern <- "(^(-|\\+)?\\d+$)|(^(-|\\+)?(\\d*)e(-|\\+)?(\\d+)$)"
            # if user wants to detect all numbers
        }else{
            #regexp_pattern <- "^(-|\\+)?\\d+(\\.?\\d+)?$"
            regexp_pattern <- "(^(-|\\+)?((\\.?\\d+)|(\\d+\\.\\d+)|(\\d+\\.?))$)|(^(-|\\+)?((\\.?\\d+)|(\\d+\\.\\d+)|(\\d+\\.?))e(-|\\+)?(\\d+)$)"
        }

        # perform the regexp
        output <- grepl(pattern = regexp_pattern, x = v)

        # check for existance of exceptions
        exception_index <- is.element(v, exceptions)
        # if there are is exception detected
        if (any(exception_index)) {
            # turn their output value to TRUE
            output[exception_index] <- TRUE
        }

        # if user wants to keep NA
        if (!na.rm) {
            # NAs are marked as FALSE by grepl and we replace it with TRUE instead
            output[is.na(v)] <- TRUE
        }


        # return the result
        return(output)
    }
}



# ----------- view utils ----------------

#' Set the size for the following plots
#' @description  basically, this function affects the following plot function which will use global var WIDTH and HEIGHT
#'   to determine the display repr.plot.width and repr.plot.height, set global variables depend on the third parameter
#'   `.set_global`
#' @param w set global WIDHT,  in inch, if length(w) == 2
#' @param h set global HEIGHT, in inch
#' @param .set_global TRUE, set global vars, else no global variables set
#'
#' @return Nothing is returned
#' @export
#'
#' @examples
make.custom = function (w=NULL, h=NULL, .set_global = TRUE) {
    if (is.null(w) && is.null(h)){
        return(c(WIDTH,HEIGHT))
    }
    l = length(w)
    if (l==2){
        h = w[[2]]
        w = w[[1]]
    }else if(l == 1){
        stopifnot(!is.null(h))
    }else{
        stop('wrong w,h value')
    }
    options(repr.plot.width = w, repr.plot.height = h)
    if (.set_global == TRUE) {
        assign("WIDTH", value = w, envir = .GlobalEnv)
        assign("HEIGHT", value = h, envir = .GlobalEnv)
    }
}


#' Title
#'
#' @param .set_global
#'
#' @return Nothing is returned
#' @export
#'
#' @examples
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
#' Title
#'
#' @param .set_global
#'
#' @return
#' @export
#'
#' @examples
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

#' Generate color palette for annotation
#'
#' @param annotation
#' @param annotation_colors
#' @param palette
#' @param show_all_ggsci_palette
#'
#' @return list of colors
#' @export
#'
#' @examples
ygenerate_annotation_colors = function(annotation,
                                       annotation_colors=NA,
                                       palette="npg",
                                       show_all_ggsci_palette=FALSE){
    ggsci_palette = c(npg=10,aaas=10,nejm=8,
                      lancet=9,jama=7,jco=10,
                      ucscgb=26,d=310,locuszoom=7,
                      igv=50,uchicago=9,startrek=7,
                      tron=7,futurama=12,rickandmorty=12,
                      simpsons=16,gsea=12,material=10)
    if (show_all_ggsci_palette==TRUE){
        print('ALL available ggsci_palette:')
        return(ggsci_palette)
    }
    if (!is.na(annotation_colors) && !is.null(annotation_colors)){
        already_set = names(annotation_colors)
    }
    else already_set = c()
    if (palette %in% names(ggsci_palette)){
        fun = do.call(what = `::`
                      ,args = list("ggsci",str_flatten(c('pal_',palette))))
        pal = fun()(ggsci_palette[palette])
    }
    # stat total subgroups
    n = 0
    grps = list()
    for (colname in colnames(annotation)){
        grps[[colname]] = annotation[,colname] %>% unique
        n = n + length(grps[[colname]])
    }
    anno_colors = list()
    pal = rep(pal,n %/% (length(pal)+1) +1)
    s = 1
    for (g in names(grps)){
        if (g %in% already_set){
            anno_colors[[g]] = annotation_colors[[g]]
        }else{
            grp = grps[[g]]
            e = length(grp)
            cl = pal[s:(s+e-1)]
            names(cl) = grp
            anno_colors[[g]] = cl
            s = s + e
        }
    }
    anno_colors
}



#' @export
print.y.GSEA.plots <- function(x) {
    # Display only the ggplot object
    for (n in names(x)){
        IRdisplay::display(x[[n]]$gg)
    }
}

#' @export
print.y.GSEA.plot <- function(x) {
    # Display only the ggplot object
    make.custom(x$figsize)
    IRdisplay::display(x$gg)
}

#' @export
print.y.GSEA.resx <- function(x){
    print(paste(length(x),'y.GSEA.res objects.','You can use yplot_GSEA_res to generate and view the plots.'))
    for (n in names(x)){
        print(x[[n]])
    }
}

#' @export
print.y.GSEA.res <- function(x){
    msg = str_flatten(c('GSEA of ',x$gsea.set,' with ranked ',length(x$gene_list),' genes resulting in ',nrow(x$res.tb.gsea),'significant enrichments under filters:'))
    msg = c(msg, paste0('p.adjust.t = ',x$p.adjust.t),
            paste0('p.t = ',x$p.t),
            paste0('pvalueCutoff = ', x$pvalueCutoff)
    )
    print(str_flatten(msg,collapse = '  '))
}

#' Title
#'
#' @param pvalues a numeric vector
#'
#' @return character, * formatted pvalues
#' @export
#'
#' @examples
ypvalue_add_asterisks = function(pvalues) {
    pvalues %>% purrr::map_chr(function(p){
        if (p < 0.0001){
            return('****')
        }else if (p < 0.001) {
            return("***")
        } else if (p < 0.01) {
            return("**")
        } else if (p < 0.05) {
            return("*")
        } else {
            return("ns") # not significant
        }
    })
}



