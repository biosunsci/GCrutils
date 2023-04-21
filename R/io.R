#' Load data.frame sheet from files
#'
#' @details 依赖于 yslice, ypush, ylog, ysplit_path, ysplit_file_name, yfile_path, <stringr> if `li` && `ext` both NULL
#'   and `frm` is not NULL, read all files from `frm`. Depends on: INPUTROOT, yfunc_args, read.table, %>%,
#'   modifyList, yslice, str_flatten, str_detect, fixed, ypush, str_sub, ypath_join
#'
#' @param li overwritting <pattern> settings a character vector, if length(<li>)>1 generate a named list else generate
#'   the dataframe read, return is determined further by <retmode> a character vector whose element <i> refers to dfx
#'   names, if <i> is not started with c('/','~','./'), find <i>.<ext> in folder <frm> else use <i> as the path, ignore
#'   <frm> for the current <i> if input li is a data.frame it will be returned immediately
#' @param frm a folder.path, all the data will be read from this path
#' @param pattern a RegExpr to glob, return matched files, only acts when `li` is NULL
#' @param ext basicly, <li> members do not include an extention name, ext is used to set the <li> extention name. set
#'   ext=NULL to only load exactly the name <li> provided. ext不会修改已经显式声明的扩展名，如li=c('abx.txt'),ext='dfx'
#'   则读入仍是'abx.txt'而不是'abx.dfx'或'
#' @param worker
#' @param verbose
#' @param row.names
#' @param ...
#' @param retmode
#'
#' @return df or list of dfs if retmode=='local', assign every df to .GlobalEnv if retmode=='global'
#' @export
#' @import dplyr
#' @import stringr
#' @examples
#'
yload_dfx = function(li = NULL,
                     frm = INPUTROOT,
                     pattern = NULL,
                     ext = NULL,
                     worker = NULL,
                     retmode = 'local',
                     verbose = TRUE,
                     row.names = NULL,
                     ...) {
    # retmode in c('local','global')
    # ext in c(str,NULL)
    # default parameters with custom one
    std_args = yfunc_args(utils::read.table)
    argv = list(
        sep = "\t",
        header = T,
        quote = "",
        stringsAsFactors = F,
        comment.char = "",
        na.strings = "",
        row.names = row.names
    ) %>% utils::modifyList(yfunc_args(), keep.null = TRUE) %>% ysubset_list_named(std_args)

    if ('data.frame' %in% (li %>% class))
        return(li)
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
            if (length(tmp) == 1) {
                ffname = tmp[[1]]
                myfrm = frm
            } else if (length(tmp) >= 2) {
                ffname = tmp %>% yslice(-1)
                myfrm = tmp %>% yslice(1:-2)
                if (myfrm[[1]] %in% c('/', '~', '.','..')) {
                    myfrm = myfrm %>% str_flatten(collapse = '/')
                } else{
                    myfrm = c(frm, myfrm) %>% str_flatten(collapse = '/')
                }
            } else{
                stop("item in @li is wrong")
            }
            # 检查i是否能作为变量名
            if (retmode == 'global' && make.names(i) != i) {
                # 不能
                ylog('[WARN] ignore ',
                     i,
                     ' for its not a legal var name',
                     .addtime = F)
            } else{
                tmp = ysplit_file_name(ffname, ext = ext, mode = 'auto')
                fname = tmp[[1]]
                if (length(tmp) == 1) {
                    ext = 'dfx'
                } else if (length(tmp) == 2) {
                    ext = tmp[[2]]
                } else{
                    stop('wrong tmp length')
                }
                if (fname %>% str_detect(fixed('.'))) {
                    write(
                        paste(
                            '[INFO] skip ',
                            i,
                            ' for its not right format: multi-dot in filename'
                        ),
                        2
                    )
                } else{
                    files = files %>% ypush(c(myfrm, fname, ext))
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
        if (vname %>% str_sub(1, 3) %>% str_detect('[A-Z][0-9]\\_'))
            worker = NULL
        if (!is.null(worker)) {
            ffname = str_flatten(c(worker, '_', vname, '.', ext))
            input = yfile_path(myfrm, ffname)
            if (!file.exists(input)) {
                ylog(
                    '[WARN] prefix with <worker> file not found, retry without <worker> one',
                    .addtime = FALSE
                )
                ffname = str_flatten(c(vname, '.', ext))
            }
        } else{
            ffname = str_flatten(c(vname, '.', ext))
        }
        input = yfile_path(myfrm, ffname)

        if (verbose == TRUE)
            ylog("read [", ffname , "] as '", vname, "' from > ", myfrm)

        if (ext == 'gmt') {
            # library(GSVA)
            original_gmt_GSVA <- readLines(input)
            strsplit_no_name <- function(gmt.list_layer) {
                    as.character(unlist(
                        strsplit(
                            gmt.list_layer,
                            split = '\t',
                            fixed = T
                        )
                    ))[-2]
                }
            ref_df <- lapply(original_gmt_GSVA, strsplit_no_name)
            for (layers in 1:length(ref_df)) {
                names(ref_df)[layers] <- ref_df[layers][[1]][1]
                ref_df[layers][[1]] <- ref_df[layers][[1]][-1]
            }
            res[[vname]] = ref_df
        } else{
            "ext == {txt,dfx,...} "
            argv$file = input
            res[[vname]] = utils::read.table %>% do.call(args = argv)
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
    } else if (retmode == 'global') {

    } else{
        # invalid retmode value
        stop('invalid retmode value')
    }
}

#' @export
db.importing = yload_dfx



#' Save plots to pdf, and data.frame list to xlsx
#'
#' @description auto save Robject into corresponding file format output file type is auto determined
#'   by the class(fc), currently support classes are 'ggplot', 'pheatmap','Heatmap' -> pdf files
#'   'call' -> eval(fc) -> pdf files 'data.frame', 'matrix' -> xlsx files 'list' of
#'   'data.frame'/'matrix's -> xlsx files with each element in the list -> a sheet in the xlsx
#'
#' @param fc the obj to export
#' @param table_title_lines a list of title lines, element i will be the title line in the
#'   corresponding table `fc[[i]]`, make sure the length of fc and table_title_lines are the same
#' @param outputdir output file into these dir
#' @param family used in grDevices::cairo_pdf(family=family)
#' @param ...
#' @param sep
#'
#' @return NULL
#' @export
#' @examples
#' \dontrun{
#'  data(mtcars)
#'  list(mtcars) %>% ysave('data','mtcars',outputdir='./') # file <'./data_mtcars.xlsx'> is created
#'  gg = ggplot(mtcars,aes(mpg,gear)) + geom_point()
#'
#'  # file <'./mtcars.pdf'> which is 10 inches wide and 12 inches tall is created
#'  make.custom(10,12)
#'  gg %>% ysave('mtcars',outputdir='./')
#' }
#'
#'
ysave = function (fc = NULL,
                  ...,
                  table_title_lines = list(),
                  outputdir = OUTPUTROOT,
                  sep = "_",
                  family = "sans")
{
    t = class(fc)
    l = length(t)
    name = str_flatten(c(...), collapse = sep)
    if (l == 1 && t == 'NULL') {
        return(c(
            outname = name,
            outputdir = outputdir,
            width = WIDTH,
            height = HEIGHT,
            sep = sep
        ))
    }
    if (name == "") {
        name = "default"
    }
    ext = 'pdf'
    if (l == 1 && t == 'list') {
        ext = 'xlsx'
        fpath = file.path(outputdir, paste0(name , '.' , ext))

        # save list of data.frames into excel
        #         file <- file.path('report',"data_titv.xlsx")
        # wb <- openxlsx::createWorkbook()
        #if (names(fc) %>% is.null)
        seqnames <-
            paste0("Sheet", seq_along(fc)) # or names(datas) if provided
        # fill '' names
        sheetnames = c()
        N = names(fc)
        for (i in seq_along(fc)) {
            n = N[[i]]
            if (n == "") {
                sheetnames = sheetnames %>% append(seqnames[[i]])
            } else{
                sheetnames = sheetnames %>% append(n) %>% str_sub(1, 25)
            }
        }
        sheetnames = sheetnames %>% make.unique

        wb = openxlsx::createWorkbook()
        for (i in seq_along(fc)) {
            sheetname = sheetnames[[i]]
            openxlsx::addWorksheet(wb, sheetname)
            tt = class(fc[[i]])
            title_line = table_title_lines[i]
            rowNames = yhas_rownames(fc[[i]])
            if (tt %>% intersect(c('tibble', 'data.frame', 'matrix')) %>% length > 0) {
                if (!is.null(title_line) && title_line != "") {
                    openxlsx::writeData(
                        wb,
                        sheet =  sheetname,
                        x =  paste0('## ', title_line),
                        startRow = 1,
                        startCol = 1
                    )
                    openxlsx::writeData(
                        wb,
                        sheet =  sheetname,
                        x =  fc[[i]],
                        startRow = 2,
                        startCol = 1,
                        rowNames = rowNames
                    )
                } else{
                    openxlsx::writeData(
                        wb,
                        sheet =  sheetname,
                        x =  fc[[i]],
                        startRow = 2,
                        startCol = 1,
                        rowNames = rowNames
                    )
                }
            } else{
                ylog('Skip [',
                     tt,
                     ']',
                     sheetname,
                     ', cause it is in an unsupport format')
            }
        }
        openxlsx::saveWorkbook(wb, file = fpath, overwrite = TRUE)

    } else if (l == 1 && t == "call") {
        ext = 'pdf'
        fpath = file.path(outputdir, paste0(name , '.', ext))

        tryCatch(expr = {
            grDevices::cairo_pdf(fpath, width = WIDTH, height = HEIGHT)
            eval(fc, envir = .GlobalEnv)
        }, finally = {
            dev.off()
        })
    } else if (c('ggplot', 'pheatmap', 'Heatmap') %>% intersect(t) %>% length > 0) {
        ext = 'pdf'
        fpath = file.path(outputdir, paste0(name , '.', ext))

        tryCatch(expr = {
            grDevices::cairo_pdf(
                fpath,
                width = WIDTH,
                height = HEIGHT,
                family = family
            )
            print(fc)
        }, finally = {
            dev.off()
        })
    } else if (c('tibble', 'data.frame', 'matrix') %>% intersect(t) %>% length > 0) {
        ext = 'xlsx'
        fpath = file.path(outputdir, paste0(name , '.', ext))
        title_line = table_title_lines[[1]]
        rowNames = yhas_rownames(fc)
        wb = openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, 'Sheet1')
        if (!is.null(title_line) && title_line != "") {
            openxlsx::writeData(
                wb,
                sheet =  'Sheet1',
                x =  paste0('## ', title_line),
                startRow = 1,
                startCol = 1
            )
            openxlsx::writeData(
                wb,
                sheet =  'Sheet1',
                x =  fc,
                startRow = 2,
                startCol = 1,
                rowNames = rowNames
            )
        } else{
            openxlsx::writeData(
                wb,
                sheet =  'Sheet1',
                x =  fc,
                startRow = 2,
                startCol = 1,
                rowNames = rowNames
            )
        }
        openxlsx::saveWorkbook(wb, file = fpath, overwrite = TRUE)
    }
    ylog("write 1", ext, "at", fpath, verbose = TRUE)
}



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
#' @importFrom grDevices dev.off
#' @export
#'
#' @examples
ydumpto = function(x,
                   fname = NULL,
                   export = NULL,
                   outputdir = OUTPUTROOT,
                   ext = NULL,
                   prefix = NULL,
                   flag = NULL,
                   worker = NULL,
                   verbose = T,
                   suffix = NULL,
                   mkdir = F,
                   ...) {
    # OUTER_USAGE = c('fname',"fname",'export',
    #                 'outputdir','ext','prefix',
    #                 'flag','worker','verbose',
    #                 'suffix','mkdir')

    if (is.null(fname) &&
        is.null(export))
        stop('fname export both NULL')
    if (!is.null(export)) {
        if (export == TRUE) {
            fname = 'plot'
        } else if (export == FALSE) {
            return("export==FALSE, exporting aborted")
        } else if (typeof(export) == 'character') {
            fname = export %>% str_flatten()
        } else {
            stop('export must be boolean or string')
        }
    }

    # fname contains path '/'
    if (fname %>% str_detect('/')) {
        tmp = fname %>% ysplit_path()
        if (length(tmp) > 1) {
            outputdir = tmp[[1]]
            fname = tmp[[2]]
        } else{
            fname = tmp[[2]]
        }
    }

    # fname contains ext i.e. '.'
    if (fname %>% str_detect('.')) {
        tmp = fname %>% ysplit_file_name(mode = 'auto', ext = ext)
        if (length(tmp) > 1) {
            fname = tmp[[1]]
            ext = tmp[[2]]
        } else{
            fname = tmp[[1]]
            ext = NULL
        }
    }

    # labels order : [worker,prefix,tag,export(fname),suffix]
    fname = c(worker, flag, prefix, fname, suffix)  %>%
        str_flatten(collapse = '_') %>%
        str_replace_all(pattern = fixed('__'), '_') %>%
        str_remove_all(pattern = '^_+|_+$')

    if ((!is.null(ext)) && ext == 'json') {
        attr(x, 'class') = c('json', class(x))
    }
    # out is generated
    #return(list(fname=fname,ext=ext,outputdir=outputdir,verbose=verbose))
    # if (class(x)=='function') .ydumpto(x,fname,ext,outputdir,verbose=verbose,args=args,...)
    if (mkdir == TRUE && !file.exists(outputdir)) {
        dir.create(outputdir, recursive = TRUE, mode = '0755')
    } else if (mkdir == FALSE && !file.exists(outputdir)) {
        ylog(outputdir, 'does not exist.', .addtime = FALSE)
        stop('!')
    }
    argv = yget_args(...)
    if (is.null(x)) {
        return(argv)
    } else{
        .ydumpto %>%  do.call(argv)
    }
}

# -----------utils-----------------------------------



#' yfile_path
#'
#' join given list of folders, likely the R version of os.path.join
#' may cause some unpredictable return values at certain circumstances
#'
#' @param ... parameters of path to join
#'
#' @return character represent the full path
#' @export
#'
#' @examples
yfile_path = function(...){
    args = list(...)
    args = as.character(args) %>% str_replace('\\\\','/')
    res = paste(args, sep = "/", collapse = "/")
    res = gsub("/[.]{0,1}/", "/", res)
    #     res = gsub("/[\\\\/]{2}/", "/", res)
    res
}

#' Extract parameters for a given function
#'
#' Get the arguments whose name match the given function(.f)'s formal arg names in the input args lists
#'
#' @param ... filter all the args by the formals of .f, return all matched,
#' @param .f is.function() and is the filter function, the returned args list will match the formals of the .f
#'
#' @return like `list(...)[intersect(list(...)%>%names,formals(.f)%>%names)]`
#' @export
#' @example
#' data(mtcars)
#' yget_args(data=mtcars,a='a',b='b',.f=ggplot) # list(data=mtcars)
yget_args = function (..., .f = NULL) {
    # if ((!is.null(.filter))&&(.f %>% is.null)){
    #     .f = .filter
    # }
    dots = list(...)
    pf <- parent.frame()
    pf = pf %>% as.list %>% utils::modifyList(dots)
    if (!is.null(.f)) {
        formal_args = formals(.f) %>% names
        pf = pf %>% ysubset_list_named(formal_args)
    }
    pf
}

#' Get the formal args of a function name, Only used inside a function
#'
#' @param func c(character,NULL,symbol), character/symbol represents a function, if NULL get the formal args of the
#'   parent function
#' @param get_default_values FALSE, only used when func is NULL, also get the default values of a given function?
#'
#' @return character vector if get_default_values==FALSE NOT recommended to use, cause can not determine the exact
#'   behavior at current at 2023-04-19 seems to get the list with names are the arg names and values are the default
#'   values at the run time
#' @export
#'
#' @examples
#' \dontrun{
#'   yfunc_args(utils::read.table)
#' }
#'
yfunc_args <- function(func=NULL, get_default_values = FALSE) {
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
    hasdots = '...' %in% fnames
    fnames <- fnames[fnames != '...']

    if (!is.null(func)) return(fnames)

    # Get currently set values for named variables in the parent frame
    args <- evalq(as.list(environment()), envir = parent.frame())
    if (hasdots){
        # Get the list of variables defined in '...'
        args <- c(args[fnames], evalq(list(...), envir = parent.frame()))
    }
    if(get_default_values) {
        # get default values
        defargs <- as.list(parent_formals)
        defargs <- defargs[unlist(lapply(defargs, FUN = function(x) class(x) != "name"))]
        args[names(defargs)] <- defargs
        setargs <- evalq(as.list(match.call())[-1], envir = parent.frame())
        args[names(setargs)] <- setargs
    }
    return(args)
}




#' Split a file name string into fname, extName
#'
#' @param str
#' @param ext set extention name of
#' @param mode in ['auto','keep','add','replace'] default 'auto', behave, see example
#' @param ... diff depend on str and ext values
#' @param str
#' @param ext
#' @param mode
#' @param ...
#'
#' @return character of length 2
#' @export
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
#'
ysplit_file_name = function(str, ext = NULL, mode = 'auto', ...) {
    if (!is.null(ext)) {
        # 判断ext首字符是不是点“.” 有点则去掉
        if (ext %>% str_sub(1, 1) == '.')
            ext = ext %>% str_sub(2, )
        .count = str_count(str, fixed('.'))
        if (.count == 0)
            res = c(str)
        else if (.count >= 1)
            res = str_split(str, fixed('.'))[[1]]
    } else{
        res = str_split(str, fixed('.'))[[1]]
    }
    l = length(res)
    if (l == 1) {
        # 'abc'
        file_name = res[[1]]
        ext_name = NULL
    } else if (l == 2) {
        # 'abc.txt'
        file_name = res[[1]]
        ext_name = res[[2]]
    } else if (l > 2) {
        # 'a.b.c.txt'
        file_name = res %>% yslice(1:-2) %>% str_flatten(collapse = '.')
        ext_name = res %>% yslice(-1)
    }
    if (mode == 'auto') {
        if (ext |> is.null()) {
            res = c(file_name, ext_name)
        } else{
            res = c(file_name, ext)
        }
    } else if (mode == 'keep')
        res = c(file_name, ext_name)
    else if (mode == 'add')
        res = c(file_name, ext_name, ext)
    else if (mode == 'replace')
        res = c(file_name, ext)
    res
}



#' split a string into a list of dirnames
#'
#' if path start with /(root), the 1st element of returned list
#' is '/', multiple '/' like '//' will be treated as single '/'
#'
#' @param path
#'
#' @return
#' @export
#'
ysplit_path = function(path) {
    path.list = str_split(path, pattern = '/+')[[1]]
    if (str_sub(path, 1, 1) == '/') {
        path.list[[1]] = '/'
    }
    path.list
}

#' Slice a iterable like `[:]` in a python indexing style
#'
#' @description  slice a iterable like `[:]`, but use a python-like style, negative-n means the
#'   last-n-th(python style) but not drop the n-th(R style)
#' @details if .style_negative_index == 'py' (default) negative values will be treated as to get the
#'   last N like python does except index is from 1 to length(@iterable) and the end of @seqs is
#'   included, too if .style_negative_index == 'r' (optional) negative values will be treated as to
#'   remove those at the positive position just like in r code []
#'
#' @param seqs
#' @param .style_negative_index
#' @param iterable vector or list
#'
#' @return sliced iterable
#' @export
#'
#' @examples
#' \dontrun{
#' a = c(0,1,2,3,'b')
#' a %>% yslice(-1) # 'b'
#' a %>% yslice(-3:-1) # c('2','3','b')
#' a %>% yslice(1:3) # c('0','1','2')
#' a %>% yslice(3) #  c('2','3','b')
#' }
yslice = function(iterable, seqs=NULL, .style_negative_index='py'){
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



# generic function for export
.ydumpto = function(x, ...) {
    UseMethod('.ydumpto')
}


#' @importFrom grDevices dev.off
.ydumpto.call = function(x,
                         fname,
                         ext = NULL,
                         outputdir = OUTPUTROOT,
                         verbose = TRUE,
                         row.names = FALSE,
                         ...) {
    # if (row.names %>% is.null) row.names= FALSE
    if (is.null(ext))
        ext = 'pdf'
    file = yfile_path(outputdir, paste(fname, ext, sep = '.'))

    if (ext %in% c('png', 'svg')) {
        dev = get(ext, .GlobalEnv)
    } else if (ext == 'pdf') {
        dev = get('grDevices::cairo_pdf', .GlobalEnv)
    }

    tryCatch({
        dev(file, ...)
        eval(x)
        ylog('[plot args]', args %>% names, ',', .addtime = F)
    }, finally = {
        dev.off()
    })

}


#' @importFrom grDevices dev.off
.ydumpto.data.frame = function(x,
                               fname,
                               ext = NULL,
                               outputdir = OUTPUTROOT,
                               verbose = TRUE,
                               ...) {
    # if (row.names %>% is.null) row.names= TRUE
    if (is.null(ext))
        ext = 'dfx'
    file = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    argv =  yget_args(..., .f = utils::write.table)
    if (argv$quote |> is.null())
        argv$quote = FALSE
    if (argv$sep |> is.null())
        argv$sep = '\t'
    #     print(c("[argv used]",names(argv)))
    utils::write.table %>% do.call(args = argv)
    ylog('write 1', ext, 'at', file, verbose = verbose)
    argv$file
}


#' @importFrom grDevices dev.off
.ydumpto.ggplot = function(x,
                           fname,
                           ext = NULL,
                           outputdir = OUTPUTROOT,
                           verbose = TRUE,
                           ...) {
    if (is.null(ext))
        ext = 'pdf'
    filename = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    if (ext %in% c('png', 'svg')) {
        dev = get(ext, .GlobalEnv)
    } else if (ext == 'pdf') {
        dev = get('grDevices::cairo_pdf', .GlobalEnv)
    }
    argv = yget_args(..., .filter = dev)
    tryCatch({
        dev |> do.call(argv)
        print(x)
    }, finally = {
        dev.off()
    })
    ylog('write 1', ext, 'at', filename, verbose = verbose)
}


#' @importFrom grDevices dev.off
.ydumpto.pheatmap = function(x,
                             fname,
                             ext = NULL,
                             outputdir = OUTPUTROOT,
                             verbose = TRUE,
                             ...) {
    if (is.null(ext))
        ext = 'pdf'
    #     filename = yfile_path(outputdir,paste(fname,ext,sep='.'))
    filename = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    if (ext %in% c('png', 'svg')) {
        dev = get(ext, .GlobalEnv)
    } else if (ext == 'pdf') {
        dev = get('grDevices::cairo_pdf', .GlobalEnv)
    }
    argv = yget_args(..., .filter = dev)
    tryCatch({
        dev |> do.call(argv)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
    }, finally = {
        dev.off()
    })
    ylog('write 1', ext, 'at', filename, verbose = verbose)
}

.ydumpto.matrix = function(x,
                           fname,
                           ext = NULL,
                           outputdir = OUTPUTROOT,
                           verbose = TRUE,
                           ...) {
    # if (row.names %>% is.null) row.names= TRUE
    if (is.null(ext))
        ext = 'dfx'
    file  = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    argv = yget_args(..., .filter = utils::write.table)
    path = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    if (argv$row.names %>% is.null())
        argv$row.names = TRUE
    if (argv$quote |> is.null())
        argv$quote = FALSE
    if (argv$sep |> is.null())
        argv$sep = '\t'
    ylog('[argv used]', names(argv), ',', .addtime = F)
    utils::write.table |> do.call(args = argv)
    ylog('write 1', ext, 'at', path, verbose = verbose)
}

.ydumpto.json = function(x,
                         fname,
                         ext = NULL,
                         outputdir = OUTPUTROOT,
                         verbose = TRUE,
                         ...) {
    if (is.null(ext))
        ext = 'json'
    path = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    argv = yget_args(..., .f = jsonlite::write_json)
    if (argv$auto_unbox |> is.null())
        argv$auto_unbox = TRUE
    ylog('[argv used]', names(argv), ',', .addtime = F)
    if (typeof(x) == 'character') {
        fileConn <- file(path)
        writeLines(x, fileConn)
        close(fileConn)
    } else{
        jsonlite::write_json |> do.call(args = argv)
    }
    ylog('write 1 json at', path, verbose = verbose)
}

# args is the list of args used by the plotting function
.ydumpto.function = function(x, args, verbose, ext, outputdir, fname, ...) {
    " return plot_func and its args used for plotting named 'plot_func' & 'plot_args' "
    argv = yget_args(...)
    if ('args' %in% (argv %>% names)) {
        print('⚠ run .ydumpto.function with no @args argument supplied')
        args = argv %>% ygetlast('args', .squeeze = T)
        argv = argv %>% yrmlast(c('args', 'x'))
    } else{
        args = NULL
        argv = argv %>% yrmlast(c('x'))
    }

    if (verbose == TRUE)
        ylog('you are dumping a function, assuming its a plotting one',
             .addtime = F)
    if (is.null(ext)) {
        ext = 'pdf'
    } else if (nchar(as.character(ext)) > 20) {
        print(
            paste0(
                'There maybe somthing wrong with <ext>, it is to long: "',
                str_sub(as.character(ext), 1, 20),
                '..." total ',
                nchar(as.character(ext)),
                ' chars'
            )
        )
        stop('!')
    }
    filename = yfile_path(outputdir, paste(fname, ext, sep = '.'))
    argv$filename = filename
    if (ext %in% c('png', 'svg')) {
        dev = get(ext, .GlobalEnv)
    } else if (ext == 'pdf') {
        dev = get('grDevices::cairo_pdf', .GlobalEnv)
    }
    tryCatch({
        argv1 = argv %>% ygetlast(formals(dev) %>% names)
        ylog('[dev args]', argv1 %>% names, ',', .addtime = F)
        dev %>% do.call(args = argv1)
        ylog('[plot args]', args %>% names, ',', .addtime = F)
        res = x %>% do.call(args)
    }, finally = {
        dev.off()
    })
    ylog('write 1', ext, 'at', filename, verbose = verbose)
    list(plot_func = x, plot_args = args)
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



