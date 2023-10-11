tryCatch(
    expr = {
        cairo_pdf('data/test.pdf',5,5)
        plot(c(1,2,3),c(4,5,6))
        PDF_DEVICE <<- 'cairo_pdf'
        # assign("PDF_DEVICE", 'cairo_pdf',envir = .GlobalEnv)
        dev.off()
    },
    warning = function(w){
        PDF_DEVICE <<- 'pdf'
        # assign("PDF_DEVICE", 'pdf',envir = .GlobalEnv)
    },finally = {
        tryCatch({ dev.off() }, error = function(e) {})
        if (file.exists("data/test.pdf"))
            file.remove('data/test.pdf')
        message('PDF_DEVICE is: ',PDF_DEVICE)
    }
)

# ------------- IO ---------------------

#' Load files into data.table or list of data.tables
#'
#' @description use xlsx/openxlsx to load xls, xlsx and use data.table::fread to read plain text
#'   format tables like tsv,csv,txt,dfx,etc.
#'
#' @details 依赖于 yslice, ypush, ylog, ysplit_path, ysplit_file_name, yfile_path, <stringr> if `li` &&
#'   `fextname` both NULL and `frm` is not NULL, read all files from `frm`. Depends on: INPUTROOT,
#'   yfunc_args, read.table, %>%, modifyList, yslice, str_flatten, str_detect, fixed, ypush,
#'   str_sub, ypath_join
#'
#' @param li A character vector, whose element <i> will be interpretated as the value of  dfx names,
#'   if <i> starts with c('/','~','../'), find `<i>.<fextname>` in the root ('/'), home ('~') or
#'   parent ('../') folder. if <i> is a file name with no path info (like 'abc.txt' or 'abc'), then
#'   file <i> will be search in the <frm> folder. if <i> contains a folder level, (like 'abc/tmp',
#'   'report/abc.txt'), ignore <frm> for the current <i>, search files in the <i> given folder.
#'   NOTE: If you want to refer to the files in the subfolder in the <frm> path, add './' in the
#'   beginning of the <i>: if `frm = 'export', li = 'folder/data.txt'`, search the data.txt in the
#'   current folder (`./folder/data.txt` is searched), `frm = 'export', li = './folder/data.txt'`,
#'   search data.txt in the <frm> folder (`./export/folder/data.txt` is searched). if input li is a
#'   data.frame, a matrix, or a data.table, it will be returned as it is. `li` Overwrites <pattern>
#'   settings, if `length(<li>) > 1` generate a named list, else generate the data.frame read.
#'   Returned value is determined further by <retmode>.
#' @param frm a folder.path, all the data will be read from this path
#' @param pattern a RegExpr to glob, return matched files, only acts when `li` is NULL
#' @param fextname basicly, <li> members do not include an extention name, fextname is used to set
#'   the <li> extention name. set fextname=NULL to only load exactly the name <li> provided.
#'   fextname不会修改已经显式声明的扩展名，如li=c('abx.txt'),fextname='dfx' 则读入仍是'abx.txt'而不是'abx.dfx'或'
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
#' @seealso \code{\link[data.table]{fread}}
#' @examples
#'
yload_dfx = function(li = NULL,
                     frm = INPUTROOT,
                     pattern = NULL,
                     fextname = NULL,
                     engine = 'auto',
                     worker = NULL,
                     retmode = 'local',
                     verbosing = TRUE,
                     sheet_name = 1,
                     ...) {
    # retmode in c('local','global')
    # fextname in c(str, NULL)
    # default parameters with custom one
    user_args = yfunc_args()
    if (engine=='auto'){
        ns = names(user_args)
        if (length(intersect(c("row.names",'comment.char','fileEncoding'),ns))>0){
            default_load_func = utils::read.table
        }else if (length(intersect(c("skip",'sep2','select','drop'),ns))>0){
            default_load_func = data.table::fread
        }else{
            default_load_func = data.table::fread
        }
    }else if (engine == 'fread'){
        default_load_func = data.table::fread
    }else if(engine == 'read.table'){
        default_load_func = utils::read.table
    }else{
        stopifnot(is.function(engine))
    }
    std_args = yfunc_args(default_load_func)
    argv = list(
        sep = "\t",
        quote = "",
        na.strings = "",
        header = TRUE,
        check.names = FALSE,
        strip.white = TRUE,
        fill = FALSE,
        blank.lines.skip = FALSE,
        stringsAsFactors = FALSE
    ) %>% utils::modifyList(user_args, keep.null = TRUE) %>% ysubset_list_named(std_args)

    if ( length(intersect(c('data.frame', 'data.table', 'matrix'), class(li))) > 0 )
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
                if (myfrm[[1]] %in% c('/', '~','..')) {
                    myfrm = myfrm %>% str_flatten(collapse = '/')
                } else if(myfrm[[1]] %in% c('.')){
                    myfrm = c(frm, myfrm) %>% str_flatten(collapse = '/') %>% str_replace(fixed('/./'),'/')
                } else{
                    myfrm = myfrm %>% str_flatten(collapse = '/')
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
                tmp = ysplit_file_name(ffname, fextname = fextname, mode = 'auto')
                fname = tmp[[1]]
                if (length(tmp) == 1) {
                    fextname = 'dfx'
                } else if (length(tmp) == 2) {
                    fextname = tmp[[2]]
                } else{
                    stop('wrong tmp length')
                }
                files = files %>% ypush(c(myfrm, fname, fextname))
            }
        }
    }

    res = list()
    for (file in files) {
        # files is list of ('path','file_name','ext_name')
        myfrm = file[[1]]
        vname = file[[2]]
        fextname = file[[3]]
        if (vname %>% str_sub(1, 3) %>% str_detect('[A-Z][0-9]\\_'))
            worker = NULL
        if (!is.null(worker)) {
            ffname = str_flatten(c(worker, '_', vname, '.', fextname))
            input = yfile_path(myfrm, ffname)
            if (!file.exists(input)) {
                ylog(
                    '[WARN] prefix with <worker> file not found, retry without <worker> one',
                    .addtime = FALSE
                )
                ffname = str_flatten(c(vname, '.', fextname))
            }
        } else{
            ffname = str_flatten(c(vname, '.', fextname))
        }

        input = yfile_path(myfrm, ffname)

        if (verbosing == TRUE)
            ylog("read [", ffname , "] as '", vname, "' from > ", myfrm)

        if (fextname == 'gmt') {
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
        } else if (fextname =='xls'){
            res[[vname]] <- xlsx::read.xlsx(input, sheetIndex = sheet_name,...)
        } else if (fextname =='xlsx') {
            res[[vname]] <- openxlsx::read.xlsx(input, sheet = sheet_name,...)
        } else {
            "fextname == {txt,dfx,...} "
            argv$file = input
            res[[vname]] = default_load_func %>% do.call(args = argv)
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

#' Alias to \code{yload_dfx}
#'
#' @export
#' @seealso \code{\link{yload_dfx}}
#'
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
                  table_title_lines = NULL,
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
    fextname = 'pdf'
    if (l == 1 && t == 'list') {
        fextname = 'xlsx'
        fpath = yfile_path(outputdir, paste0(name , '.' , fextname))

        # save list of data.frames into excel
        #         file <- yfile_path('report',"data_titv.xlsx")
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
            if (!is.null(table_title_lines)){
                title_line = table_title_lines[i]
            }else{
                title_line = NULL
            }

            rowNames = yhas_rownames(fc[[i]])
            if (tt %>% intersect(c('data.table','tibble', 'data.frame', 'matrix')) %>% length > 0) {
                sheet  = as.data.frame(fc[[i]])
                if (!is.null(title_line)) {
                    if (title_line != ""){
                            openxlsx::writeData(
                            wb,
                            sheet =  sheetname,
                            x =  paste0('## ', title_line),
                            startRow = 1,
                            startCol = 1
                        )
                    }else{
                        openxlsx::writeData(
                            wb,
                            sheet =  sheetname,
                            x =  "",
                            startRow = 1,
                            startCol = 1
                        )
                    }
                    openxlsx::writeData(
                        wb,
                        sheet =  sheetname,
                        x =  sheet,
                        startRow = 2,
                        startCol = 1,
                        rowNames = rowNames
                    )
                } else {
                    openxlsx::writeData(
                        wb,
                        sheet =  sheetname,
                        x =  sheet,
                        startRow = 1,
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
        fextname = 'pdf'
        fpath = yfile_path(outputdir, paste0(name , '.', fextname))

        tryCatch(expr = {
            grDevices::cairo_pdf(fpath, width = WIDTH, height = HEIGHT)
            eval(fc, envir = .GlobalEnv)
        }, finally = {
            dev.off()
        })
    } else if (c('ggplot', 'pheatmap', 'Heatmap') %>% intersect(t) %>% length > 0) {
        fextname = 'pdf'
        fpath = yfile_path(outputdir, paste0(name , '.', fextname))

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
        fextname = 'xlsx'
        fpath = yfile_path(outputdir, paste0(name , '.', fextname))
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
    ylog("write 1", fextname, "at", fpath, verbose = TRUE)
}



#' Dump objects into types of files depend on object class
#'
#' @param x
#' @param fname the file name to put,
#' @param outputdir the folder to put the files, default is the global var [OUTPUTROOT]
#' @param fextname NULL, if NULL or '', auto add extension name to the fname, else set the extension name
#'   to [fextname]
#' @param prefix, @param flag, @param worker, @param suffix tags add to the filename part of the
#'   output file
#' @param verbose TRUE
#' @param mkdir FALSE, if the destination dir does not exist, raise Error when FALSE, else create
#'   dirs to the destination dir
#' @param ...
#'
#' @details
#' Please specify fextname='txt' if your fname contains '.', or the last part of .* will be used
#' as extension names
#'
#' @return
#' @importFrom grDevices dev.off
#' @export
#'
#' @examples
ydumpto = function(x,
                   fname = NULL,
                   outputdir = OUTPUTROOT,
                   ext = NULL,
                   to_plain_txt = FALSE,
                   prefix = NULL,
                   flag = NULL,
                   worker = NULL,
                   suffix = NULL,
                   verbose = T,
                   mkdir = FALSE,
                   ...) {

    # assert non-Empty objs
    if (is.null(fname))
        stop('fname is NULL')
    if (is.null(x))
        stop('x is NULL')
    if (tibble::is_tibble(x)){
        x = as.data.frame(x)
    }
    # fpath is the final outputdir to use
    # ffname is the final filename to use
    # fextname is the final extentsion name to use

    # determine fpath base on input fname and input outputdir
    # if fname start with ../ or ./ or / or ~/ ignore outputdir else use outputdir



    if (stringr::str_starts(fname,stringr::fixed('./'))){
        fpath = fs::path_dir(fname)
    }else{
        if (basename(dirname(fname)) == ".") {
            # path contains only one level
            fpath = outputdir
        } else {
            # path contains more than one level
            fpath = fs::path_dir(fname)
        }
    }

    fname = fs::path_file(fname)
    ffname = tools::file_path_sans_ext(fname)

    # determine the fextname based on input fextname and input fname
    # fextname will be '' if fextname is NULL and not provide fname
    if (!is.null(ext))
        fextname = ext
    else
        fextname = fs::path_ext(fname)


    # add worker, flag, prefix, suffix labels to ffname, trim additional '_'
    ffname = c(worker, flag, prefix, ffname, suffix)  %>%
        str_flatten(collapse = '_') %>%
        str_replace_all(pattern = fixed('__'), '_') %>%
        str_remove_all(pattern = '^_+|_+$')

    if (to_plain_txt == FALSE){
        if (fextname == 'json') {
            attr(x, 'class') = c('json', class(x))
        }else if (fextname == 'RDS' && is.list(x)){
            attr(x,'class') = c('RDS', class(x))
        } else if (fextname == 'xlsx' && is.list(x)){
            attr(x, 'class') = c('excelObj', class(x))
        }
    }else{
        attr(x,'class') = c('PlainTXT',class(x))
    }


    if (mkdir == TRUE && !file.exists(outputdir)) {
        dir.create(outputdir, recursive = TRUE, mode = '0644')
    } else if (mkdir == FALSE && !file.exists(outputdir)) {
        ylog(outputdir, 'does not exist.', .addtime = FALSE)
        stop('!')
    }

    .ydumpto(x,fpath=fpath,ffname=ffname,fextname=fextname,verbose=verbose,...)

}

# -------------------------- utils --------------------------------



#' Return nice format path values of file.path
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
    p = file.path(...)
    p %>% str_replace(fixed('\\'),'/') %>% str_replace(fixed('//'),'/') %>% str_replace(fixed('/./'),'/')
    # args = list(...)
    # args = as.character(args) %>% str_replace('\\\\','/')
    # res = paste(args, sep = "/", collapse = "/")
    # res = gsub("/[.]{0,1}/", "/", res)
}

#' Extract vars in the parent frame, then filtered vars by a given function's formals
#'
#' Get the arguments whose name match the given function(.f)'s formal arg names in the input args
#' lists
#'
#' @param ... filter all the args by the formals of .f, return all matched,
#' @param .f .filter, is.function() and is the filter function, the returned args list will match
#'   the formals of the .f
#'
#' @return like `list(...)[intersect(list(...)%>%names,formals(.f)%>%names)]`
#' @export
#' @example data(mtcars) yget_args(data=mtcars,a='a',b='b',.f=ggplot) # list(data=mtcars)
yget_args = function (..., .f = NULL) {
    stopifnot(is.null(.f) || is.function(.f))
    dots = list(...)
    # get all the vars including arguments and vars define in the parent function
    pf = parent.frame()
    pf = pf %>% as.list %>% utils::modifyList(dots)
    if (!is.null(.f)) {
        formal_args = formals(.f) %>% names
        pf = ysubset_list_named(pf,formal_args)
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
#' @examples
#' \dontrun{
#'   yfunc_args(utils::read.table)
#' }
#'
yfunc_args = function(func=NULL, get_default_values = FALSE) {
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
#' @param fextname set extention name of
#' @param mode in ['auto','keep','add','replace'] default 'auto', behave, see example
#' @param ... diff depend on str and fextname values
#' @param str
#' @param fextname
#' @param mode
#' @param ...
#'
#' @return character of length 2
#'
#' @example
#' ysplit_file_name('a.b.c.txt',fextname=NULL,mode='auto') # c('a.b.c', 'txt')
#' ysplit_file_name('abc',fextname='dfx',mode='auto') # c('abc', 'dfx')
#' ysplit_file_name('abc.txt',fextname='dfx',mode='auto') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',fextname='dfx',mode='auto') # c('a.b.c', 'dfx')
#'
#' ysplit_file_name('a.b.c.txt',fextname=NULL,mode='replace') # c('a.b.c')
#' ysplit_file_name('abc',fextname='dfx',mode='replace') # c('abc', 'dfx')
#' ysplit_file_name('abc.txt',fextname='dfx',mode='replace') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',fextname='dfx',mode='replace') # c('a.b.c', 'dfx')
#'
#' ysplit_file_name('a.b.c.txt',fextname=NULL,mode='keep') # c('a.b.c', 'txt')
#' ysplit_file_name('abc',fextname='dfx',mode='keep') # c('abc')
#' ysplit_file_name('abc.txt',fextname='dfx',mode='keep') # c('abc', 'txt')
#' ysplit_file_name('a.b.c.txt',fextname='dfx',mode='keep') # c('a.b.c', 'txt')
#'
#' ysplit_file_name('abc',fextname='dfx',mode='add') # c('abc', 'dfx')
#' ysplit_file_name('a.b.c.txt',fextname=NULL,mode='add') # c('a.b.c', 'txt')
#' ysplit_file_name('abc.txt',fextname='dfx',mode='add') # c('abc', 'txt', 'dfx')
#' ysplit_file_name('a.b.c.txt',fextname='dfx',mode='add') # c('a.b.c', 'txt', 'dfx')
#'
ysplit_file_name = function(str, fextname = NULL, mode = 'auto', ...) {
    if (!is.null(fextname)) {
        # 判断fextname首字符是不是点“.” 有点则去掉
        if (fextname %>% str_sub(1, 1) == '.')
            fextname = fextname %>% str_sub(2, )
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
        if (fextname %>% is.null()) {
            res = c(file_name, ext_name)
        } else{
            res = c(file_name, fextname)
        }
    } else if (mode == 'keep')
        res = c(file_name, ext_name)
    else if (mode == 'add')
        res = c(file_name, ext_name, fextname)
    else if (mode == 'replace')
        res = c(file_name, fextname)
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

.ydumpto.call = function(x,
                         fpath,
                         ffname,
                         fextname,
                         verbose=TRUE,
                         ...) {
    # if (row.names %>% is.null) row.names= FALSE
    if (fextname == '')
        fextname = 'pdf'
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
    }

    tryCatch({
        dev(file, ...)
        eval(x)
    }, finally = {
        dev.off()
    })
    ylog('write 1 <', fextname, '> at', file, verbose = verbose)
    file
}


.ydumpto.data.frame = function(x,
                               fpath,
                               ffname,
                               fextname,
                               verbose=TRUE,
                               ...) {
    # if (row.names %>% is.null) row.names= TRUE

    if (fextname == '')
        fextname = 'dfx'
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))

    argv = yget_args(..., .f = utils::write.table)
    if (argv$quote %>% is.null())
        argv$quote = FALSE
    if (argv$sep %>% is.null())
        argv$sep = '\t'
    if (argv$row.names %>% is.null())
        argv$row.names = yhas_rownames(x)

    #     print(c("[argv used]",names(argv)))
    utils::write.table %>% do.call(args = argv)
    ylog('write 1 <', fextname, '> at', file, verbose = verbose)
    argv$file
}

.ydumpto.matrix = function(x,
                           fpath,
                           ffname,
                           fextname,
                           verbose=TRUE,
                           ...) {
    if (fextname == '')
        fextname = 'dfx'
    file  = yfile_path(fpath, paste(ffname, fextname, sep = '.'))

    argv = yget_args(..., .f = utils::write.table)
    if (argv$row.names %>% is.null())
        argv$row.names = !(x %>% rownames %>% is.null)
    if (argv$quote %>% is.null())
        argv$quote = FALSE
    if (argv$sep %>% is.null())
        argv$sep = '\t'

    utils::write.table %>% do.call(args = argv)
    ylog('write 1 <', fextname, '> at', file, verbose = verbose)
    argv$file
}

.ydumpto.gg = function(x,
                           fpath,
                           ffname,
                           fextname,
                           verbose=TRUE,
                           ...) {
    if (fextname == '')
        fextname = 'pdf'
    # used by pdf
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    # used by cairo_pdf, png, svg
    filename = file

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
    }

    argv = yget_args(..., .f = dev)
    if (argv$width %>% is.null)
        argv$width = get("WIDTH",envir = .GlobalEnv)
    if (argv$height %>% is.null)
        argv$height = get("HEIGHT",envir = .GlobalEnv)

    tryCatch({
        dev %>% do.call(argv)
        print(x)
    }, finally = {
        dev.off()
    })
    ylog('write 1 <', fextname,'> in [',argv$width,',',
         argv$height, '] inches at', filename, verbose = verbose)
    argv[[1]]
}

.ydumpto.ggsurvplot = function(x,
                       fpath,
                       ffname,
                       fextname,
                       verbose=TRUE,
                       ...) {
    if (fextname == '')
        fextname = 'pdf'
    # used by pdf
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    # used by cairo_pdf, png, svg
    filename = file

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
    }

    argv = yget_args(..., .f = dev)
    if (argv$width %>% is.null)
        argv$width = get("WIDTH",envir = .GlobalEnv)
    if (argv$height %>% is.null)
        argv$height = get("HEIGHT",envir = .GlobalEnv)

    tryCatch({
        dev %>% do.call(argv)
        print(x)
    }, finally = {
        dev.off()
    })
    ylog('write 1 <', fextname,'> in [',argv$width,',',
         argv$height, '] inches at', filename, verbose = verbose)
    argv[[1]]
}

.ydumpto.Heatmap = function(x,
                            fpath,
                            ffname,
                            fextname,
                            verbose=TRUE,
                            ...){
    if (fextname == '')
        fextname = 'pdf'

    # used by pdf
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    # used by cairo_pdf, png, svg
    filename = file

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
    }

    argv = yget_args(..., .f = dev)
    if (argv$width %>% is.null)
        argv$width = get("WIDTH",envir = .GlobalEnv)
    if (argv$height %>% is.null)
        argv$height = get("HEIGHT",envir = .GlobalEnv)

    tryCatch({
        dev %>% do.call(argv)
        ComplexHeatmap::draw(x,...)
    }, finally = {
        dev.off()
    })
    ylog('write 1 <', fextname,'> in [',argv$width,',',
         argv$height, '] inches at', filename, verbose = verbose)

}



.ydumpto.pheatmap = function(x,
                             fpath,
                             ffname,
                             fextname,
                             verbose=TRUE,
                             ...) {
    if (fextname == '')
        fextname = 'pdf'

    # used by pdf
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    # used by cairo_pdf, png, svg
    filename = file

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
    }

    argv = yget_args(..., .f = dev)
    if (argv$width %>% is.null)
        argv$width = get("WIDTH",envir = .GlobalEnv)
    if (argv$height %>% is.null)
        argv$height = get("HEIGHT",envir = .GlobalEnv)

    tryCatch({
        dev %>% do.call(argv)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
    }, finally = {
        dev.off()
    })
    ylog('write 1 <', fextname,'> in [',argv$width,',',
         argv$height, '] inches at', filename, verbose = verbose)
    argv[[1]]
}

.ydumpto.PlainTXT = function(x,
                             fpath,
                             ffname,
                             fextname,
                             verbose=TRUE,
                             ...){
    if (fextname == '')
        fextname = 'txt'

    path = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    tryCatch(expr = {
        fileConn <- file(path)
        writeLines(x, fileConn)
    },finally = {
        close(fileConn)
    })
    ylog('write 1 plain text <',fextname, '> at', path, verbose = verbose)
}

.ydumpto.json = function(x,
                         fpath,
                         ffname,
                         fextname,
                         verbose=TRUE,
                         ...) {
    if (fextname == '')
        fextname = 'json'
    path = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    argv = yget_args(..., .f = jsonlite::write_json)
    if (argv$auto_unbox %>% is.null())
        argv$auto_unbox = TRUE
    if (typeof(x) == 'character') {
        fileConn <- file(path)
        writeLines(x, fileConn)
        close(fileConn)
    } else{
        jsonlite::write_json %>% do.call(args = argv)
    }
    ylog('write 1 <json> at', path, verbose = verbose)
}


.ydumpto.RDS = function(x,
                        fpath,
                        ffname,
                        fextname,
                        verbose=TRUE,
                        ...){
    if (fextname == '')
        fextname = 'RDS'

    path = yfile_path(fpath, paste(ffname, fextname, sep = '.'))

    saveRDS(x,path)

    ylog('write 1 <',fextname, '> at', path, verbose = verbose)
}

#
#' Title
#'
#' return plot_func and its args used for plotting named 'plot_func' & 'plot_args'
#' args is the list of args used by the plotting function
#'
#' @param x
#' @param fpath
#' @param ffname
#' @param fextname
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
.ydumpto.function = function(x,
                             args,
                             fpath,
                             ffname,
                             fextname,
                             verbose=TRUE,
                             ...) {
    "  "
    if (fextname == '')
        fextname = 'pdf'
    # used by pdf
    file = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    # used by cairo_pdf, png, svg
    filename = file

    if (fextname %in% c('png', 'svg')) {
        dev = get(fextname, asNamespace('grDevices'))
        argv = yget_args(..., .f = dev)
        argv$filename = filename
    } else if (fextname == 'pdf') {
        dev = get(PDF_DEVICE, asNamespace('grDevices'))
        argv = yget_args(..., .f = dev)
        if (PDF_DEVICE=='cairo_pdf')
            argv$filename = filename
        else if (PDF_DEVICE =='pdf')
            argv$file = filename
    }


    if (argv$width %>% is.null)
        argv$width = get("WIDTH",envir = .GlobalEnv)
    if (argv$height %>% is.null)
        argv$height = get("HEIGHT",envir = .GlobalEnv)

    tryCatch({
        dev %>% do.call(args = argv)
        res = x %>% do.call(args)
    }, finally = {
        dev.off()
    })

    ylog('write 1 <', fextname,'> in [',argv$width,',',
         argv$height, '] inches at', filename, verbose = verbose)
    # argv = yget_args(...)
    # if ('args' %in% (argv %>% names)) {
    #     print('⚠ run .ydumpto.function with no @args argument supplied')
    #     args = argv %>% ygetlast('args', .squeeze = T)
    #     argv = argv %>% yrmlast(c('args', 'x'))
    # } else{
    #     args = NULL
    #     argv = argv %>% yrmlast(c('x'))
    # }
    #
    # if (verbose == TRUE)
    #     ylog('you are dumping a function, assuming its a plotting one',
    #          .addtime = F)
    # if (fextname=='') {
    #     fextname = 'pdf'
    #
    # } else if (nchar(as.character(fextname)) > 20) {
    #     print(
    #         paste0(
    #             'There maybe somthing wrong with <fextname>, it is to long: "',
    #             str_sub(as.character(fextname), 1, 20),
    #             '..." total ',
    #             nchar(as.character(fextname)),
    #             ' chars'
    #         )
    #     )
    #     stop('!')
    # }
    # filename = yfile_path(fpath, paste(ffname, fextname, sep = '.'))
    #
    #
    #
    # ylog('write 1', fextname, 'at', filename, verbose = verbose)
    filename
}






