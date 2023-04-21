#' Set the size for the following plots
#' @description  basically, this function affects the following plot function which will use global var WIDTH and HEIGHT
#'   to determine the display repr.plot.width and repr.plot.height, set global variables depend on the third parameter
#'   `.set_global`
#' @param w set global WIDHT,  in inch
#' @param h set global HEIGHT, in inch
#' @param .set_global TRUE, set global vars, else no global variables set
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param .set_global
#'
#' @return
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
#' @return
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
