#' Title
#'
#' @param name
#'
#' @return
#' @export
#'
#' @examples
ygen_set_color_scheme = function(name='wjw'){
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

#' Title
#'
#' @param df
#' @param palette
#' @param cols
#' @param .assign_global
#'
#' @return
#' @export
#'
#' @examples
yset_random_anno_colors = function(df, palette=NULL, cols = c('Clin_classification'), .assign_global=FALSE){
    if (is.null(palette)){
        palette = ggsci::pal_npg()(10)
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

#' @title ygen_anno_color
#' @description random colors for annotation columns (set `colors=NULL`)or set with a list of index of colors in palette
#'  used for generate COL_ANNO for maftools::oncoplot(annotationColor = COL_ANNO)
#'
#' @param df data.frame
#' @param palette color palette, provide HEX color string vectors or use pal_*()(MAX_COLOR_NUM) to get a palette
#'  if palette is NULL, use palette = pal_png()(10)
#' @param cols generate colors for these columns in  `df`
#' @param colors a list of index of colors in palette `colors` is a list with the same length of cols,
#' each element will act on the corresponding element of cols, if colors is NULL, random colors will be used, the funcion will invoke yset_random_anno_colors
#' for each element of colors, if  NA or NULL, use palette(1:length()) to generate colors,
#'  if 'random', use sample(palette, length()) to generate colors,
#'  if a integer vector, use palette$color to generate colors,
#' @param .autoNAcolor if TRUE, set color of NA value to '#BBBBBB'
#' @param .assign_global assign to global env or not, if colors is NULL, .assign_global will be set to FALSE
#'
#' @return list of colors
#' @export
#' @examples
#' df = data.frame(a = c('a','a','b'), b = c('a','c','c'), c = c('a','b','c'))
#' ygen_anno_color(df, cols = c('a','b','c'),)
ygen_anno_color = function(df, palette=NULL, colors=NULL, cols = c('Clin_classification'), .autoNAcolor = TRUE, .assign_global=TRUE){
    if (is.null(palette)){
        palette = ggsci::pal_npg()(10)
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
#'
#' @param ... allowed value are c(labels = FALSE,borders = FALSE,cex_label = ?)
#' @param col_anno
#'
#' @export
#'
show_anno_col = function(col_anno,...){
    res = c()
    l = col_anno %>% sapply(length) %>% max
    for (i in col_anno){
        res = c(res,i,rep('#FFFFFF',l-length(i)))
    }
    print(col_anno)
    scales::show_col(res,ncol = l,...)
}
