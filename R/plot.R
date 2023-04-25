
#' Plot heatmap using ComplexHeatmap
#'
#' @param mat matrix, usually columns is tsb, if not matrix, will try to coerce using [as.matrix()]
#' @param colData None, or the data.frame with rownames (is tsb) must be accord with colnames(`mat`) and columns of
#'   colData will be draw at the top annotation of Heatmap, each column a line
#' @param annotation_colors define the colors used for draw top annotation, representing colData
#' @param scale default 'auto', use [yscale_rows()] to standardize each rows
#' @param by_group if colData provided
#' @param params_top_annotation plotting params to draw top annotation
#' @param cluster_rows pass to [ComplexHeatmap::Heatmap()]
#' @param cluster_columns pass to [ComplexHeatmap::Heatmap()]
#' @param clustering_method_columns pass to [ComplexHeatmap::Heatmap()]
#' @param clustering_method_rows pass to [ComplexHeatmap::Heatmap()]
#' @param name pass to [ComplexHeatmap::Heatmap()]
#' @param row_title pass to [ComplexHeatmap::Heatmap()]
#' @param show_row_dend pass to [ComplexHeatmap::Heatmap()]
#' @param show_column_dend pass to [ComplexHeatmap::Heatmap()]
#' @param show_row_names pass to [ComplexHeatmap::Heatmap()]
#' @param show_column_names pass to [ComplexHeatmap::Heatmap()]
#' @param color in format 'Num1~Num2',default '-2~2', will be converted to [circlize::colorRamp()] function then pass to
#'   [ComplexHeatmap::Heatmap()] as col. if the two Nums span 0, then 0 will be inserted into the middle
#'   ELSE you can provide a custom [circlize::colorRamp()] object
#' @param column_split pass to [ComplexHeatmap::Heatmap()]
#' @param heatmap_legend_param pass to [ComplexHeatmap::Heatmap()]
#' @param ... other legal params pass to [ComplexHeatmap::Heatmap()]
#'
#' @return list
#' @export
#'
#' @examples
yplot_heatmap = function(
        mat
        ,colData = NULL
        ,annotation_colors = NULL
        ,scale="auto"
        ,by_group = TRUE
        ,params_top_annotation = NULL

        ,cluster_rows = TRUE
        ,cluster_columns = TRUE
        ,clustering_method_columns = 'ward.D2'
        ,clustering_method_rows = 'ward.D2'
        ,name = 'Z-Score'
        ,row_title=NA
        ,show_row_dend = TRUE
        ,show_column_dend = TRUE
        ,show_row_names = TRUE
        ,show_column_names = TRUE
        ,color = '-2~2'
        ,column_split = NULL
        ,heatmap_legend_param = NULL
        ,...){

    if (! 'matrix' %in% class(mat)){
        mat %<>% as.matrix
    }

    if (is.null(colData)){
        .ha = NULL
        if (by_group == TRUE){
            by_group = FALSE
            print('colData is None, set by_group = FALSE')
        }
    }else{
        if (! 'data.frame' %in% class(colData)){
            colData %<>% as.data.frame
        }

        .text_width = ComplexHeatmap::max_text_width(
            rownames(mat),
            gp = grid::gpar(fontsize = 12)
        )
        cn = colnames(colData)

        if (!yhas_rownames(colData)){
            can = c('Custom_Label','Tumor_Sample_Barcode','tsb','a','ap')
            x = intersect(can,cn)
            if (length(x)>0) {
                sel = x[1]
            } else stop('One of `',str_flatten_comma(can),'` must in colnames(colData)')
            colData %<>% dplyr::column_to_rownames(sel)
        }

        colData %<>% mutate(across(!where(is.factor),factor))

        .group = ydeframe(colData)
        if (!is.null(annotation_colors) && !is.na(annotation_colors)){
            .group$col = annotation_colors
        }
        if (!is.null(params_top_annotation)){
            stopifnot(is.list(params_top_annotation))
            .group = utils::modifyList(.group,params_top_annotation)
        }
        .ha = do.call(ComplexHeatmap::HeatmapAnnotation, args = .group)
    }

    if (by_group==TRUE){
        .column_split = .group[[1]]
    }else{
        .column_split = column_split
    }

    .text_width = ComplexHeatmap::max_text_width(
        rownames(mat),
        gp = grid::gpar(fontsize = 12)
    )

    if (is.null(scale) || is.na(scale)){
        # pass
        use_scale = FALSE
    }else{
        if (mat %>%
            apply(1,stats::sd) %>%
            bazar::almost.equal (1) %>%
            all &&
            mat %>%
            apply(1,mean) %>%
            bazar::almost.equal (0) %>%
            all){
            print('mat is already z-scored')
            use_scale = FALSE
        }else{
            use_scale = TRUE
        }
    }
    if (use_scale){
        if (scale == TRUE || scale == 'auto'){
            mat %<>% yscale_rows
            print(paste0("scale == ",scale,", use yscale_rows"))
        }else if (is.function(scale)){
            mat %<>% apply(1, scale)
            print('scale == TRUE, use apply(1, scale)')
        }else if(scale == 'none' || scale == FALSE){
            # pass
        }
    }
    if (is.null(heatmap_legend_param)){
        heatmap_legend_param = list(title=name)
    }
    if ('color' %in% names(heatmap_legend_param)){
        color = heatmap_legend_param$color
        heatmap_legend_param$color = NULL
    }

    if (is.character(color)){
        if (color == '-2~2' || color == 'default'){
            color = circlize::colorRamp2(breaks = c(-2, 0, 2),colors =  c("blue", "white", "red"), space = 'RGB')
        }else if (str_detect(color,'~')){
            tryCatch(expr = {
                x = str_split(color,'~')[[1]] %>% as.numeric
            },warning = function(w){
                stop(w)
            },error=function(e){
                stop(e)
            }
            )
            if (length(x)<=3){
                if (x[[1]] * x[[2]] < 0){
                    x = c(x[[1]],0,x[[2]])
                }
                color = circlize::colorRamp2(breaks = x, colors =  c("blue", "white", "red")[1:length(x)], space = 'RGB')
            }else{
                color = circlize::colorRamp2(breaks = x, colors = rainbow(length(x)), space = 'RGB')
            }
        }else{
            # pass
        }
    }else if (is.function(color)){

    }else{
        stop('color must be function or a color vector')
    }
    hm = ComplexHeatmap::Heatmap(mat
                                 ,col = color
                                 ,name = name
                                 ,clustering_method_columns = clustering_method_columns
                                 ,clustering_method_rows = clustering_method_rows
                                 # ,row_names_side = 'left'

                                 ,row_title = row_title
                                 ,show_row_dend = show_row_dend
                                 ,show_column_dend = show_column_dend
                                 ,show_row_names = show_row_names
                                 ,show_column_names = show_column_names

                                 ,heatmap_legend_param = heatmap_legend_param
                                 ,top_annotation = .ha
                                 ,column_split = .column_split
                                 ,row_names_max_width = .text_width
                                 ,cluster_rows = cluster_rows
                                 ,cluster_columns = cluster_columns
                                 # ,row_split = 4
                                 # ,cluster_columns = cluster_within_group(mat, group)
                                 ,...
    )
    gg.hm = ggplotify::as.ggplot(hm)

    list(gg=gg.hm
         ,heatmap=hm
         ,mat=mat
         ,colData=colData
         ,HeatmapAnnotation=.ha
         ,params=list(mat
                      ,clustering_method_columns = clustering_method_columns
                      ,clustering_method_rows = clustering_method_rows
                      ,name = name
                      ,color = color
                      ,row_title = row_title
                      ,show_row_dend = show_row_dend
                      ,show_column_dend = show_column_dend
                      ,show_row_names = show_row_names
                      ,show_column_names = show_column_names
                      ,top_annotation = .ha
                      ,column_split = .column_split
                      ,row_names_max_width = .text_width
                      ,heatmap_legend_param = heatmap_legend_param
                      ,...)
    )
}



#' Plot volcano plot
#'
#' @param diff_expr RNA diff expr
#' @param value_var column name which used for filter
#' @param threshold vector of length 2, the log2FoldChange threshold
#' @param p
#' @param add_label whether to add the label text of genes
#' @param ...
#'
#' @return list(gg=ggplot,de=data.frame)
#' @export
#'
#' @examples
yplot_volcano_using_de = function(diff_expr, value_var = 'padj', threshold=c(-2,2),p=0.05,add_label=FALSE,...){
    if (diff_expr %>% is.data.frame){
        de = diff_expr
    }else{
        de = diff_expr %>% data.frame
    }
    de$diffexpressed <- "NO"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    de$diffexpressed[de$log2FoldChange > threshold[2] & de$padj < p] <- "UP"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    de$diffexpressed[de$log2FoldChange < threshold[1] & de$padj < p] <- "DOWN"
    de$diffexpressed = factor(de$diffexpressed,levels=c('UP','NO','DOWN'))
    de$y = -log10(de[[value_var]])
    if (add_label){
        library(ggrepel)
        de$delabel <- NA
        de$gene_symbol = de %>% rownames
        de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
        g <- ggplot(data=de, aes(x=log2FoldChange, y=y, col=diffexpressed, label=delabel))+
            geom_point(alpha=0.7)+
            geom_text_repel(max.iter=10, box.padding = 1)
    }else{
        g <- ggplot(data=de, aes(x=log2FoldChange, y=y, col=diffexpressed))+
            geom_point(alpha=0.7)
    }
    if (value_var=='padj'){
        old_value_var = value_var
        value_var = 'FDR'
    }
    g = g + theme_minimal() +
        ylab(paste0("-log10(",value_var,")")) +
        scale_color_manual(values=c("red", "gray", "blue")) +
        geom_vline(xintercept=threshold, col="gray",linetype='longdash') +
        geom_hline(yintercept=-log10(p), col="gray",linetype='longdash')

    list(gg=g,de=de)
}

#' t-SNE analysis for mat
#'
#' @param normd data.fram or matrix, only columns will be preserved, rows will be ruduced in dimensions
#' @param colData
#' @param mat
#' @param group_col
#' @param perplexity
#' @param dims
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
yplot_tsne = function (normd = NULL, colData = NULL, mat = NULL, group_col = NA,
                       perplexity = NA, dims = 2, ...) {
    if (is.null(mat)) {
        mat = t(normd)
        pids = colnames(normd)
    }
    else {
        pids = rownames(mat)
    }
    if (is.na(perplexity)) {
        perplexity = floor((nrow(mat) - 1)/3)
    }
    stopifnot(!is.null(colData) && !is.null(mat))
    res = Rtsne::Rtsne(mat, perplexity = perplexity, dims = dims, ...)
    df = res$Y %>% data.frame(row.names = pids)
    colData2 = colData[pids, , drop = FALSE]
    df = cbind(df, colData2)
    res$df = df
    col = colnames(colData)
    if (is.na(group_col)) {
        if ("Group" %in% col)
            group_col = "Group"
        else if ("Clin_classification" %in% col)
            group_col = "Clin_classification"
        else if ("g" %in% col)
            group_col = "g"
        else {
            print("can not determine group_col in colData, stop plotting")
            return(res)
        }
    }
    group_col = sym(group_col)
    gg = df %>% ggplot(aes(x = X1, y = X2, color = !!group_col)) +
        geom_point(size = 3, alpha = 0.7) + ggtitle("t-SNE")
    res$gg = gg
    res
}


#' PCA analysis for mat
#'
#' @param normd data.fram or matrix, only columns will be preserved, rows will be ruduced in dimensions
#' @param colData
#' @param color_col
#' @param add_label
#' @param tag_fontsize pass to `size` in geom_label_repel
#' @param max.overlaps pass to ggrepel
#' @param alpha pass to geom_point
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
yplot_pca = function (normd, colData, color_col = "Group", add_label = FALSE,
                      tag_fontsize = 6, max.overlaps = Inf, alpha = 0.6, ...) {
    sig_table = normd %>% t
    mask = (sig_table %>% apply(sd, MARGIN = 2)) > 0
    mat = sig_table[, mask]
    if (ncol(mat) < ncol(sig_table))
        ylog("containing ", sum(!mask), " features with sd==0, removed")
    pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
    color_col = sym(color_col)
    data = cbind(pca$x[colData %>% rownames, 1:2], colData) %>%
        rownames_to_column("Tumor_Sample_Barcode")
    gg = data %>% ggplot(aes(PC1, PC2, color = !!color_col, label = Tumor_Sample_Barcode)) +
        geom_point(size = 3, alpha = alpha) + ggtitle("PCA")
    if (add_label == TRUE) {
        gg = gg + ggrepel::geom_label_repel(size = {
            {
                tag_fontsize
            }
        }, alpha = alpha, max.overlaps = max.overlaps)
    }
    list(gg = gg, pca = pca, data = data)
}

#' Plot venn diagram for 2~5 sets
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
yplot_venns = function(...){
    sets = list(...)
    x = length(sets)
    stopifnot(x <= 5&& x >=1)

    argv = list()
    DICT = alist(VennDiagram::draw.single.venn, VennDiagram::draw.pairwise.venn, VennDiagram::draw.triple.venn,
                VennDiagram::draw.quad.venn, VennDiagram::draw.quintuple.venn)
    plot_func = DICT[[x]]
    for (l in 1:length(sets)) {
        for (i in combn(1:length(sets), l, simplify = F)) {
            s = sets[i]
            if (l > 1) {
                r = Reduce(intersect, s) %>% length
                argv[[str_flatten(c("n", i))]] = r
            }
            else {
                r = length(s[[1]])
                argv[[str_flatten(c("area", i))]] = r
            }
        }
    }
    argv$cross.area = argv$n12
    argv = argv %>% ypush(list(category = sets %>% names, fill = c("dodgerblue",
                                                                   "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(sets)],
                               cat.cex = 1.2, cex = 1, margin = 0.05, ind = TRUE))
    make.custom(6,6)
    x = do.call(plot_func,args = argv)
    list(func=plot_func,args=argv)
}


#' Title
#'
#' @param df
#' @param sel_col
#' @param export
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
yplot_venns_bydf = function (df, sel_col, export = F, ...) {
    library(VennDiagram)
    grps = df$Clin_classification %>% unique
    DICT = list(draw.single.venn, draw.pairwise.venn, draw.triple.venn,
                draw.quad.venn, draw.quintuple.venn)
    plot_func = DICT[[length(grps)]]
    sets = list()
    for (i in grps) {
        sets[[i]] = df %>% filter(g == !!i) %>% .[[sel_col]] %>%
            unique
    }
    print(sets %>% names)
    argv = list()
    for (l in 1:length(sets)) {
        for (i in combn(1:length(sets), l, simplify = F)) {
            s = sets[i]
            if (l > 1) {
                r = Reduce(intersect, s) %>% length
                argv[[str_flatten(c("n", i))]] = r
            }
            else {
                r = length(s[[1]])
                argv[[str_flatten(c("area", i))]] = r
            }
        }
    }
    argv$cross.area = argv$n12
    argv = argv %>% ypush(list(category = sets %>% names, fill = c("dodgerblue",
                                                                   "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(sets)],
                               cat.cex = 1.2, cex = 1, margin = 0.05, ind = TRUE))
    print(argv %>% names)
    args = list(...) %>% ypush(list(x = plot_func, args = argv,
                                    export = export, flag = "venn"))
    n = args %>% names
    if (!("suffix" %in% n)) {
        args$suffix = sel_col
    }
    res = ydumpto %>% do.call(args)
    res
}

# ------------- PLOT Enrichment ---------------------
#' Plot the res obj of ydo_GO*
#'
#' @param go.res
#' @param go.tb.filter
#' @param .T_ylab_wrap
#' @param ggsci_col_scheme
#' @param ...
#' @import ggplot2
#' @return
#' @export
#'
#' @examples
yplot_GO_res = function(go.res, go.tb.filter='top10',.T_ylab_wrap = 50,ggsci_col_scheme = NULL,...){
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
                data.tb = res.tb[,-9] %>% group_by(ONTOLOGY) %>% top_n(-top, wt = pvalue) %>% ungroup
            }else{
                warning("go.tb.filter<str> value is wrong, No filtered applied")
                data.tb = res.tb[,-9]
            }
        }else if(is.function(go.tb.filter)){
            print('using custom filter funcion')
            data.tb = go.tb.filter(res.tb[,-9])
        }else{
            data.tb = res.tb[,-9]
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
    res = list(plot.tb = data.tb, gg = gg.go, figsize = figure_size, w = w,
               h = h, height_ratio = height_ratio, maxlen_ylabel = maxlen_ylabel)
    attr(res,'class') = c(class(res),'y.GO.res.plot')
    res

}

