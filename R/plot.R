
#' Plot heatmap using ComplexHeatmap
#'
#' @param mat matrix, usually columns is tsb, if not matrix, will try to coerce using [as.matrix()]
#' @param colData None, or the data.frame with rownames (is tsb) must be accord with colnames(`mat`) and columns of
#'   colData will be draw at the top annotation of Heatmap, each column a line
#' @param annotation_colors define the colors used for draw top annotation, representing colData
#' @param scale default 'auto', use [yscale_rows()] to standardize each rows
#' @param by_group if colData provided
#' @param params_top_annotation plotting params to draw top annotation
#' @param cluster_rows pass to ComplexHeatmap::Heatmap()
#' @param cluster_columns pass to ComplexHeatmap::Heatmap()
#' @param clustering_method_columns pass to ComplexHeatmap::Heatmap()
#' @param clustering_method_rows pass to ComplexHeatmap::Heatmap()
#' @param name pass to ComplexHeatmap::Heatmap()
#' @param row_title pass to ComplexHeatmap::Heatmap()
#' @param show_row_dend pass to ComplexHeatmap::Heatmap()
#' @param show_column_dend pass to ComplexHeatmap::Heatmap()
#' @param show_row_names pass to ComplexHeatmap::Heatmap()
#' @param show_column_names pass to ComplexHeatmap::Heatmap()
#' @param color in format 'Num1~Num2',default '-2~2', will be converted to [circlize::colorRamp()] function then pass to
#'   ComplexHeatmap::Heatmap() as col. if the two Nums span 0, then 0 will be inserted into the middle
#'   ELSE you can provide a custom [circlize::colorRamp()] object
#' @param column_split pass to ComplexHeatmap::Heatmap()
#' @param heatmap_legend_param pass to ComplexHeatmap::Heatmap()
#' @param ... other legal params pass to ComplexHeatmap::Heatmap()
#'
#' @return list with ggplot obj, ComplexHeatmap object and other parameters' values
#' @export
#'
#' @examples
#' \dontrun{
#'
#'     mat = matrix(1:12,nrow=4)
#' }
#'
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
        if (by_group != FALSE){
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
    }else if (is.character(by_group)){
        stopifnot(by_group %in% names(.group))
        .column_split = .group[[by_group]]
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
#' @param diff_expr RNA diff expr, colnames are:
#'   'baseMean','log2FoldChange','abs_log2FC','lfcSE','stat','pvalue','padj'
#' @param value_var ='padj', column name which used for filter
#' @param add_label whether to add the label text of genes
#' @param ...
#' @param logFC_threshold vector of length 2, the log2FoldChange threshold
#' @param p.t pvalue filter threshold which apply to `value_var` column
#'
#' @return list(gg=ggplot,de=data.frame)
#' @export
#'
#' @examples
yplot_volcano_using_de = function (diff_expr,
                                   value_var = "padj",
                                   logFC_threshold = c(-2, 2),
                                   p.t = 0.05,
                                   add_label = FALSE,
                                   ...) {
    if (diff_expr %>% is.data.frame) {
        de = diff_expr
    }
    else {
        de = diff_expr %>% data.frame
    }
    if (length(logFC_threshold)==1){
        logFC_threshold = abs(logFC_threshold)
        logFC_threshold = c(-logFC_threshold,logFC_threshold)
    }
    de$diffexpressed <- "NO"
    de$diffexpressed[de$log2FoldChange > logFC_threshold[2] &
                         de[[value_var]] <
                         p.t] <- "UP"
    de$diffexpressed[de$log2FoldChange < logFC_threshold[1] &
                         de[[value_var]] <
                         p.t] <- "DOWN"
    de$diffexpressed = factor(de$diffexpressed, levels = c('DOWN', 'NO', 'UP'))
    de$y = -log10(de[[value_var]])

    if (add_label) {
        library(ggrepel)
        de$delabel <- NA
        de$gene_symbol = de %>% rownames
        de$delabel[de$diffexpressed != "NO"] <-
            de$gene_symbol[de$diffexpressed !=
                               "NO"]
        g <- ggplot(data = de,
                    aes(
                        x = log2FoldChange,
                        y = y,
                        col = diffexpressed,
                        label = delabel
                    )) + geom_point(alpha = 0.7) +
            geom_text_repel(max.iter = 10, box.padding = 1)
    }
    else {
        g <- ggplot(data = de, aes(x = log2FoldChange, y = y,
                                   col = diffexpressed)) +
            geom_point(alpha = 0.7)
    }
    if (value_var == "padj") {
        old_value_var = value_var
        value_var = "FDR"
    }
    g = g + theme_minimal() +
        ylab(paste0("-log10(", value_var,
                    ")")) +
        scale_color_manual(values = c("blue", "gray",
                                      "red")) +
        geom_vline(xintercept = logFC_threshold,
                   col = "gray",
                   linetype = "longdash") +
        geom_hline(
            yintercept = -log10(p.t),
            col = "gray",
            linetype = "longdash"
        )
    list(gg = g, de = de)
}

#' guess group column name in the data.frame
#'
#' @param df
#'
#' @return colnames vector
#' @export
#'
#' @examples
yinfer_group_col = function(df,group_col=NA,raiseError=TRUE){
    col = colnames(df)
    if (is.na(group_col)) {
        if ("Group" %in% col)
            group_col = "Group"
        else if ("Clin_classification" %in% col)
            group_col = "Clin_classification"
        else if ("g" %in% col)
            group_col = "g"
        else {
            if (raiseError){
                stop("can not determine group_col in df, stop")
            }else
                return(NULL)
        }
    }
    group_col
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
#' @return if .ret==TRUE return res object else Nothing is returned
#' @export
#'
#' @examples
yplot_tsne = function (normd = NULL, colData = NULL, mat = NULL, color_col = NULL,
                       add_label = FALSE, add_polygon=FALSE,
                       perplexity = NA, dims = 2, .ret=FALSE, ...) {
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
    res$data = df

    if (is.null(color_col)){
        group_col = yinfer_group_col(colData,raiseError = FALSE)
        if (is.null(group_col)){
            print("can not determine group_col in colData, stop plotting")
            return(res)
        }
        color_col = sym(group_col)
    }else{
        color_col = sym(color_col)
    }
    gg = df %>% ggplot(aes(x = X1, y = X2, color = !!color_col)) +
        geom_point(size = 3, alpha = 0.7) + ggtitle("t-SNE")
    if (add_polygon){
        gg = gg + stat_ellipse(aes(fill=!!color_col), show.legend = FALSE,
                     geom = "polygon",
                     alpha=0.1,
                     lwd = 0.5)
    }
    res$gg = gg
    if (.ret==TRUE){
        return(res)
    }else{
        print(gg)
    }
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
#' @return list(gg = gg, pca = pca, data = data)
#' @export
#'
#' @examples
yplot_pca = function (normd, colData, color_col = NULL, add_label = FALSE, add_polygon=FALSE,
                      tag_fontsize = 6, max.overlaps = Inf, alpha = 0.6, .ret=FALSE, ...) {
    sig_table = normd %>% t
    mask = (sig_table %>% apply(sd, MARGIN = 2)) > 0
    mat = sig_table[, mask]
    if (ncol(mat) < ncol(sig_table))
        ylog("containing ", sum(!mask), " features with sd==0, removed")
    pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)

    data = cbind(pca$x[colData %>% rownames, 1:2], colData) %>%
        rownames_to_column("Tumor_Sample_Barcode")

    res = list(gg = NULL, pca = pca, data = data)
    if (is.null(color_col)){
        group_col = yinfer_group_col(colData,raiseError = FALSE)
        if (is.null(group_col)){
            print("can not determine group_col in colData, stop plotting")
            return(res)
        }
        color_col = sym(group_col)
    }else
        color_col = sym(color_col)

    gg = data %>% ggplot(aes(PC1, PC2, color = !!color_col, label = Tumor_Sample_Barcode)) +
        geom_point(size = 3, alpha = alpha) + ggtitle("PCA")
    if (add_label == TRUE) {
        gg = gg + ggrepel::geom_label_repel(size = {
            {
                tag_fontsize
            }
        }, alpha = alpha, max.overlaps = max.overlaps)
    }
    if (add_polygon){
        gg = gg + stat_ellipse(aes(fill=!!color_col), show.legend = FALSE,
                               geom = "polygon",
                               alpha=0.1,
                               lwd = 0.5)
    }
    res$gg = gg
    if (.ret==TRUE){
        return(res)
    }else{
        print(gg)
    }

}

#' Plot venn diagram for 2~5 sets
#'
#' @param ... named params or un-named params, each represents a set to be draw, Number of params
#'   min2 max5
#' @param .title NOT implement at current
#' @param plot.argv additional params control the plot appearance
#' @param .ret default FALSE, whether to return plotting data and func, default FALSE
#'
#' @return if .ret==TRUE return list(plot_func=plot_func,args=argv,func_name=plot_func_name) else
#'   nothing is returned
#' @export
#'
#' @examples
yplot_venns = function(...
                       ,.title=NULL
                       ,plot.argv=list(cat.dist = 0.05)
                       ,.ret = FALSE
                       ){
    # if parameters is not named, get the sym name of each parameter
    sets = list(...)
    if (is.null(names(sets))) {
        sets_syms <- match.call(expand.dots = FALSE)$...
        names(sets) <- paste0("arg", seq_along(sets_syms))
        names(sets) <- sapply(sets_syms, deparse)
    }
    # fi

    x = length(sets)
    stopifnot(x <= 5&& x >=1)

    argv = list()
    DICT = list(  draw.single.venn = VennDiagram::draw.single.venn
                  , draw.pairwise.venn = VennDiagram::draw.pairwise.venn
                  , draw.triple.venn = VennDiagram::draw.triple.venn
                  , draw.quad.venn = VennDiagram::draw.quad.venn
                  , draw.quintuple.venn = VennDiagram::draw.quintuple.venn)

    plot_func = DICT[[x]]
    plot_func_name = names(DICT)[[x]]

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
    argv = argv %>% ypush(list(category = sets %>% names
                               , fill = c("dodgerblue",
                                          "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(sets)],
                               cat.cex = 1.2, cex = 1, margin = 0.05, ind = TRUE
    )%>%c(plot.argv))
    make.custom(6,6)
    x = do.call(plot_func,args = argv)

    if (.ret == TRUE){
        return(list(plot_func=plot_func,args=argv,func_name=plot_func_name))
    }
}


#' Title
#'
#' @param id_col
#' @param group_col
#' @param .ret if TRUE return the res object
#' @param df
#'
#' @return if .ret==TRUE return the res object else nothing is returned
#' @export
#'
#' @examples
yplot_venns_bydf = function (df,id_col, group_col='Clin_classification', .ret=FALSE, ...) {
    # my_list <- lapply(unique(df[[group_col]]), function(x) {
    #     df[[id_col]][df[[group_col]] == x]
    # })
    # names(my_list) <- unique(df[[group_col]])

    sym_group_col = sym(group_col)
    sym_id_col = sym(id_col)
    my_list = df %>%
        group_by(!!sym_group_col) %>%
        summarize(values = list(!!sym_id_col)) %>%
        ungroup() %>%
        {setNames(.$values, .[[group_col]])}

    res = yplot_venns %>% do.call(my_list)
    if (.ret == TRUE){
        return(res)
    }
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
#' @return an y.GO.res.plot object, list(plot.tb = data.tb, gg = gg.go, figsize = figure_size, w = w,
#  h = h, height_ratio = height_ratio, maxlen_ylabel = maxlen_ylabel)
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


#' add plots to y.GSEA.resx objects
#'
#' @param res an 'y.GSEA.resx' object, generated by ydo_GSEA()
#' @param .lineheight adjust the line-height of the left y-axis labels
#' @param .T_ylab_wrap the width of left y-axis labels
#' @param .lineheight yaxis label line heights
#' @param .legends.color.title the title in legend "direction" which is the legend title of the
#'   .legends.color.labels
#' @param .legends.color.labels the labels of the color legend, length 2 character vector, apply to
#'   values c(-1,1)
#' @param ...
#' @param .TOP plot top N features, default 30
#' @param p.adjust.t,p.t if you want to re-filter the enrichement results with different thresholds
#'   of FDR and pvalue thresholds, set the parameters to Non-NULL values
#' @param figureTitle
#'
#' @return an y.GSEA.plots object
#' @export
#'
#' @examples
yplot_GSEA_res = function(res,
                          p.adjust.t = NULL,
                          p.t = NULL,
                          .TOP = 30,
                          figureTitle = 'GSEA of DEGs',
                          .lineheight = 0.75,
                          .T_ylab_wrap = 40,
                          .legends.color.title = 'Expression',
                          .legends.color.labels = c('Down-regulated','Up-regulated'),
                          ...){
    stopifnot('y.GSEA.resx' %in% class(res))
    npg = ggsci::pal_npg()(10)[2:1]

    for (key in names(res)) {
        res.gsea = res[[key]]
        if (!is.null(p.adjust.t) || !is.null(p.t)){
            p.adjust.t = ifelse(test = is.null(p.adjust.t),yes=1,no = p.adjust.t)
            p.t = ifelse(test = is.null(p.t),yes=1,no = p.t)
            print(stringr::str_flatten(c('Re-filtering using p.t=', p.t,'p.adjust.t=',p.adjust.t)))

            res.tb.gsea = res[[key]]$res.gsea@result %>% dplyr::filter(p.adjust < p.adjust.t, pvalue < p.t)
        }
        else{
            res.tb.gsea = res[[key]]$res.tb.gsea
        }
        n.res = res.tb.gsea %>% nrow

        if (n.res > 0){
            res.tb.gsea %<>% mutate(
                Description = forcats::fct_reorder(Description,
                                                   p.adjust, .desc = TRUE),
                direction = factor(
                    sign(NES),
                    levels = c(-1, 1),
                    labels = .legends.color.labels
                ),
                EnrichedGeneNum = core_enrichment %>%
                    str_split("/") %>% sapply(length),
                EnrichedPer = round(EnrichedGeneNum / setSize *100)
            ) %>% dplyr::select(!core_enrichment)

            tab = res.tb.gsea %>% dplyr::select(Description, EnrichedGeneNum, EnrichedPer)
            maxlen_ylabel = res.tb.gsea$ID %>% sapply(str_length) %>% max

            if (maxlen_ylabel > .T_ylab_wrap) {
                ylabel_size = 10
                height_ratio = 1 + 0.15 * (maxlen_ylabel %/% .T_ylab_wrap)
                maxlen_ylabel = .T_ylab_wrap
            }
            else {
                ylabel_size = 14
                height_ratio = 1
            }
            w = ceiling((maxlen_ylabel / 10 + 7) * 2) / 2
            h = ((round(min(
                30, n.res
            ) / 4 + 3) * height_ratio) * 2) / 2

            gg.gsea = res.tb.gsea %>% head(30) %>%
                ggplot(aes(
                    x = -log10(pvalue),
                    y = Description,
                    color = direction
                )) +
                geom_point(aes(size = EnrichedPer),
                           alpha = 0.75) +
                scale_shape_manual(values = c(15, 16, 17, 18)) +
                scale_color_manual(breaks = .legends.color.labels,
                                   values = npg) +
                ylab(NULL) +
                ggtitle(figureTitle) +
                theme_bw(base_size = 16) +
                theme(axis.text.y.left =
                          element_text(size = ylabel_size,
                                       lineheight = .lineheight)) +
                scale_y_discrete(
                    position = "left",
                    labels = function(x)
                        str_wrap(str_replace_all(x, fixed("_"), " "), width = maxlen_ylabel)
                ) +
                guides(
                    y.sec = ggh4x::guide_axis_manual(
                        breaks = tab$Description,
                        labels = paste0(tab$EnrichedGeneNum, "(", tab$EnrichedPer,
                                        "%)")
                    ),
                    size = guide_legend(title = "Enriched%"),
                    color = guide_legend(title = .legends.color.title)
                )

            make.custom(w, h)
            figure_size = c(w, h)
            plotted = TRUE
        }else{
            w = 8
            h = 6
            figure_size = c(8,6)
            height_ratio = 1
            maxlen_ylabel = 10
            plotted = FALSE
        }

        res[[key]] = ypush(res[[key]],list(
            gg = gg.gsea,
            data.tb = res.tb.gsea,
            w = w,
            h = h,
            figsize = figure_size,
            height_ratio = height_ratio,
            maxlen_ylabel = maxlen_ylabel,
            plotted = plotted
        ))
        if (! 'y.GSEA.plot' %in% res[[key]]){
            attr(res[[key]], 'class') = c('y.GSEA.plot',class(res[[key]]))
        }
    }
    if (! 'y.GSEA.plots' %in% res){
        attr(res, 'class') = c('y.GSEA.plots',class(res))
    }
    res
}



