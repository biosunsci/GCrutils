#  ----------------- Private Global VAR ---------------
SSGSEA_BASE_REF_DIR = './extdata/'
COMMON_GMT_TAGS = c('hm', 'kegg', 'im', 'im_jnj')

# ---------------- Enrich Func -----------------------------
#' GSEA analysis using different expression results
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
ydo_gsea = function (diff = NULL,
                     gene_list = NULL,
                     set = "gsea",
                     p.adjust.t = 0.05,
                     p.t = NA,
                     minGSSize = 10,
                     maxGSSize = 500,
                     .lineheight = 0.75,
                     pvalueCutoff = 0.2,
                     ...)
{
    "set in c('all','all3','gsea','hm','kegg','immune','immune.jnj') "
    npg = pal_npg()(10)[2:1]
    .T_ylab_wrap = 40
    refs = list(
        hm = c(figureTitle = "DiffExpr Hallmark enricment",
               Name = "Hallmark"),
        kegg = c(figureTitle = "DiffExpr KEGG enrichment",
                 Name = "KEGG"),
        immune = c(figureTitle = "DiffExpr Immune Infiltration enrichment",
                   Name = "Immune"),
        immune.jnj = c(figureTitle = "DiffExpr Immune Infiltration2 enrichment",
                       Name = "Immune2")
    )
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
        res.gsea = GSEA(
            geneList = gene_list,
            TERM2GENE = gmt,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            ...
        )
        res.tb.gsea = res.gsea@result %>% dplyr::filter(p.adjust <
                                                            p.adjust.t)
        n.res = res.tb.gsea %>% nrow
        print(paste("after filter,", n.res, "results left"))
        if (n.res > 0) {
            res.tb.gsea %<>% mutate(
                Description = fct_reorder(Description,
                                          p.adjust, .desc = TRUE),
                EnrichedGeneNum = core_enrichment %>%
                    str_split("/") %>% sapply(length),
                direction = factor(
                    sign(NES),
                    levels = c(-1, 1),
                    labels = c("DEC Genes", "INS Genes")
                ),
                EnrichedPer = round(EnrichedGeneNum /
                                        setSize *
                                        100)
            ) %>% dplyr::select(!core_enrichment)
            tab = res.tb.gsea %>% dplyr::select(Description,
                                                EnrichedGeneNum, EnrichedPer)
            maxlen_ylabel = res.tb.gsea$ID %>% sapply(str_length) %>%
                max
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
            ) / 4 + 3) * height_ratio) *
                2) / 2
            gg.gsea = res.tb.gsea %>% head(30) %>% ggplot(aes(
                x = -log10(pvalue),
                y = Description,
                color = direction
            )) + geom_point(aes(size = EnrichedPer),
                            alpha = 0.75) + scale_shape_manual(values = c(15,
                                                                          16, 17, 18)) + scale_color_manual(breaks = c("DEC Genes",
                                                                                                                       "INS Genes"),
                                                                                                            values = npg) + ylab(NULL) + ggtitle(figureTitle) +
                theme_bw(base_size = 16) + theme(axis.text.y.left = element_text(size = ylabel_size,
                                                                                 lineheight = .lineheight)) + scale_y_discrete(
                                                                                     position = "left",
                                                                                     labels = function(x)
                                                                                         str_wrap(str_replace_all(x,
                                                                                                                  fixed("_"), " "), width = maxlen_ylabel)
                                                                                 ) +
                guides(
                    y.sec = ggh4x::guide_axis_manual(
                        breaks = tab$Description,
                        labels = paste0(tab$EnrichedGeneNum, "(", tab$EnrichedPer,
                                        "%)")
                    ),
                    size = guide_legend(title = "Enriched%"),
                    color = guide_legend(title = "Direction")
                )
            make.custom(w, h)
            figure_size = c(w, h)
            res[[key]] = list(
                gsea.set = Name,
                gene_list = gene_list,
                res.gsea = res.gsea,
                tb.gsea = res.tb.gsea,
                gg = gg.gsea,
                figsize = figure_size,
                w = w,
                h = h,
                height_ratio = height_ratio,
                maxlen_ylabel = maxlen_ylabel
            )
        }
        else {
            print(
                paste0(
                    "[!] Handling ",
                    key,
                    ": too less res got after filtering, consider to increase @p.adjust.t"
                )
            )
            res[[key]] = list(
                gsea.set = Name,
                gene_list = gene_list,
                res.gsea = res.gsea,
                res.tb.gsea = res.tb.gsea,
                gg = NA,
                figsize = NA
            )
        }
    }
    message("GSEA analysis done, now you can use yplot_")
    res
}


#' GSEA analysis using different expression results
#'
#' @param diff a data.frame which is a output df of DESeq2 different expr analysis output, with
#'   columns: symbol, log2FC, pvalue, padj, stat
#' @param gene_list list, gene symbols as names and `stat` or `p-value` as values, decreasing
#'   arranged. will overrite the diff as input
#' @param set set in c('all','all3','gsea','hm','kegg','im','im_jnj') 'all' - all the gmt
#'   sets in the  refs, currently: Hallmark, KEGG, Immune, Immune2(im_jnj) 'all3' - hm, kegg,
#'   immune do the 3 gmt sets 'gsea' - hm, kegg do the 2 gmt sets 'hm','kegg','im','im_jnj'
#'   use the corresponding single gmt set
#' @param p.t pvalue threshold, default 0.01
#' @param p.adjust.t padj threshold, default 0.05
#' @param log2fc.abs.t log2FC threshold, default 1
#' @param ...
#'
#' @return a list obj of class 'y.GO.res', used for input for yplot_GO_res
#' @export
#'
#' @examples
ydo_GSEA = function (diff = NULL,
                     gene_list = NULL,
                     set = "gsea",
                     p.adjust.t = 0.05,
                     p.t = NA,
                     minGSSize = 10,
                     maxGSSize = 500,
                     .lineheight = 0.75,
                     pvalueCutoff = 0.2,
                     ...) {
    npg = pal_npg()(10)[2:1]
    .T_ylab_wrap = 40
    refs = list(
        hm = c(figureTitle = "Hallmark enricment",
               Name = "Hallmark"),
        kegg = c(figureTitle = "KEGG enrichment",
                 Name = "KEGG"),
        im = c(figureTitle = "Immune Infiltration enrichment",
               Name = "Immune"),
        im_jnj = c(figureTitle = "Immune Infiltration enrichment",
                   Name = "Immune2")
    )
    if (set == "all") {
        run = refs
    }
    else if (set == "gsea") {
        run = refs[c("hm", "kegg")]
    }
    else if (set == "all3") {
        run = refs[c("hm", "kegg", "im")]
    }
    else if (set %in% c("hm", "kegg", "im", "im_jnj")) {
        run = refs[set]
    }
    else {
        stop("set must in in c('all','all3','gsea','hm','kegg','im','im_jnj')")
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
        res.gsea = GSEA(
            geneList = gene_list,
            TERM2GENE = gmt,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalueCutoff,
            ...
        )
        res.tb.gsea = res.gsea@result %>% dplyr::filter(p.adjust <
                                                            p.adjust.t)
        n.res = res.tb.gsea %>% nrow
        print(paste("after filter,", n.res, "results left"))
        if (n.res > 0) {
            res.tb.gsea %<>% mutate(
                Description = fct_reorder(Description,
                                          p.adjust, .desc = TRUE),
                EnrichedGeneNum = core_enrichment %>%
                    str_split("/") %>% sapply(length),
                direction = factor(
                    sign(NES),
                    levels = c(-1, 1),
                    labels = c("DEC Genes", "INS Genes")
                ),
                EnrichedPer = round(EnrichedGeneNum /
                                        setSize *
                                        100)
            ) %>% dplyr::select(!core_enrichment)
            tab = res.tb.gsea %>% dplyr::select(Description,
                                                EnrichedGeneNum, EnrichedPer)
            maxlen_ylabel = res.tb.gsea$ID %>% sapply(str_length) %>%
                max
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
            ) / 4 + 3) * height_ratio) *
                2) / 2
            gg.gsea = res.tb.gsea %>% head(30) %>%
                ggplot(aes(
                    x = -log10(pvalue),
                    y = Description,
                    color = direction
                )) +
                geom_point(aes(size = EnrichedPer),
                           alpha = 0.75) +
                scale_shape_manual(values = c(15,
                                              16, 17, 18)) +
                scale_color_manual(breaks = c("DEC Genes",
                                              "INS Genes"),
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
                    color = guide_legend(title = "Direction")
                )
            make.custom(w, h)
            figure_size = c(w, h)
            res[[key]] = list(
                gsea.set = Name,
                gene_list = gene_list,
                res.gsea = res.gsea,
                tb.gsea = res.tb.gsea,
                gg = gg.gsea,
                figsize = figure_size,
                w = w,
                h = h,
                height_ratio = height_ratio,
                maxlen_ylabel = maxlen_ylabel
            )
        }
        else {
            print(
                paste0(
                    "[!] Handling ",
                    key,
                    ": too less res got after filtering, consider to increase @p.adjust.t"
                )
            )
            res[[key]] = list(
                gsea.set = Name,
                gene_list = gene_list,
                res.gsea = res.gsea,
                res.tb.gsea = res.tb.gsea,
                gg = NA,
                figsize = NA
            )
        }
    }
    message("GSEA analysis done, now you can use yplot_")
    res
}




#'
#' Do GO analysis using different expression results
#'
#' @details  如果genes is.null, use p.t or p.adjust.t && log2fc.abs.t (as FILTERS) to filter diff to
#'   get genes ' diff必须包含全部基因(as
#'   ALL)而不是基因子集,GO分析才准确(https://mp.weixin.qq.com/s/zxiMXqRriGHcGAlCNk0e7A) '
#'   若genes(格式为symbols)不为空,使用genes/ALL进行GO分析,若genes为空, 使用FILTERS从diff中筛选出genes"
#'
#' @param diff a data.frame with columns: symbol, log2FC, pvalue, padj, stat
#' @param ont Ontology, in c('BP', 'MF', 'CC', 'ALL')
#' @param by.expr.direction =TRUE, if TRUE genes in diff is group_by change direction before doing
#'   GO analysis
#' @param p.adjust.t padj threshold, default 0.05
#' @param ...
#' @param go.tb.filter
#' @param pvalueCutoff
#' @param t.p.val
#' @param t.log.fc
#' @param t.p.adj
#'
#' @return a list obj of class 'y.GO.res', used for input for yplot_GO_res
#' @export
#'
#' @examples
ydo_GO_bydiff = function(diff,
                         ont = "ALL",
                         by.expr.direction = FALSE,
                         go.tb.filter = 'top10',
                         pvalueCutoff = 0.2,
                         p.adjust.t = 0.05,
                         t.p.val = NA,
                         t.log.fc = 1,
                         t.p.adj = 0.05,
                         ...) {
    ont = toupper(ont)
    sig_genes = yutils_filter_stat_res(diff,
                                       t.p.val = t.p.val,
                                       t.log.fc = t.log.fc,
                                       t.p.adj = t.p.adj)

    if (!yhas_rownames(diff)) {
        diff = diff %>% rownames_to_column("symbol")
    }
    all_genes = rownames(diff)

    Name = "GO_enrich"
    tryCatch(
        expr = {
            GENE.REF = get("GENE.REF", envir = .GlobalEnv)
            print("Got Global var GENE.REF")
        },
        error = function(err) {
            assign(
                "GENE.REF",
                value = bitr(
                    all_genes,
                    fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db
                ) %>%
                    arrange(SYMBOL) %>% distinct(SYMBOL, .keep_all = TRUE) %>%
                    column_to_rownames("SYMBOL"),
                envir = .GlobalEnv
            )
            print("Global var GENE.REF assigned")
        }
    )

    if (by.expr.direction == TRUE) {
        res.tb = list()
        res.ego = list()
        for (direction in c('Up', 'Down')) {
            genes = sig_genes %>% filter(FC_Ins == (direction == "Up")) %>% rownames
            ego = enrichGO(
                gene = GENE.REF[genes, "ENTREZID"],
                universe = GENE.REF$ENTREZID,
                ont = ont,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                pAdjustMethod = "BH",
                pvalueCutoff = pvalueCutoff,
                qvalueCutoff = p.adjust.t,
                readable = TRUE,
                ...
            )
            res.tb[[direction]] = ego@result %>% mutate(Change = factor(!!direction))
            res.ego[[direction]] = ego
        }
        res.tb = rbind(res.tb[[1]], res.tb[[2]]) %>%
            mutate(
                Description = fct_reorder(Description, p.adjust, .desc = TRUE),
                ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
                EnrichedGeneNum = Count,
                setSize = BgRatio %>%
                    str_split("/", simplify = TRUE) %>% .[, 1] %>% as.integer,
                EnrichedPer = round(EnrichedGeneNum / setSize * 100)
            )
    } else{
        genes = sig_genes %>% rownames
        res.ego = enrichGO(
            gene = GENE.REF[genes, "ENTREZID"],
            universe = GENE.REF$ENTREZID,
            ont = ont,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            pAdjustMethod = "BH",
            pvalueCutoff = p.t,
            qvalueCutoff = p.adjust.t,
            readable = TRUE,
            ...
        )
        res.tb = res.ego@result %>%
            mutate(
                Description = fct_reorder(Description, p.adjust, .desc = TRUE),
                ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
                EnrichedGeneNum = Count,
                setSize = BgRatio %>%
                    str_split("/", simplify = TRUE) %>% .[, 1] %>% as.integer,
                EnrichedPer = round(EnrichedGeneNum / setSize * 100),
                Change = 'DiffExpr'
            ) %>% dplyr::select(!geneID)
    }
    attr(res.tb, 'class') = c(class(res.tb), 'y.GO.res.dtable')
    print('GO analysis Done')
    res = list(
        gsea.set = Name,
        gene_list = sig_genes,
        res.tb = res.tb,
        res.go = res.ego,
        args = list(
            ont = ont,
            by.expr.direction = by.expr.direction,
            p.t = p.t,
            p.adjust.t = p.adjust.t,
            log2fc.abs.t = t.log.fc
        )
    )
    attr(res, 'class') = c(class(res), 'y.GO.res')
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
ydo_GO = function(sig_genes,
                  all_genes = NULL,
                  ont = 'ALL',
                  p.t = 0.01,
                  p.adjust.t = 0.05,
                  ggsci_col_scheme = NULL,
                  ...) {
    .T_ylab_wrap = 50
    ont = toupper(ont)
    if (is.null(ggsci_col_scheme)) {
        npg = pal_npg()(10)[c(3, 6, 9)]
    } else{
        npg = ggsci_col_scheme
    }
    if (is.null(all_genes)) {
        all_genes = readRDS('./extdata/all_genes_symbol.RDS')
        print(
            'all_genes is NULL, Read 19039 gene symbols from extdata/all_genes_symbol.RDS -> all_genes'
        )
    }
    if (!is.character(sig_genes)) {
        print('Try to Convert sig_genes into charactor')
        sig_genes = as.character(sig_genes)
    }

    figureTitle = paste("DiffExpr GO enrichment", ont)
    Name = 'GO_enrich'

    # 生成symbol-> entrezID的map：GENE.REF
    tryCatch(
        expr = {
            GENE.REF = get("GENE.REF", envir = .GlobalEnv)
            stopifnot(is.data.frame(GENE.REF))
            stopifnot(nrow(GENE.REF) > 0)
            print('Got Global var GENE.REF')
        },
        error = function(err) {
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
                    column_to_rownames('SYMBOL'),
                envir = .GlobalEnv
            )
            GENE.REF = get("GENE.REF", envir = .GlobalEnv)
            print('Global var GENE.REF assigned')
        }
    )

    rm.genes = setdiff(sig_genes, rownames(GENE.REF))

    if (length(rm.genes) > 0) {
        sig_genes = intersect(sig_genes, rownames(GENE.REF))
        print(paste0(
            'remove ',
            length(rm.genes),
            ' genes symbol not in GENE.REFß:'
        ))
        print(rm.genes)
    }
    # sig_genes/diff为GO分析的主体数据, 使用ENTREZID做分析, 使用全局变量GENE.REF作为map转换symbols->ENTREZID
    ego = enrichGO(
        gene       = GENE.REF[sig_genes, 'ENTREZID'],
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
    if (nrow(ego@result) == 0) {
        print('GO analysis result in 0 terms enriched')
        res.tb = data.frame()
    } else{
        res.tb = ego@result %>% mutate(
            Description = fct_reorder(Description, p.adjust, .desc = TRUE),
            ONTOLOGY = factor(ONTOLOGY, levels = c('BP', 'MF', 'CC')),
            EnrichedGeneNum = GeneRatio %>% str_split('/', simplify = TRUE) %>% `[`(, 1) %>% as.integer,
            setSize = BgRatio %>% str_split('/', simplify = TRUE) %>% `[`(, 1) %>% as.integer,
            EnrichedPer = round(EnrichedGeneNum / setSize * 100)
        )
    }

    res = list(
        gsea.set = Name,
        gene_list = sig_genes,
        res.tb = res.tb,
        res.go = ego,
        args = list(
            ont = ont,
            by.expr.direction = FALSE,
            p.t = p.t,
            p.adjust.t = p.adjust.t,
            log2fc.abs.t = NA
        )
    )
    attr(res, 'class') = c(class(res), 'y.GO.res')
    message('GO analysis done, now you can use yplot_GO_res to get the figure')

    return(res)
}

#' do_ssGSEA analysis
#'
#' @description depends on yload_list_gmt
#'  use the GSVA::gsea method which need 2 critical parameters
#'  a <expr_matrix> is the numeric matrix of expression, normed count of DSeq2 or tpm
#'  a <ref_gmt_list> which is the returns of yload_list_gmt, see details of yload_list_gmt documentation
#'
#' @param expr_matrix a numeric matrix witch rows are genes and columns are patients
#' @param todo list, combinations of elements in list('hm','kegg','im','im_jnj')
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
ydo_ssGSEA = function(expr_matrix,
                      todo = list('hm', 'kegg', 'im', 'im_jnj'),
                      ref.gmt = NULL,
                      ref.gmt.name = NULL,
                      mx.diff = FALSE,
                      verbose = FALSE,
                      method = 'ssgsea',
                      kcdf = 'Gaussian',
                      parallel.sz = 1,
                      ...) {
    res = list()
    if (!is.null(ref.gmt)) {
        if (is.character(ref.gmt)) {
            ref_df = yload_list_gmt(gmt_txt_file_path = ref.gmt)
        }
        if (is.null(ref.gmt.name)) {
            ref.gmt.name = 1
        }
        es = GSVA::gsva(
            expr_matrix,
            ref_df,
            mx.diff = mx.diff,
            verbose = verbose,
            method = method,
            kcdf = kcdf,
            parallel.sz = 1
        )
        res[[ref.gmt.name]] = es
    } else{
        for (i in 1:length(todo)) {
            k = todo[[i]]
            if (!is.null(ref.gmt)) {
                ref_df = ref.gmt
            } else{
                ref_df = yload_list_gmt(k)
            }
            es = GSVA::gsva(
                expr_matrix,
                ref_df,
                mx.diff = mx.diff,
                verbose = verbose,
                method = method,
                kcdf = kcdf,
                parallel.sz = 1,
                ...
            )
            # es %>% ydumpto(flag='ssGSEA_result_',export=k,outputdir=OUTPUTROOT,row.names=T)
            res[[k]] = es
        }
    }
    res
}


# ------------- helpers --------------

#' load symbol.gmt files using [clusterProfile::read.gmt()] into data.frame format obj
#'
#' @description see docs of yload_list_gmt, yload_gmt
#'
#' @details The source of the gmt file is:
#' hm = 'h.all.v7.5.symbols.gmt',
#' kegg = 'c2.cp.kegg.v7.5.symbols.gmt',
#' im = 'immune.cell.rep.symbols.gmt',
#' im_jnj = 'immune.jnj.immunity.symbols.gmt'
#'
#' ref = list( hm='hm.df', kegg='kegg.df', im='im.df',im_jnj='im_jnj.df')
#'
#'
#' @param key one of names(`ref`), default one of c('hm','kegg','im','im_jnj'),
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
#'\dontrun{
#' gmt = yload_gmt('hm')
#' }
#'
yload_gmt = function(key = NULL,
                     base_ref_dir = NULL) {
    if (key %in% COMMON_GMT_TAGS) {
        f = file.path('extdata/', paste0(key, '.symbol.df.RDS'))
        f = system.file(f, package = 'GCrutils')
        gmt = readRDS(f)
    } else{
        print('else')
        if (is.null(base_ref_dir)) {
            base_ref_dir = SSGSEA_BASE_REF_DIR
        }
        path = file.path(base_ref_dir, key)
        stopifnot(file.exists(path))
        gmt = clusterProfiler::read.gmt(path)
    }
    gmt
}

#' read gmt file as a list
#'
#' read Gene Matrix Transposed (gmt) file and output a list with the the first
#' column as the names of items in the list. see
#' \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{Gene Matrix Transposed file format}
#' for more details.
#' @seealso yload_gmt
#'
#' @param gmtfile a gmt file. Examples are from MSigDB Collections.
#' A list of gene set could be find in the vignette of cogena
#' @return a gmt list
#' @export
#'
#' @seealso gmtlist2file
#' @examples
#' \dontrun{
#' anno <- "c2.cp.kegg.v7.01.symbols.gmt.xz"
#' gmtfile <- system.file("extdata", anno, package="cogena")
#' gl <- gmt2list(gmtfile)
#' }
#'
yload_gmtfile2list <- function(gmtfile) {
    if (!file.exists(gmtfile)) {
        stop("There is no such a gmt file!")
    }
    if (tools::file_ext(gmtfile) == "xz") {
        gmtfile <- xzfile(gmtfile)
        x <- scan(gmtfile,
                  what = "",
                  sep = "\n",
                  quiet = TRUE)
        close(gmtfile)
    } else if (tools::file_ext(gmtfile) == "gmt") {
        x <- scan(gmtfile,
                  what = "",
                  sep = "\n",
                  quiet = TRUE)
    } else {
        stop ("Only gmt and gmt.xz are accepted for gmt2list")
    }

    y <- strsplit(x, "\t")
    names(y) <- sapply(y, `[[`, 1)

    annoList <- lapply(y, `[`, c(-1, -2))
    return(annoList)
}

#' load gmt files using read.table into list format obj
#'
#' @description see docs of yload_list_gmt, yload_gmt
#' @seealso yload_gmt
#' @param key one of names(`ref`), default one of c('hm','kegg','im','im_jnj')
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
#' \dontrun{
#' gmt.list = yload_list_gmt('hm')
#' }
yload_list_gmt = function(key = NULL,
                          base_ref_dir = NULL) {
    if (key %in% COMMON_GMT_TAGS) {
        f = file.path('extdata/', paste0(key, '.symbol.list.RDS'))
        f = system.file(f, package = 'GCrutils')
        gmt_list = readRDS(f)

    } else{
        if (is.null(base_ref_dir)){
            base_ref_dir = SSGSEA_BASE_REF_DIR
        }
        path = file.path(base_ref_dir, key)
        stopifnot(file.exists(path))
        # tmp = read.table(
        #     path,
        #     quote = '"',
        #     header = T,
        #     sep = ',',
        #     row.names = 1
        # )
        # gmt_list = split(tmp[, 1], tmp[, 2])
        gmt_list = yload_gmtfile2list(path)

    }
    gmt_list
}
