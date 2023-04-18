yget_subset_maf = function(laml=NULL,maf_=NULL,gl_=NULL,export=F,include_total=F,retmode='maf.list',...){
    # retmode in ['maf']
    if (retmode=='maf.list'){
        laml = yload_laml_maf(maf_,gl_,laml)
        glx = laml %>% getClinicalData
        sub_mafs = list()
        grps = glx$Clin_classification %>% unique()
        
        if (length(grps)==1){
            warning('Only 1 subgroup detected, set include_total <- FALSE')
            include_total = FALSE
        }
        if (include_total==TRUE){
            # all together    
            sub_mafs[[length(sub_mafs)+1]] = laml
            names(sub_mafs)[[sub_mafs %>% length]] = 'ALL'
        }
        # sub groups
        for (g in grps){
            submaf = subsetMaf(laml,clinQuery = paste0("Clin_classification == '",g,"'"))   
            sub_mafs[[length(sub_mafs)+1]] = submaf
            names(sub_mafs)[[sub_mafs %>% length]] = g
        }    
        return(sub_mafs)
    }
    
}


yload_laml_maf = function(maf_=NULL,gl_=NULL,laml=NULL,only_12_cols=TRUE,...){
    if (is.null(laml) && !is.null(maf_) && !is.null(gl_)) {
        maf_ = yload_dfx(maf_,...)
        if (only_12_cols==TRUE) maf_ = maf_[,1:12]
        gl_ = yload_dfx(gl_,...)
        laml = read.maf(maf_,clinicalData = gl_)
    }else if(!is.null(laml)){
        #pass
    }else stop('laml,maf_,gl_ all NULL')
    laml
}
#-----------------

# for repel labeling
# load showtext if you want to plot Chinese charactors 
# library(showtext)
# showtext_auto()

# yplot_pca = function(cnt_file_=NULL,glx_=NULL
#                       ,normd =NULL
#                       ,colData=NULL
#                       ,frm=INPUTROOT
#                       ,export=FALSE
#                       ,outputdir=OUTPUTROOT
#                       ,export_matrix = FALSE
#                       ,flag = 'PCA'
#                       ,add_label = F
#                       ,tag_fontsize = 6
#                       ,max.overlaps=Inf
#                       ,alpha =0.6
#                       ,...){
#     # 
#     if (!is.null(cnt_file_) && !(is.null(glx_))){
#         r = ynormalize_count_bydeseq2(cnt_file_,glx_,frm=frm,diff=F)
#         normd = r$normd
#         colData = r$normd_count
#     }else if (!is.null(normd) && !(is.null(colData))){
#         # pass
#     }else{
#         stop('c(cnt_file_,glx_) && c(normd,colData) both null')
#     }
#     sig_table = normd %>% t
#     mask = (sig_table %>% apply(sd,MARGIN = 2)) > 0
#     mat = sig_table[,mask]
#     if (ncol(mat) < ncol(sig_table)) ylog('containing ',sum(!mask),' features with sd==0, removed')
#     pca <- prcomp(mat, center = TRUE, scale. = TRUE)
#     data = cbind(pca$x[colData%>%rownames,1:2],
#                colData) %>% 
#         rownames_to_column("Tumor_Sample_Barcode") 

#     # mutate(a=(rowname %>% str_sub(1,6))) %>%
#     gg =  data %>% 
#         ggplot(aes(PC1,PC2,color=Clin_classification,label=Tumor_Sample_Barcode)) +
#         geom_point()
#     if (add_label==TRUE){
#         gg = gg + geom_label_repel(size={{tag_fontsize}},alpha=alpha,max.overlaps=max.overlaps)
#     }        
#     argv = yget_args(.f = ydumpto) %>% ypush(list(x=gg,...),.expand = T) 
#     # argv =  environment() %>% as.list() %>% ygetlast(formals(ydumpto) %>% names) %>% ypush(list(x=gg))
#     ydumpto %>% do.call(args=argv)
#     print(argv %>% names)
#     list(gg=gg,pca=pca,data=data,argv=argv)
# }

yplot_pca = function(normd,colData,color_col = 'Group',add_label=FALSE,tag_fontsize = 6, 
    max.overlaps = Inf, alpha = 0.6,...){
    sig_table = normd %>% t
    mask = (sig_table %>% apply(sd, MARGIN = 2)) > 0
    mat = sig_table[, mask]
    if (ncol(mat) < ncol(sig_table)) 
        ylog("containing ", sum(!mask), " features with sd==0, removed")
    pca <- prcomp(mat, center = TRUE, scale. = TRUE)
    color_col = sym(color_col)
    data = cbind(pca$x[colData %>% rownames, 1:2], colData) %>% 
        rownames_to_column("Tumor_Sample_Barcode")
    gg = data %>% ggplot(aes(PC1, PC2, color = !!color_col, 
        label = Tumor_Sample_Barcode)) + 
        geom_point(size=3,alpha=0.7) + 
        ggtitle('PCA')
    if (add_label == TRUE) {
        library(ggrepel)
        gg = gg + geom_label_repel(size = {{ tag_fontsize }}, alpha = alpha, max.overlaps = max.overlaps)
    }
    list(gg = gg, pca = pca, data = data)
}
#-----------------

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


yget_grps = function(glx,order=NULL){
    cols = colnames(glx)
    if ('Group' %in% cols){
        g = 'Group'
    }else if ('Clin_classification' %in% cols){
        g = 'Clin_classification'
    }else if ('g' %in% cols){
        g = 'g'
    }else{
        stop("None of c('Group','Clin_classification','g') not in colnames(glx) ")
    }
    grps = colData[[g]] %>% unique
    if (!is.null(order)){
        grps = factor(grps,levels=order)
    }
    list(col=g,grps=grps)
}


# library('pheatmap')
# library('RColorBrewer')
#' force_exe if nrow(mat) > 2000 or ncol(mat) > 1000, stop plotting unless force_exe is TRUE
#' clustering_method same as in hclust 
#' the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#' "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
yplot_heatmap = function(mat
                        ,colData=NULL
                        ,bygroup=FALSE
                        ,anno.palette = 'npg'
                        ,zscore.rize=TRUE
                        ,force_exe=FALSE
                     
                        ,clustering_method = 'ward.D2' #
                        ,cluster_rows = TRUE
                        ,cluster_cols = TRUE
                         
                        ,color = colorRampPalette(c('blue', "white",'red'))(51)
                        ,fontsize = 14
                        ,cellwidth = 3
                        ,cellheight = 0.2
                        ,border_color = 'white'
                        ,treeheight_row = 0
                        ,annotation_colors = NA
                        ,show_rownames = F
                        ,show_colnames = F
                        ,breaks = seq(-2,2,length.out=51)
                         # ,legend_breaks = c(-2,-1,0,1,2)
                        ,fontsize_col = 8
                          ,...){
    # generate annotation_colors for each group
    if (!is.null(colData)){
        colorx = list()
        anno_color = list()
        if (colData %>% is.data.frame){
            anno_color = ygenerate_annotation_colors(colData
                                                     ,annotation_colors = annotation_colors
                                                     ,palette=anno.palette
                                                    )
        }else{
            stop('unsupported colData type, must data.frame')
            # ! may cause exception
        }
    }else{
        anno_color = NULL
        bygroup = FALSE
        print('colData is NULL, set bygroup=FALSE')
    }
    
    if ('data.frame' %in% class(mat)){
        mat = as.matrix(mat)
    }
    if('matrix' %in% class(mat)){
        stopifnot(is.numeric(mat))
    }else{
        stop('Unsupported mat type, must be matrix or numeric data.frame')
    }
    
    if (((nrow(mat) > 2000) || (ncol(mat) > 1000)) && (force_exe==FALSE)){
        stop('too many rows, set force_exe = T to continue!')
    }
    if (mat %>% apply(1,sd) %>% bazar::almost.equal (1) %>% all && mat %>% apply(1,mean) %>% bazar::almost.equal (0) %>% all){
        'already z-scored'
    }else if(zscore.rize == TRUE){
        mat = mat %>% yscale_rows()
        print('Standarizing input mat -> N(0,1)')
    }
        
    x = yget_grps(colData)
    grps = x[['grps']] # all groups
    col = x[['col']] # column name
    
    sub.heatmap = list()
    if (bygroup == TRUE){
        col_orders = c()
        for (g in grps){
            mask = colData[,col] == g
            # annotation_colors_subgrp = annotation_colors[g]
            submat = mat[,mask]
            subcolData = colData %>% dplyr::filter(mask)
            res = pheatmap::pheatmap(submat
                                    ,color = color
                                    ,annotation = colData
                                    ,cluster_rows = cluster_rows
                                    ,cluster_cols = cluster_cols
                                    ,clustering_method = clustering_method
                                    ,fontsize = fontsize
                                    ,cellwidth = cellwidth
                                    ,cellheight = cellheight
                                    ,border_color = border_color
                                    ,treeheight_row = treeheight_row
                                    ,annotation_colors = anno_color
                                    ,show_rownames = show_rownames
                                    ,show_colnames = show_colnames
                                    ,breaks = breaks
                                    ,fontsize_col = fontsize_col
                                    ,...
                                    )    
            .labels = res$tree_col$labels 
            sub.heatmap[[g]] = res
            col_orders = c(col_orders,.labels)
            
        }
        
        cluster_cols = FALSE
        print('bygroup is TRUE, set cluster_cols=FALSE')
        mat = mat[,col_orders]
    }

    res=pheatmap::pheatmap(mat
            ,color = color
            ,annotation = colData
              
            ,cluster_rows = cluster_rows
            ,cluster_cols = cluster_cols
            ,clustering_method = clustering_method
            ,fontsize = fontsize
            ,cellwidth = cellwidth
            ,cellheight = cellheight
            ,border_color = border_color
            ,treeheight_row = treeheight_row
            ,annotation_colors = anno_color
            ,show_rownames = show_rownames
            ,show_colnames = show_colnames
            ,breaks = breaks
            ,fontsize_col = fontsize_col
            ,...
        )
    
    if (bygroup == FALSE){
        col_orders = res$tree_col$labels
    }
    return(list(mat=mat,gg=res,col.orders=col_orders,sub.heatmap=sub.heatmap))

}
#-----------------

yplot_volcano_using_de = function(diff_expr, value_var = 'padj',export=FALSE,threshold=c(-2,2),p=0.05,add_label=FALSE,...){
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
    
    # value_var = old_value_var
    g %>% ydumpto(export = export, ...)
    
    list(gg=g,de=de)
}
#-----------------

yrun = function(key=NULL, args = NULL, p=NULL, .para_list=para_list){
    if (!is.null(key)) p = .para_list[[key]]
    stopifnot(!(is.null(key)&&is.null(p))) # stop if you do not provide key and p

    argv = p$args
    if (!is.null(args)){
        argv = argv %>% ypush(args)
    }
    for (i in p$requires){
        library(package = i, character.only = T)
    }
    
    func = get(p$func,envir = .GlobalEnv)    
    
    func %>% do.call(args=argv)
}
#-----------------

library(maftools)
yplot_oncoplot = function(laml=NULL, maf_=NULL, gl_=NULL
                        ,frm='export' #
                        ,only_12_cols = TRUE
                        
                        ,worker=NULL #
                        ,flag = 'oncoplot' #
                        ,outputdir='./export' #
                        ,export=FALSE #
                        ,putdatatable=F
                        
                        ,unique_on='default'
                        ,color_palette_type='snv'
                        ,height=NULL,width=NULL
                        ,retmode= 'laml'

                        ,top = 20
                        ,clinicalFeatures=c('Clin_classification')
                        ,annotationColor=NULL
                        ,removeNonMutated=F
                        ,...){
    # retmode if 'laml' return laml = read.maf(maf_,gl_)
    #         if NULL, or others return NULL
    # using laml instead of maf_,gl_ on cnv plot will cause an error @2022-03-08
    if (is.null(laml) && (!is.null(maf_)&&!is.null(gl_))){
        maf_ = yload_dfx(maf_, frm = frm)
        if (only_12_cols==TRUE){
            maf_ = maf_[,1:12]
        }
        gl_ = yload_dfx(gl_,frm = frm)
        laml = read.maf(maf=maf_,clinicalData=gl_)
    }else if(!is.null(laml)){
        stopifnot(laml %>% class == "MAF")
        gl_ = laml %>% getClinicalData
        
    }else{
        stop ("maf_,gl_,laml all NULL")
    }
    # <drop> drop duplicated rows to eliminate false 'Multi_Hit'
    # NOTE: <drop> will affect <getname> results, makes it to put contents instead of name, so must do <drop> after <getname>
    if (unique_on=='default' && color_palette_type=='cnv' ){
        unique_on = c('Hugo_Symbol','Tumor_Sample_Barcode','Variant_Classification')
        mask = duplicated(maf_[,unique_on])
        maf_ = maf_[!mask,]
        # especially removes false 'Multi_Hit' created by paralog genes like Histone genes
        # which have the same name but locate at different genome locations treated by maftools as same gene with Multi_Hit
    }
    # </drop>


    if (is.null(height)){
        height = round(top/3,1.5)
    }
    if (is.null(width)){
        width = round(nrow(gl_)/5,1.5)
    }
    if (color_palette_type == 'cnv'){
        colors=COL_CNV
    }else{
        colors=COL_SNV
    }
       
    argv = list(maf=laml,top=top, colors = colors
             ,removeNonMutated=removeNonMutated
             ,clinicalFeatures = clinicalFeatures
             ,sortByAnnotation = TRUE
             ,annotationColor = annotationColor
             ,writeMatrix = putdatatable
            )
    argv <- argv %>% ypush(list(...),.expand = TRUE) 
    

    oncoplot %>% do.call(args=argv)
    if (export != FALSE){
        oncoplot %>% ydumpto(export=export,outputdir = outputdir,
                         flag=flag,worker=worker,
                         args=argv,
                         width=width,height=height,
                         ...)
    }
    if (putdatatable==TRUE){
        onco_matrix = yload_dfx('./onco_matrix.txt')
        onco_matrix %>% ydumpto(export=putdatatable, suffix=export,flag='oncoplot_data',
#                               prefix = color_palette_type,
                              worker=worker,
                              outputdir=outputdir,
                            ext='txt',
                            row.names=T
                             )
    }
    if (retmode =='laml') return(laml)
    else return (argv)
}

#-----------------

yplot_oncoplot2 = function(laml=NULL, maf_=NULL, gl_=NULL
                        # custom func params
                        ,color_palette='snv'
                        ,unique_on='default'
                        ,putdatatable=F
                        ,na_color = NULL
                        ,wt_color = NULL
                        ,only_12_cols = TRUE
                        ,random_cc_color = FALSE
                        ,retmode= 'laml'
                        # output params
                        ,export=FALSE #
                        ,worker=NULL #
                        ,flag = 'oncoplot' #
                        ,outputdir=OUTPUTROOT
                        ,height=NULL,width=NULL
                        # set default params
                        ,top = 20
                        #,genes = genes
                        ,sortByAnnotation=TRUE#
                        ,clinicalFeatures=c('Clin_classification')
                        ,removeNonMutated=FALSE
                        ,...){
    # retmode if 'laml' return laml = read.maf(maf_,gl_)
    #         if NULL, or others return NULL
    # using laml instead of maf_,gl_ on cnv plot will cause an error @2022-03-08
    
    laml = yload_laml_maf(maf_ = maf_,gl_ = gl_,laml = laml,only_12_cols = T)
    
    gl_ = laml %>% getClinicalData
    # <drop> drop duplicated rows to eliminate false 'Multi_Hit'
    # NOTE: <drop> will affect <getname> results, makes it to put contents instead of name, so must do <drop> after <getname>
    if (unique_on=='default' && color_palette=='cnv' ){
        unique_on = c('Hugo_Symbol','Tumor_Sample_Barcode','Variant_Classification')
        mask = duplicated(maf_[,unique_on])
        maf_ = maf_[!mask,]
        # especially removes false 'Multi_Hit' created by paralog genes like Histone genes
        # which have the same name but locate at different genome locations treated by maftools as same gene with Multi_Hit
    }
    # </drop>
    argv = yget_args(...)
    formal = argv %>% names
    
    # set default dot color palette
    if (!'colors' %in% formal){
        cols = yget_gc_anno_colors(color_palette)
        argv$colors = cols
    }
    
    if (is.null(height)){
        argv$height = round(top/3,1.5)
    }
    if (is.null(width)){
        argv$width = round(nrow(gl_)/5,1.5)
    }

    if (!('annotationColor' %in% formal)){
        lengend_color = list()
        for (g in clinicalFeatures){
            clsnames = gl_[[g]] %>% unique
            print(clsnames)
            cols = yroll_colors_n(clsnames %>% length, 2,random = random_cc_color)
            lengend_color[[g]] = cols
            names(lengend_color[[g]]) = clsnames
            if ((NULL %in% clsnames) && (!is.null(na_color))){
                lengend_color[[g]][lengend_color[[g]] %>% names %>% is.null] = na_color
            }
            # lengend_color[[g]][lengend_color[[g]] %>% names %>% is.null] = wt_color
            if (('WT' %in% clsnames) && (!is.null(wt_color))){
                lengend_color[[g]][['WT']] = wt_color
            }
            argv$annotationColor = lengend_color
        }
    }
    argv$maf = laml
    argv = argv %>% ygetlast(c('maf',formals(oncoplot) %>% names))
    gg = oncoplot %>% do.call(args=argv)
    oncoplot %>% ydumpto(export=export,outputdir = outputdir,
                         flag=flag,worker=worker,
                         args=argv,
                         width=width,height=height,
                         ...)
    if (putdatatable==TRUE){
#         r = ydumpto(NULL,worker=worker,outputdir=outputdir,export=export,...)
#         r$ext='txt'
        onco_matrix = yload_dfx('./onco_matrix.txt')
        onco_matrix %>% ydumpto(export=putdatatable, suffix=export,flag='oncoplot_data',
#                               prefix = color_pallete_type,
                              worker=worker,
                              outputdir=outputdir,
                            ext='txt',
                            row.names=T
                             )
#         out = file.path(outputdir,paste0('onco_matrix_',export,'.dfx'))
#         write(paste0('[R] write 1 dfx at ',out),2)
#         file.rename(from = './onco_matrix.txt',  to = out)
    }
    if (retmode =='laml')  return(laml)
    else return (argv)
}

#-----------------

yplot_forestplot = function(laml=NULL,maf_=NULL,gl_=NULL,
                             include_total=F,
                             pVal = 0.1,
                             export = F,
                             outputdir='export',
                             minMut = 1,
                             worker=NULL,
                             putdatatable=F,
                             ret_data = T,
                             ret_data_type ='list' # {'df','list'}
                             ,...){
    res = list()
    laml = yload_laml_maf(maf_,gl_,laml)
    laml %>% yget_subset_maf(include_total = include_total) %>% combn(2, simplify = F) %>% map(
        function(i){
            n = names(i)        
            xn = n[[1]]
            yn = n[[2]]
            x = i[[xn]]
            y = i[[yn]]
            vs = mafCompare(x,y,xn,yn,minMut=minMut)
            argv = list(mafCompareRes=vs, pVal = pVal,...)
            label = paste0(xn,'_vs_',yn)
            ylog('[INFO] plot ',label)
            if (putdatatable==TRUE){
                vs$results %>% ydumpto(export=export
                                   ,flag='forest_data'
                                   ,suffix=paste0(xn,'_vs_',yn)
                                   ,worker=worker
                                   ,outputdir=outputdir
                                   ,ext='txt'
                                   ,...
                                  )
            }
            forestPlot %>% do.call(args=argv)
            forestPlot %>% ydumpto(export=export
                                   ,flag='forest'
                                   ,suffix=paste0(xn,'_vs_',yn)
                                   ,args=argv
                                   ,worker=worker
                                   ,outputdir=outputdir
                                   ,...
                                  )
            
            res[[label]] <<- vs$results 
        }
    )
    if (ret_data==TRUE){
        return(res)
    }
}

#-----------------

yplot_coOncoplot = function(laml=NULL,maf_=NULL,gl_=NULL,
                             include_total=F,
                             pVal = 0.1,export = F,
                             genes = NULL, removeNonMutated = FALSE,
                              worker='A0'
                             ,...){
    laml = yload_laml_maf(maf_,gl_,laml)
    laml %>% yget_subset_maf(include_total = include_total) %>% combn(2, simplify = F) %>% map(
        function(i){
            n = names(i)        
            xn = n[[1]]
            yn = n[[2]]
            x = i[[xn]]
            y = i[[yn]]
            ylog('[INFO] plot ',xn,' vs ',yn)
            argv = list(m1=x,m2=y,m1Name=xn,m2Name=yn,genes=genes,removeNonMutated=removeNonMutated,...)
            res = coOncoplot %>% do.call(argv)
            res = coOncoplot %>% ydumpto(export=export
                                         ,flag = 'coOncoplot'
                                         ,args=argv
                                         ,worker=worker
                                        )
        }
    )
}
#-----------------

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
#-----------------

yplot_TiTvplot = function (laml = NULL, maf_ = NULL, gl_ = NULL, export = F, plot_subgrp = T, include_total = T, 
    outputdir = outputdir, worker = NULL, putdatatable = F, sample_order_df = NULL, 
    ...) {
    # sample_order_matrix, like glx/tmb generated by tmb barplot, 
    # contains ['g','Custom_Label'] ordered sample id named 'Custom_Label', where 'g' is used to apply order to subgroups
    
    laml = yload_laml_maf(maf_, gl_, laml)
    
    res =  laml %>% yget_subset_maf(include_total = include_total) %>% 
        ymap(function(x, n) {
            print(n)
            laml.titv = titv(maf = x, plot = FALSE, useSyn = TRUE)
            ylog(sample)
            ylog("[INFO] plot ", n)
            argv = list(res = laml.titv)
            if (!(sample_order_df %>% is.null)) {
                sampleOrder = sample_order_df[sample_order_df$g == 
                  n, "Custom_Label"]
            }
            else {
                sampleOrder = NULL
            }
            plotTiTv %>% ydumpto(export = n, flag = "TiTv", 
                outputdir = outputdir, args = argv, worker = worker, 
                sampleOrder = sampleOrder, ...)
            if (putdatatable == TRUE) {
                for (i in laml.titv %>% names) {
                  laml.titv[[i]] %>% ydumpto(export = i %>% str_replace(fixed("."), 
                    "_"), ext = "txt", flag = "titv_data", prefix = n, 
                    outputdir = outputdir, worker = worker)
                }
            }
            plotTiTv %>% R.utils::doCall(args = argv, sampleOrder = sampleOrder, 
                ...)
            laml.titv
        })
    res
}

#-----------------

yplot_lollipopplot = function(laml=NULL,maf_=NULL,gl_=NULL,export = F
        ,include_total=T
        ,AACol ='HGVSp_Short'
        ,labelPos=NULL
        ,top=5
        ,genes=NULL
        ,showMutationRate =T
        ,outputdir = OUTPUTROOT
        ,worker=WORKER
        ,...){
    laml = yload_laml_maf(maf_,gl_,laml)
    plot_func = lollipopPlot
    argv =  as.list(environment()) %>% ygetlast(formals(plot_func) %>% names) %>% ypush(list(maf = laml))

    if (is.null(genes)){
        laml %>% getGeneSummary %>% head(top) %$% Hugo_Symbol %>% yloop(function(g){
            argv$gene = g
            res = plot_func %>% ydumpto(export=export
                                        ,suffix=g
                                         ,flag = 'lollipop'
                                         ,args=argv
                                        ,worker=worker
                                        )
            res = plot_func %>% do.call(argv)
        })
    }else{
        genes %>% yloop(function(g){
            argv$gene = g
            res = plot_func %>% ydumpto(export=export
                                        ,suffix=g
                                         ,flag = 'lollipop'
                                        ,outputdir=outputdir
                                         ,args=argv
                                        ,worker=worker
                                        )
            res = plot_func %>% do.call(argv)
            NULL
        })
    }
}
#-----------------

yplot_somaticInteractions = function(laml=NULL,maf_=NULL,gl_=NULL,export = F
                                        ,outputdir = 'export'
                                        ,include_total=T
                                        ,top=25
                                        ,pvalue = c(0.05, 0.1) 
                                      ,worker=WORKER
        ,...){
    
    laml = yload_laml_maf(maf_,gl_,laml)
    plot_func = somaticInteractions
    argv =  as.list(environment()) %>% ygetlast(formals(plot_func) %>% names) 
    laml %>% yget_subset_maf(include_total = include_total) %>% yloop(function(x,n){
        argv$maf = x
        res = plot_func %>% ydumpto(export=export
                                ,prefix=n
                                 ,flag = 'somaticInteractions'
                                ,outputdir=outputdir
                                 ,args=argv
                                ,worker=worker
                                )
        res = plot_func %>% do.call(argv)
        ylog('plot ',n)
        
        NULL 
    })
}
#-----------------

# functions suitable:  [somaticInteractions, lollipopPlot, oncoplot, drugInteractions, OncogenicPathways,
#               PlotOncogenicPathways, ]

yplot_maf_byfunc = function(laml=NULL,maf_=NULL,gl_=NULL,export = F
                    ,plot_func=NULL
                    ,outputdir = 'export'
                    ,include_total = T
                    ,only_total = T    
                    ,worker="A0"
        ,...){
    # plot_func in [somaticInteractions, lollipopPlot, oncoplot, drugInteractions, OncogenicPathways,
    #               PlotOncogenicPathways, ]
    laml = yload_laml_maf(maf_,gl_,laml)
    
    if (is.name(plot_func)) {
        func_name = deparse(substitute(plot_func))
    } else if (plot_func %>% is.character) {
        func_name = plot_func
        plot_func = get(func_name,envir = .GlobalEnv)
    } else if (is.function(plot_func)){
        func_name = 'plot_func'
    }
    
    argv =  as.list(match.call()) %>% ygetlast(formals(plot_func) %>% names) 
    
    if (is.name(plot_func)) {
        func_name = deparse(substitute(plot_func))
    } else if (plot_func %>% is.character) {
        func_name = plot_func
        plot_func = get(func_name,envir = .GlobalEnv)
    } else if (is.function(plot_func)){
        func_name = 'plot_func'
    }
    
    if (only_total==TRUE){
        argv$maf = laml
        argv = argv %>% ygetlast(names(formals(plot_func)))
#         return(list(plot_func,argv,names(formals(plot_func)))) 
        res = plot_func %>% ydumpto(export=export
                                 ,flag = func_name
                                    ,suffix="ALL"
                                ,outputdir=outputdir
                                 ,args=argv
                                ,worker=worker
                                )
        res = plot_func %>% do.call(argv)
        ylog('plot ',func_name,', group: ALL')
    }else{
        laml %>% yget_subset_maf(include_total = include_total) %>% yloop(function(x,n){
            argv$maf = x 
            argv = argv %>% ygetlast(names(formals(plot_func)))
#             return(list(plot_func,argv))            
            res = plot_func %>% ydumpto(export=export
                                    ,suffix=n
                                     ,flag = func_name
                                    ,outputdir=outputdir
                                     ,args=argv
                                    ,worker=worker
                                    )
            res = plot_func %>% do.call(argv)
            ylog('plot ',func_name,' group: ',n)
            NULL
        })
    }
    
}
#-----------------

# functions suitable:  [coBarplot, coOncoplot, lollipopPlot2, mafCompare]
yplot_vs_maf_byfunc = function(laml=NULL,maf_=NULL,gl_=NULL,export = F
                    ,plot_func = NULL
                    ,outputdir = OUTPUTROOT
                    ,include_total = F
                    ,top=NULL
                    ,worker=NULL
        ,...){
    # plot_func in [somaticInteractions, lollipopPlot, oncoplot, drugInteractions, OncogenicPathways,
    #               PlotOncogenicPathways, ]
    laml = yload_laml_maf(maf_,gl_,laml)
    
    if (plot_func %>% is.character) {
        func_name = plot_func
        plot_func = get(func_name,envir = .GlobalEnv)
    } else if (is.name(plot_func)) {
        func_name = deparse(substitute(plot_func))
    } else if (is.function(plot_func)){
        func_name = substitute(plot_func) |> as.character()
    } else{
        stop('wrong argument @func_name type')
    }
    argv =  yget_args(...)
    res = laml %>% 
            yget_subset_maf(include_total = include_total) %>% 
            combn(2,simplify = F) %>% 
            yloop(function(li){
                print(li %>% class)
                assign(x = 'li',value = li, envir = .GlobalEnv)
                ns = names(li)
                xn = ns[[1]]
                yn = ns[[2]]
                x = li[[xn]]
                y = li[[yn]]
                tag = ns |> str_flatten(collapse = '_vs_')

                if (func_name =='coOncoplot' && !(is.null(top))){
                    m1.genes = getGeneSummary(x)[1:top]
                    m2.genes = getGeneSummary(y)[1:top]
                    mdt = merge(m1.genes[, .(Hugo_Symbol, MutatedSamples)], 
                        m2.genes[, .(Hugo_Symbol, MutatedSamples)], by = "Hugo_Symbol", 
                        all = TRUE)
                    mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
                    mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
                    mdt$max = apply(mdt[, .(MutatedSamples.x, MutatedSamples.y)], 
                        1, max)
                    mdt = mdt[order(max, decreasing = TRUE)]
                    genes = mdt[, Hugo_Symbol]
                    argv = argv |> ypush(list(genes=genes))
                }
                argv1 = argv
                argv = argv %>% 
                        ypush(list(m1=x,m2=y,m1Name=xn,m2Name=yn,m1_name=xn,m2_name=yn),.expand = T) %>% 
                        ygetlast(plot_func %>% formals %>% names)
                # print(names(argv))
                print(paste0('Included groups ',tag))
                # return (list(plot_func,argv))
                res = plot_func %>% ydumpto(export=export
                                        ,suffix= tag
                                         ,flag = func_name
                                        ,outputdir=outputdir
                                         ,args=argv
                                        ,worker=worker
                                        )
                res = plot_func %>% do.call(argv)
        #         ylog('plot ',func_name,' group: ',xn,' vs ',yn)
            })
#     submafs
    res
}
#-----------------

plot_cosmic_inner<-function(laml,sig_ver='v2',export=FALSE,putdatatable=FALSE,outputdir='export',worker=NULL,...){
# official example has an additional parameter prefix='chr' which must be removed on our dataset
    laml.tnm = trinucleotideMatrix(maf = laml, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
#    laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
#     plotCophenetic(res = laml.sign)
    laml.tnm$nmf_matrix = laml.tnm$nmf_matrix + 0.001
    laml.sig = extractSignatures(mat = laml.tnm, n = 3)
    # sometimes after extractSignatures the laml.sig$signatures loss rownames which will cause the raise of
    # Error in crossprod(sig, x): non-conformable arguments
    rownames(laml.sig$signatures)<-colnames(laml.tnm$nmf_matrix)
    if (sig_ver == 'v2'){
        ## Compate against original 30 signatures 
        laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
        laml.cosm = laml.og30.cosm
    }else if (sig_ver == 'v3'){
         ## Compate against updated version3 60 signatures 
        laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")
        laml.cosm = laml.v3.cosm
    }else{
        stop('[sig_ver] Value Error!')
    }
    
   #     xx<-pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
    xx<-pheatmap(mat = laml.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
    ## export cosine_similarity matrix
    if (putdatatable!=FALSE){
        laml.cosm$cosine_similarities %>% ydumpto(export='data',
                                                  outputdir=outputdir,worker=worker,...)
    }
    xx %>% ydumpto(export=export,
                   outputdir=outputdir,
                   worker=worker,
                   flag='sig_vs_cosmic',
                   prefix=sig_ver,
                  )
    if (sig_ver == 'v2'){
        # plot against original 30 signatures 
        plotSignatures(nmfRes = laml.sig,title_size = 1.2, sig_db = "legacy")
        argv = list(nmfRes = laml.sig,
                   title_size = 1.2, 
                   sig_db = "legacy")
        plotSignatures %>% ydumpto(export=export,
                                   outputdir=outputdir,
                                   worker=worker,
                                   flag='signature',
                                   prefix=sig_ver,
                                   args=argv
                                  )
    }else 
    if (sig_ver == 'v3'){
         # plot against updated version3 60 signatures 
        plotSignatures(nmfRes = laml.sig,title_size = 1.2, sig_db = "SBS")
        argv = list(nmfRes = laml.sig,
                   title_size = 1.2, 
                   sig_db = "SBS")
        plotSignatures %>% ydumpto(export=export,
                                   outputdir=outputdir,
                                   worker=worker,
                                   flag='signature',
                                   prefix=sig_ver,
                                   args=argv
                                  )
    }
    
#     put.func(paste0(fname,'_similarity'),pheatmap::pheatmap,
#              mat = laml.og30.cosm$cosine_similarities, 
#              cluster_rows = FALSE, main = "cosine similarity against validated signatures")
#     dump.dfx(laml.sig,paste0('cosmic_',fname,'_signatures'))
    
    
}

yplot_cosmic <- function(laml=NULL,maf_=NULL,glx_=NULL
                          ,only_plot_subgroups=FALSE
                          ,putdatatable=TRUE
                          ,sig_ver='v2'
                          ,export=FALSE
                          ,outputdir=OUTPUTROOT,worker=NULL,
                          ...){
    if (is.null(laml)){
        laml = yload_laml_maf(maf_,glx_)
    }else{
        glx_ = laml %>% getClinicalData
    }
    grps = glx_$Clin_classification %>% unique
    print(grps)
    if (length(grps) ==1){
        only_plot_subgroups = FALSE
    }
    if (only_plot_subgroups){
        res = list()
        for (g in grps){
            res[[g]]= subsetMaf(laml,clinQuery=paste0('Clin_classification=="',g,'"'))  
            if (export!=FALSE) export=names(res[g])
            plot_cosmic_inner(res[[g]],export=export,flag='CosSig'
                              ,prefix=sig_ver
                              ,putdatatable=putdatatable
                              ,suffix=g
                              ,sig_ver=sig_ver
                              ,outputdir=outputdir
                              ,worker=worker
                             )
            
        }
    }else{
        plot_cosmic_inner(laml,export=export
                          ,prefix=sig_ver
                          ,flag='CosSig'
                          ,sig_ver=sig_ver
                          ,putdatatable=putdatatable
                          ,outputdir=outputdir
                          ,worker=worker)
        
    }
       
}

#-----------------

library("BSgenome.Hsapiens.UCSC.hg19")
library('MutationalPatterns')

#-----------------

yget_cosmic_signature_ref = function(cosmic_ver='v2'){
    if (cosmic_ver=='v2'){
        # get cosimic_cancer_signatrues v2 30 sigs
        sp_url <- '/GCI/jup/A_TSF/libs/signatures_probabilities.txt'
        cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
        # Match the order of the mutation types to MutationalPatterns standard
        new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
        # Reorder cancer signatures dataframe
        cancer_signatures = cancer_signatures[as.vector(new_order),]
        # Add trinucletiode changes names as row.names
        row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
        # Keep only 96 contributions of the signatures in matrix
        cancer_signatures = as.matrix(cancer_signatures[,4:33])
    }else if (cosmic_ver =='v3'){
        # get cosmic_caner_signatures v3 78 sigs
        sp_url = '/GCI/jup/A_TSF/libs/COSMIC_v3.2_SBS_GRCh37.txt'
        cancer_signatures2 = read.table(sp_url, sep = "\t", header = TRUE)

        new_order = match(row.names(mut_mat), cancer_signatures2$Type)
        # Reorder cancer signatures dataframe
        cs = cancer_signatures2[as.vector(new_order),]
        #         cs = cancer_signatures2
        # Keep only 96 contributions of the signatures in matrix
        row.names(cs) = row.names(mut_mat)
        cancer_signatures = cs[,-1]
    }else{
        stop("wrong COSMIC_VER value, must in ['v2','v3']")
    }
    cancer_signatures
}

yget_signature_contribution <- function(laml, cancer_signatrues){
        laml.tnm = trinucleotideMatrix(maf = laml
    #                                , prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
        mut_mat = t(laml.tnm$nmf_matrix)

        mut_mat = data.matrix(mut_mat)
        cancer_signatrues = data.matrix(cancer_signatrues)

        fit_res <- fit_to_signatures(data.matrix(mut_mat), cancer_signatrues)
        fit_res
        # list("contribution"=fit_res$contribution,"plot"=g)
}

yplot_sig_contribution = function(laml=NULL
                                   ,maf_=NA,glx_=NA
                                   ,export=FALSE
                                   ,outputdir=OUTPUTROOT
                                   ,width=12,height=6
                                   ,plot=TRUE
                                   ,cosmic_ver='v2'
                                  ,...){

    laml = yload_laml_maf(maf_ =maf_,gl_ = glx_,laml=laml)
    
    if (!is.null(laml)){
        glx = getClinicalData(laml)
    }
    grps = glx$Clin_classification %>% unique
    print(paste(' ----- grps is', str_flatten(grps,'_'),' -------'))

    # <getsig> get cosmic signature by v2 or v3
    laml.tnm = trinucleotideMatrix(maf = laml
    #                                , prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = t(laml.tnm$nmf_matrix)
    if (cosmic_ver=='v2'){
        # get cosimic_cancer_signatrues v2 30 sigs
        sp_url <- '/GCI/jup/A_TSF/libs/signatures_probabilities.txt'
        cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
        # Match the order of the mutation types to MutationalPatterns standard
        new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
        # Reorder cancer signatures dataframe
        cancer_signatures = cancer_signatures[as.vector(new_order),]
        # Add trinucletiode changes names as row.names
        row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
        # Keep only 96 contributions of the signatures in matrix
        cancer_signatures = as.matrix(cancer_signatures[,4:33])
    }else if (cosmic_ver =='v3'){
        # get cosmic_caner_signatures v3 78 sigs
        sp_url = '/GCI/jup/A_TSF/libs/COSMIC_v3.2_SBS_GRCh37.txt'
        cancer_signatures2 = read.table(sp_url, sep = "\t", header = TRUE)

        new_order = match(row.names(mut_mat), cancer_signatures2$Type)
        # Reorder cancer signatures dataframe
        cs = cancer_signatures2[as.vector(new_order),]
        #         cs = cancer_signatures2
        # Keep only 96 contributions of the signatures in matrix
        row.names(cs) = row.names(mut_mat)
        cancer_signatures = cs[,-1]
    }else{
        stop("wrong COSMIC_VER value, must in ['v2','v3']")
    }
    # </getsig>


    get_signature_contribution <- function(laml,cancer_signatrues,plot=TRUE){
        laml.tnm = trinucleotideMatrix(maf = laml
    #                                , prefix = 'chr'
                                   , add = TRUE
                                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
        mut_mat = t(laml.tnm$nmf_matrix)

        mut_mat = data.matrix(mut_mat)
        cancer_signatrues = data.matrix(cancer_signatrues)

        fit_res <- fit_to_signatures(data.matrix(mut_mat), cancer_signatrues)
        if (plot==TRUE){
            gg = plot_contribution(fit_res$contribution,
              coord_flip = FALSE,
              mode = "relative"
            )
            print(gg)
        }else{
            gg=NA
        }
        fit_res$plot = gg
        fit_res
        # list("contribution"=fit_res$contribution,"plot"=g)
    }

    grps = glx[,'Clin_classification'] %>% unique()
    
    res = list()
    
    for (g in grps){
        sub_maf = subsetMaf(laml,clinQuery = paste0("Clin_classification == '",g,"'"))
        ## get absolute contribution by MutationalPatterns::fit_to_signatures
        resg = get_signature_contribution(sub_maf, cancer_signatures,plot=plot)
        ## gen relative contribution over each patient(column)
        resg$contribution_rel = resg$contribution %>% 
                                data.frame %>% 
                                mutate(across(everything(), ~ . * 100 / sum(.))) %>% 
                                as.matrix
        res[[g]] = resg
    }
    res
}
#-----------------

yget_known_cancer_signatures = function(cosmic_ver='v2'){
    tnn_order = c('A[C>A]A','A[C>A]C','A[C>A]G','A[C>A]T','C[C>A]A','C[C>A]C','C[C>A]G','C[C>A]T','G[C>A]A','G[C>A]C','G[C>A]G','G[C>A]T','T[C>A]A','T[C>A]C','T[C>A]G','T[C>A]T','A[C>G]A','A[C>G]C','A[C>G]G','A[C>G]T','C[C>G]A','C[C>G]C','C[C>G]G','C[C>G]T','G[C>G]A','G[C>G]C','G[C>G]G','G[C>G]T','T[C>G]A','T[C>G]C','T[C>G]G','T[C>G]T','A[C>T]A','A[C>T]C','A[C>T]G','A[C>T]T','C[C>T]A','C[C>T]C','C[C>T]G','C[C>T]T','G[C>T]A','G[C>T]C','G[C>T]G','G[C>T]T','T[C>T]A','T[C>T]C','T[C>T]G','T[C>T]T','A[T>A]A','A[T>A]C','A[T>A]G','A[T>A]T','C[T>A]A','C[T>A]C','C[T>A]G','C[T>A]T','G[T>A]A','G[T>A]C','G[T>A]G','G[T>A]T','T[T>A]A','T[T>A]C','T[T>A]G','T[T>A]T','A[T>C]A','A[T>C]C','A[T>C]G','A[T>C]T','C[T>C]A','C[T>C]C','C[T>C]G','C[T>C]T','G[T>C]A','G[T>C]C','G[T>C]G','G[T>C]T','T[T>C]A','T[T>C]C','T[T>C]G','T[T>C]T','A[T>G]A','A[T>G]C','A[T>G]G','A[T>G]T','C[T>G]A','C[T>G]C','C[T>G]G','C[T>G]T','G[T>G]A','G[T>G]C','G[T>G]G','G[T>G]T','T[T>G]A','T[T>G]C','T[T>G]G','T[T>G]T')
    if (cosmic_ver=='v2'){
        # get cosimic_cancer_signatrues v2 30 sigs
        sp_url <- '/GCI/jup/A_TSF/libs/signatures_probabilities.txt'
        cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
        # Match the order of the mutation types to MutationalPatterns standard
        new_order = match(tnn_order, cancer_signatures$Somatic.Mutation.Type)
        # Reorder cancer signatures dataframe
        cancer_signatures = cancer_signatures[as.vector(new_order),]
        # Add trinucletiode changes names as row.names
        row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
        # Keep only 96 contributions of the signatures in matrix
        cancer_signatures = as.matrix(cancer_signatures[,4:33])
    }else if (cosmic_ver =='v3'){
        # get cosmic_caner_signatures v3 78 sigs
        sp_url = '/GCI/jup/A_TSF/libs/COSMIC_v3.2_SBS_GRCh37.txt'
        cancer_signatures2 = read.table(sp_url, sep = "\t", header = TRUE)

        new_order = match(tnn_order, cancer_signatures2$Type)
        # Reorder cancer signatures dataframe
        cs = cancer_signatures2[as.vector(new_order),]
        #         cs = cancer_signatures2
        # Keep only 96 contributions of the signatures in matrix
        row.names(cs) = tnn_order
        cancer_signatures = cs[,-1]
    }else{
        stop("wrong COSMIC_VER value, must in ['v2','v3']")
    }
    cancer_signatures
}


yget_signatures_related_all = function(laml,cosmic_ver = 'v2',rank=3,nrun=3,single_core=TRUE,...){
    "require MutationalPatterns::extract_signatures, yget_known_cancer_signatures, MutationalPatterns::fit_to_signatures"
    # get mut_mat of 96 tri-nucleotide
    laml.tnm = trinucleotideMatrix(maf = laml
                    #, prefix = 'chr'
                   , add = TRUE
                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    mut_mat = laml.tnm$nmf_matrix %>% t %>% as.matrix
    # get custom signatures
    custom_sigs = MutationalPatterns::extract_signatures(mut_mat, rank = rank, nrun = nrun, single_core = single_core,...)
    # get contributions
    cancer_signatures = yget_known_cancer_signatures(cosmic_ver=cosmic_ver) %>% data.matrix
    fit_res <- fit_to_signatures(data.matrix(mut_mat), cancer_signatures)
    list(mut_mat=mut_mat,extracted_n_signatures_suit=custom_sigs,contribution=fit_res$contribution,tnm=laml.tnm)
}

yget_maftools_signatures = function(laml_or_tnm, rank = NA, cosmic_ver = 'v2', plot=FALSE, search_rank = FALSE, nTry =6 ){
    "extract N/rank signatures of group, the input patients are treated as a whole, and calculate the cosine similarity vs COSMIC v2/v3 signatures."
    if (is.na(rank) && search_rank == FALSE ){
        print('please run with search_rank=TRUE to find the best rank first')
        return (NA)
    }
    cls = class(laml_or_tnm)
    if (length(cls)==1 && cls == 'MAF'){
        laml.tnm = trinucleotideMatrix(maf = laml_or_tnm
                    #, prefix = 'chr'
                   , add = TRUE
                   , ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
        
    }else if (cls == 'list' ){
        mask = names(laml_or_tnm) == c('nmf_matrix','APOBEC_scores')
        if (length(mask)==2 && all(mask)==TRUE){
            tnm = laml_or_tnm
        }
    }else{
        stop('wrong type input laml_or_tnm, must be a maftools::maf or a tnm by trinucleotideMatrix(maftools::maf)')
    }
    if (is.na(rank) && search_rank == TRUE){
        laml.sign = estimateSignatures(mat = tnm, nTry = nTry)
        plotCophenetic(res = laml.sign)
        return (laml.sign)
    }
    if (cosmic_ver=='v2'){
        sig_db = 'legacy'
    }else if (cosmic_ver=='v3'){
        sig_db = 'SBS'
    }
    tnm$nmf_matrix[tnm$nmf_matrix==0] = 0.01
    laml.sig = extractSignatures(mat = tnm, n = rank)
    laml.cosm = compareSignatures(nmfRes = laml.sig, sig_db = sig_db)
    if (plot ==TRUE){
        pheatmap::pheatmap(mat = laml.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
        maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = sig_db)
    }
    c(laml.sig,laml.cosm)
}

#-----------------

yplot_gistic_chrom <- function(dir,conf=95,export=FALSE,outputdir='export',worker=NULL,markBands = 'all',...){
    '@dir dir to the gistic2.0 output group folder
     markBands {"all",NULL},NULL is top5
    '
    # dir = '/GCI/home/ysun/GISTIC/result-CCA-haibo/chol_msk_2018_GISTIC2_result'
    # dir = '/GCI/home/ysun/GISTIC/result-CCA-haibo/Genomicare_GISTIC2_result'
    confidence = conf
    env = environment()

    for (name in c('all_lesions','amp_genes','del_genes')){
        assign(name,file.path(dir,paste0(name,'.conf_',confidence,'.txt')))
    }
    scores.gistic = file.path(dir,'scores.gistic')

    laml.gistic = readGistic(gisticAllLesionsFile = all_lesions
                             , gisticAmpGenesFile = amp_genes
                             , gisticDelGenesFile = del_genes
                             , gisticScoresFile = scores.gistic
                             , isTCGA = F)
    argv = yget_args(...)
    argv$gistic = laml.gistic
    argv = argv %>% ygetlast(gisticChromPlot %>% formals %>% names)
    gisticChromPlot %>% do.call(argv)
    gisticChromPlot %>% ydumpto(args=argv,export=export,outputdir=outputdir,worker=worker,...)

}
#-----------------

yplot_venns = function(df,sel_col,export=F,...){
    library(VennDiagram)
    grps = df$Clin_classification %>% unique
    DICT = list(draw.single.venn,
                draw.pairwise.venn,
                draw.triple.venn,
                draw.quad.venn,
                draw.quintuple.venn)
    plot_func = DICT[[length(grps)]]

    sets = list()
    for (i in grps){
        sets[[i]] = df %>% filter(g == !!i) %>% .[[sel_col]]%>% unique
    }
    print(sets %>% names)

    # generate params for venn plot function 
    argv = list()
    for (l in 1:length(sets)){
        for (i in combn(1:length(sets),l,simplify = F)){
            s = sets[i]
            if (l>1){
                r = Reduce(intersect,s) %>% length
                argv[[str_flatten(c('n',i))]] = r

            }else{            
                r = length(s[[1]])
                argv[[str_flatten(c('area',i))]] = r
            }
        }
    }
    argv$cross.area = argv$n12
    argv = argv %>% ypush(list(category=sets %>% names,
            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:length(sets)],
            cat.cex = 1.2,
            cex=1,
            margin = 0.05,
            ind = TRUE
        ))
    print(argv %>% names)
    
    args = list(...) %>% 
            ypush(list(x =plot_func,
                       args=argv,
                       export=export,
                        flag='venn'))
    n = args %>% names
    if (!('suffix' %in% n)){
        args$suffix = sel_col
    }
    
    res = ydumpto %>% do.call(args)            
    res
}

#-----------------

library(Rtsne)
yplot_tnse = function(normd,colData,group_col=NA,perplexity=NA,dims=2,...){
    pids = colnames(normd)
    mat = t(normd)
    if (is.na(perplexity)){
        perplexity = floor((ncol(normd) - 1) / 3)
    }
    
    res = Rtsne(t(normd),perplexity = perplexity,dims=dims,...)
    df = res$Y %>% data.frame(row.names = pids)
    colData2 = colData[pids,,drop=FALSE]
    df = cbind(df,colData2)
    res$df = df
    
    col = colnames(colData)
    if (is.na(group_col)){
        if ('Group' %in% col) group_col='Group'
        else if ('Clin_classification' %in% col) group_col='Clin_classification'
        else if ('g' %in% col) group_col = 'g'
        else {
            print("can not determine group_col in colData, stop plotting")
            return (res)
        }
    }
    group_col = sym(group_col)
    gg = df %>% ggplot(aes(x=X1,y=X2,color=!!group_col)) +
        geom_point(size = 3,alpha = 0.7)+
        ggtitle('t-SNE')

    res$gg = gg
    res
}
#-----------------

