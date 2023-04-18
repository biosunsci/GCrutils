# The final cluster number is K + 1
# GEN CUSTOM DATA

gbm.mut = read.table('wes_snv_all.dfx',sep='\t',header=TRUE) %>%
    mutate(val = 1) %>%
    select(Hugo_Symbol,Tumor_Sample_Barcode,val) %>%
    distinct() %>%
    pivot_wider(id_cols = Tumor_Sample_Barcode
                , names_from = Hugo_Symbol
                , values_from = val
                , values_fill = 0) %>%
    column_to_rownames('Tumor_Sample_Barcode')
mut.rate=apply(gbm.mut,2,mean)
gbm.mut2 = gbm.mut[,which(mut.rate>0.03)]
gbm.mut2 = gbm.mut2[order(row.names(gbm.mut2)), ] %>% as.matrix()

glx = read.table('wes_glx.dfx',sep='\t',header=TRUE)
m2id = glx %>% pull(Custom_Label,name = ap)
# read seg files remove chr in the chromosome column
gbm.seg = read.table('wes_segs.dfx',sep='\t',header=TRUE) %>% filter(chr!='chrX' & chr!='chrY') %>% mutate(chr=as.integer(str_sub(chr,4,)))
gbm.seg$ap = m2id[gbm.seg$ap]
colnames(gbm.seg) = c("sample"   ,  "chromosome", "start" ,     "end",        "num.mark",   "seg.mean"  )
gbm.cn=CNregions(seg=gbm.seg
                 ,epsilon=0
                 ,adaptive=FALSE
                 # ,rmCNV=TRUE
                 # ,cnv=variation.hg18.v10.nov.2010[,3:5]
                 ,frac.overlap=0.5, rmSmallseg=TRUE,nProbes=1000)

normd = read.table('rna_normd.dfx',row.names = 1,sep='\t',header=TRUE) %>% t
expr.sd = normd %>% apply(MARGIN = 2,sd) %>% sort(decreasing = TRUE)
gbm.exp = normd[,expr.sd %>% head(5000) %>% names]

print(c(dim(gbm.mut2),dim(gbm.cn),dim(gbm.exp)))

# --------------------
# USE DEFAULT DATA
# generate table
data(gbm) # -> gbm.mut, gbm.seg, gbm.exp
mut.rate=apply(gbm.mut,2,mean)
gbm.mut2 = gbm.mut[,which(mut.rate>0.03)]
gbm.cn=CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,
                 cnv=variation.hg18.v10.nov.2010[,3:5],
                 frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
gbm.exp[1:5,1:3]

print(c(dim(gbm.mut2),dim(gbm.cn),dim(gbm.exp)))

# --------------------
# search Ks
for(k in 1:4){
    cv.fit = tune.iClusterPlus(cpus=3,dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
                               type=c("binomial","gaussian","gaussian"),K=k,n.lambda=185,
                               scale.lambda=c(1,1,1),maxiter=20)
    save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

# --------------------
# handle search results, find best K
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
    load(dir()[files[i]])
    output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
    devRatMinBIC[i] = devR[minBICid[i],i]
}

plot(1:(nK+1),
     c(0,devRatMinBIC),
     type="b",
     xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")

# determine K by highest y minial x point in the plot

# --------------------
# choose K manually
k=3
# get best clusters then plot heatmap
clusters=getClusters(output)
rownames(clusters)=rownames(gbm.exp)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]

# Select the top features based on lasso coefficient estimates for the K-cluster solution. here K+1=3
features = alist()
features[[1]] = colnames(gbm.mut2)
features[[2]] = colnames(gbm.cn)
features[[3]] = colnames(gbm.exp)
sigfeatures=alist()
for(i in 1:3){
    rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
    upper=quantile(rowsum,prob=0.75)
    sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("mutation","copy number","expression")
#print a few examples of selected features
head(sigfeatures[[1]])
head(sigfeatures[[2]])
head(sigfeatures[[3]])







