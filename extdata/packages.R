setwd('~/Downloads/')
# install.packages(
#     'BiocManager',
#     repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
#     dependencies = NA,
#     quiet = TRUE
# )

# install.packages(c('remotes', 'languageserver'),
#                  repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
#                  dependencies = NA)

# install.packages('rJava',
#                  repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
#                  dependencies = NA,
#                  quiet = TRUE)

devtools::install_git('https://github.com/rpkgs/Ipaper', upgrade = 'always')

remotes::install_local(c('./pkgs/NORMT3_1.0.4.tar.gz','./pkgs/bmm_github.tar.gz'),dependencies = NA, upgrade = FALSE)
remotes::install_github('kassambara/ggpubr')

devtools::install_github(
    c(
        'renozao/NMF@devel',
        'kassambara/survminer',
        'sinhrks/ggfortify',
        'eclarke/ggbeeswarm',
        'jokergoo/circlize',
        'genome/bmm',
        'genome/sciClone',
        'hdng/clonevol',
        'hdng/trees',

    ),
    build_vignettes = FALSE,
    upgrade = TRUE
)


install.packages(
    c(
        'rstatix',
        'showtext',
        'survivalAnalysis',
        'forestplot',
        'ggsci',
        'viridis',
        'mclust',
        'R.utils',
        'tidyverse',
        'metafor',
        'Cairo',
        'installr',
        'packcircles',
        'pheatmap',
        'htmlwidgets',
        'VennDiagram',
        'RColorBrewer'
    ),
    repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
    dependencies = NA,
    quiet = FALSE
)
deps = c(
    "rjson",
    "GenomicAlignments",
    "BiocIO",
    "restfulr",
    "rtracklayer",
    "Rsamtools",
    "GenomicFeatures",
    "rhdf5",
    "rhdf5filters",
    "ScaledMatrix",
    "irlba",
    "rsvd",
    "beachmat",
    "BSgenome",
    "VariantAnnotation",
    "pracma",
    "ggdendro",
    "ggalluvial",
    "locfit",
    "geneplotter",
    "edgeR",
    "GSEABase",
    "SingleCellExperiment",
    "sparseMatrixStats",
    "DelayedMatrixStats",
    "HDF5Array",
    "BiocSingular"
)
yinstall_pkgs = function(pkgs){
    ref = list.files('pkgs')
    unfound = c()
    for (i in pkgs){
        fd = ref[str_starts(ref,i)]
        if (length(fd)==1){
            pkg = fd[[1]]
            print(paste('----------- Installing source pkg:',pkg,'-----------'))
            remotes::install_local(file.path('pkgs',pkg), dependencies = FALSE, upgrade = FALSE)
        }else{
            unfound = c(unfound,i)
        }
    }
    unfound
}
yinstall_pkgs(deps)
yinstall_pkgs(deps[22:27])

BiocManager::install(
    c(
        #'ggplot2',
        'MutationalPatterns',
        'DESeq2',
        'pRoloc',
        'AnnotationHub',
        'GOSemSim',
        'AnnotationDbi',
        'org.Hs.eg.db',
        'hgu95av2.db',
        'clusterProfiler',
        'ConsensusClusterPlus',
        'preprocessCore',
        'sva',
        'IRanges',
        'timescape',
        'GSVA'
    )
    ,ask = TRUE
    ,force = FALSE
    ,update = FALSE
)


install.packages('openxlsx',
                 repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
                 dependencies = NA,
                 quiet = FALSE)

# remotes::install_github('ebecht/MCPcounter',ref='master', subdir='Source', build_vignettes = FALSE,  upgrade='always');
# install.packages(c('./BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz', './DPpackage_1.1-7.tar.gz', './nloptr_2.0.0.tar.gz'),  repos = NULL, ask=FALSE, type='source')


pkg2 = c(
    'minfi',
    'Illumina450ProbeVariants.db',
    'sva',
    'IlluminaHumanMethylation450kmanifest',
    'limma',
    'RPMM',
    'DNAcopy',
    'preprocessCore',
    'impute',
    'marray',
    'wateRmelon',
    'goseq',
    'plyr',
    'GenomicRanges',
    'RefFreeEWAS',
    'qvalue',
    'isva',
    'doParallel',
    'bumphunter',
    'quadprog',
    'shiny',
    'shinythemes',
    'plotly',
    'RColorBrewer',
    'DMRcate',
    'dendextend',
    'IlluminaHumanMethylationEPICmanifest',
    'FEM',
    'matrixStats',
    'missMethyl',
    'combinat',
    'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'
)

failed = yinstall_pkgs(pkg2)
print(failed)

BiocManager::install(
    c(
        'minfi',
        'Illumina450ProbeVariants.db',
        'sva',
        'IlluminaHumanMethylation450kmanifest',
        'limma',
        'RPMM',
        'DNAcopy',
        'preprocessCore',
        'impute',
        'marray',
        'wateRmelon',
        'goseq',
        'plyr',
        'GenomicRanges',
        'qvalue',
        'isva',
        'doParallel',
        'bumphunter',
        'quadprog',
        'shiny',
        'shinythemes',
        'plotly',
        'RColorBrewer',
        'DMRcate',
        'dendextend',
        'IlluminaHumanMethylationEPICmanifest',
        'FEM',
        'matrixStats',
        'missMethyl',
        'combinat'
    )
)
BiocManager::install('RefFreeEWAS')
# install.packages(
#     c(
#         './geneLenDataBase_1.30.0.tar.gz',
#         './IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz',
#         './ChAMPdata_2.26.0.tar.gz'
#     ),
#     repos = NULL,
#     ask = FALSE,
#     type = 'source'
# )
pkg3 = c('boot', 'brew', 'broom', 'bslib', 'callr', 'car', 'caret', 'class', 'cli', 'clipr', 'clue', 'cluster',
         'clusterProfiler', 'collections', 'commonmark', 'conquer', 'cpp11', 'crayon', 'curl', 'data.table', 'DBI', 'dbplyr',
         'dendextend', 'DEoptimR', 'desc', 'devtools', 'DiagrammeR', 'digest', 'doParallel', 'dplyr', 'DT', 'dtplyr', 'e1071',
         'enrichplot', 'evaluate', 'exactRankTests', 'ff', 'flexmix', 'FNN', 'fontawesome', 'forcats', 'foreach', 'foreign',
         'formatR', 'future', 'future.apply', 'gargle', 'gbm', 'gdata', 'GenomeInfoDb', 'gert', 'ggforce', 'ggfortify', 'ggfun',
         'ggnewscale', 'ggplot2', 'ggpubr', 'ggraph', 'ggrepel', 'ggsignif', 'ggtext', 'gh', 'gitcreds', 'glmnet', 'globals',
         'googlesheets4', 'gower', 'graphlayouts', 'gridtext', 'gtable', 'gtools', 'haven', 'heatmaply', 'highr', 'hms',
         'htmltools', 'htmlwidgets', 'httpuv', 'httr', 'hwriter', 'igraph', 'ipred', 'isoband', 'iterators', 'jpeg', 'jsonlite',
         'kernlab', 'km.ci', 'knitr', 'languageserver', 'lava', 'lifecycle', 'limma', 'lintr', 'listenv', 'lme4', 'lobstr',
         'lpSolve', 'lubridate', 'maftools', 'MALDIquant', 'maptools', 'markdown', 'MASS', 'mathjaxr', 'Matrix', 'MatrixModels',
         'matrixStats', 'mclust', 'metafor', 'mgcv', 'minqa', 'mixtools', 'modelr', 'MsCoreUtils', 'MSnbase', 'ncdf4', 'nlme',
         'nloptr', 'nnet', 'openssl', 'parallelly', 'patchwork', 'pbapply', 'pillar', 'pkgbuild', 'pkgload', 'plotly', 'pls',
         'plyr', 'png', 'polyclip', 'polynom', 'processx', 'progressr', 'proxy', 'ps', 'purrr', 'qap', 'quantreg', 'R.methodsS3',
         'R.oo', 'R.utils', 'randomForest', 'Rcpp', 'RcppArmadillo', 'RcppEigen', 'RCurl', 'readr', 'readxl', 'recipes', 'reprex',
         'rlang', 'rmarkdown', 'robustbase', 'roxygen2', 'rpart', 'RSQLite', 'rstatix', 'rstudioapi', 'rversions', 'rvest',
         'S4Vectors', 'sass', 'scales', 'scatterpie', 'segmented', 'seriation', 'sfsmisc', 'shadowtext', 'shiny', 'sp', 'spatial',
         'stringr', 'styler', 'survival', 'survMisc', 'sys', 'testthat', 'tidygraph', 'tidyr', 'tidyselect', 'tidytree',
         'tidyverse', 'timeDate', 'tinytex', 'TSP', 'tweenr', 'tzdb', 'uuid', 'vctrs', 'viridisLite', 'visNetwork', 'vroom',
         'webshot', 'whisker', 'xfun', 'XML', 'yaml', 'yulab.utils', 'zip', 'zoo')
failed = yinstall_pkgs(pkg3)
print(failed)
BiocManager::install(c( 'RefFreeEWAS', 'FEM', 'missMethyl', 'goseq', 'DMRcate'),ask=TRUE)

# install.packages(c('./ChAMP_2.24.0.tar.gz'),repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/',ask=FALSE,type='source')
BiocManager::install('M3C', ask = TRUE)

remotes::install_github(
    c('r-lib/systemfonts', 'teunbrand/ggh4x'),
    build_vignettes = FALSE,
    upgrade = 'always'
)
remotes::install_github('biosunsci/maftools'
                        ,build_vignettes = FALSE
                        ,upgrade = FALSE
                        ,force = TRUE)
install.packages(
    c('magick', 'bazar'),
    repos = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/',
    dependencies = NA,
    quiet = TRUE
)


#2023-01-03
install.packages('fastcluster')
install.packages('rlist')
