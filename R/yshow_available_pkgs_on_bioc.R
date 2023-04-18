yshow_available_pkgs_on_bioc = function(){
    df <- read.dcf(
        url("https://bioconductor.org/packages/release/bioc/src/contrib/PACKAGES")
    )
    df = as.data.frame(df)
    df
}
