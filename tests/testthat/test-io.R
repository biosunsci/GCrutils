data(iris)
data(mtcars)

mat = as.matrix(mtcars)
colnames(iris)
gg = ggplot(iris,aes(Sepal.Length,Sepal.Width)) + geom_point()
make.custom(5,5)
hm = ComplexHeatmap::Heatmap(mat,cluster_rows = F,cluster_columns = F)

OUTPUTROOT = 'tmp/'


test_that("matrix, data.frame, ggplot, Heatmap ydumpto works", {
    expect_equal(OUTPUTROOT,'tmp/')


    expect_true(is.data.frame(mtcars))
    mtcars %>% ydumpto('df',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/df.dfx"))

    expect_true(is.ggplot(gg))
    gg %>% ydumpto('gg',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/gg.pdf"))

    # expect .dfx with rownames
    expect_true(is.matrix(mat))
    mat %>% ydumpto('matrix',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/matrix.dfx"))

    expect_true('Heatmap' %in% class(hm))
    hm %>% ydumpto('hm',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/hm.pdf"))

    # expect .dfx without rownames
    iris %>% ydumpto('iris',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/iris.dfx"))

})
