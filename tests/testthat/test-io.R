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




    expect_true(is.ggplot(gg))
    gg %>% ydumpto('gg',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/gg.pdf"))

    expect_true('Heatmap' %in% class(hm))
    hm %>% ydumpto('hm',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/hm.pdf"))

    expect_true(is.data.frame(mtcars))
    mtcars %>% ydumpto('df.dfx',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/df.dfx"))

    # expect .dfx with rownames
    expect_true(is.matrix(mat))
    mat %>% ydumpto('matrix.dfx',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/matrix.dfx"))

    # expect .dfx without rownames
    iris %>% ydumpto('iris',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/iris.dfx"))

    # expect .dfx with rownames
    expect_true(is.matrix(mat))
    mat %>% ydumpto('matrix.txt',outputdir = OUTPUTROOT)
    expect_true(file.exists("tmp/matrix.txt"))

    file.remove('tmp/gg.pdf')
    file.remove('tmp/hm.pdf')
    file.remove('tmp/df.dfx')
    file.remove('tmp/matrix.dfx')
    file.remove('tmp/iris.dfx')
    file.remove('tmp/matrix.txt')

    lis = c('hm','kegg','im','im_jnj')
    for (i in lis){
        message('gmt',i)
        li = yload_gmt(i)
        message('list gmt',i)
        li = yload_list_gmt(i)
    }
    expect_true(TRUE)
})
