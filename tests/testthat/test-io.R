
test_that("matrix, data.frame, ggplot, Heatmap ydumpto works", {
    name = yload_symbols_all_genes_18669()
    expect_true(length(name)==18669 && is.character(name))


    cr = yload_symbols_cancer_related_genes_1190()
    expect_true(length(cr)==1190 && is.character(cr))

    data(iris)
    target = './report/iris.dfx'
    if (file.exists(target)){
        file.remove(target)
    }
    iris %>% ydumpto('./report/iris',outputdir = 'export')
    expect_true(file.exists(target))
})
