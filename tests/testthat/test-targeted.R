test_that("we can grab gene names from go terms", {

    expected.go.09 <- c("PIGV","ALG12")
    actual.go.09 <- suppressMessages(get_go_term_human_genes("GO:0000009"))

    expect_equal(actual.go.09, expected.go.09)
})


test_that("we can grab KOs from GMM matrix", {
    GMM.matrix <- get_GMM_matrix()
    mtx.feature.names <- c("k02005_hly_d_family_secretion_protein",
                           "k02007_cobalt_nickel_transport_system_permease_protein",
                           "k00248_butyrate_gene_1",
                           "K00634_butyrate_gene_2",
                           "K00929_butyrate_gene_3")

    actual.features <- features_from_gmm_df(GMM="butyrate production I",
                                            GMM.kos.df=GMM.matrix,
                                            mtx.feature.names=mtx.feature.names)
    expected.features <- c("k00248_butyrate_gene_1",
                           "K00634_butyrate_gene_2",
                           "K00929_butyrate_gene_3")

    expect_equal(actual.features, expected.features)

})
