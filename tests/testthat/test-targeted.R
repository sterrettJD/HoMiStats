test_that("we can grab gene names from go terms", {

    expected.go.09 <- c("PIGV","ALG12")
    actual.go.09 <- suppressMessages(get_go_term_human_genes("GO:0000009"))

    expect_equal(actual.go.09, expected.go.09)
})


test_that("we can grab KOs from GMM matrix", {
    GMM.matrix <- suppressMessages(get_GMM_matrix())
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


test_that("we can run DESeq2 for a KO", {
    mtx.feature.names <- c("k02005_hly_d_family_secretion_protein",
                           "k02007_cobalt_nickel_transport_system_permease_protein",
                           "k00248_butyrate_gene_1",
                           "K00634_butyrate_gene_2",
                           "K00929_butyrate_gene_3")

    # simulate some host data
    relevant.host.genes <- c("PIGV","ALG12")
    nonrelevant.host.genes <- c("goofballgene", "sillygoosegene")
    all.host.genes <- c(relevant.host.genes, nonrelevant.host.genes)
    # certain seeds were causing dispersion issues, so just manually setting that here
    set.seed(2)
    host.gene.counts <- data.frame(vapply(X = all.host.genes,
                                    FUN = function(x) {
                                                        rnbinom(n = 1000,
                                                                size = 400000,
                                                                prob = runif(1))
                                                        },
                                    FUN.VALUE = numeric(1000)
                                      )
                               )

    # simulate some mtx data
    microbial.gene <- "k00248_butyrate_gene_1"
    microbial.gene.counts <- data.frame(vapply(X=mtx.feature.names,
                                            FUN=function(x){
                                                rnorm(n=1000, mean=100)},
                                            FUN.VALUE=numeric(1000))
                                     )

    results <- suppressMessages(
               suppressWarnings(
                   go_targeted_diffex(targeted.genes=relevant.host.genes,
                                  host.genes=t(host.gene.counts),
                                  microbial.gene=microbial.gene,
                                  microbial.genes=microbial.gene.counts,
                                  verbose=FALSE)
                   )
                )

    # There should be 2 rows with no significant results
    expect_equal(nrow(results), 2)
    expected.colnames <- c("baseMean", "log2FoldChange", "lfcSE",
                           "stat", "pvalue", "padj")
    expect_equal(colnames(results), expected.colnames)
    expect_true(sum(results$padj < 0.05) == 0)
})


test_that("we can run DESeq2 for KOs within a GMM", {
    GMM.matrix <- suppressMessages(get_GMM_matrix())
    mtx.feature.names <- c("k02005_hly_d_family_secretion_protein",
                           "k02007_cobalt_nickel_transport_system_permease_protein",
                           "k00248_butyrate_gene_1",
                           "K00634_butyrate_gene_2",
                           "K00929_butyrate_gene_3")

    # simulate some host data
    relevant.host.genes <- c("PIGV","ALG12")
    nonrelevant.host.genes <- c("goofballgene", "sillygoosegene")
    all.host.genes <- c(relevant.host.genes, nonrelevant.host.genes)
    set.seed(1)
    host.gene.counts <- data.frame(vapply(X = all.host.genes,
                                    FUN = function(x) {
                                                        rnbinom(n = 1000,
                                                                size = 400000,
                                                                prob = runif(1))
                                                        },
                                    FUN.VALUE = numeric(1000)
                                      )
                               )
    # simulate some mtx data
    microbial.gene <- "k00248_butyrate_gene_1"
    microbial.gene.counts <- data.frame(vapply(X=mtx.feature.names,
                                            FUN=function(x){
                                                rnorm(n=1000, mean=100)},
                                            FUN.VALUE=numeric(1000))
                                     )

    mtx.features.to.test <- features_from_gmm_df("butyrate production I",
                                                 GMM.matrix,
                                                 mtx.feature.names=mtx.feature.names)
    results <- suppressMessages(
        suppressWarnings(
            GO_targeted_for_each_KO_within_GMM(go.terms=c("GO:0000009"),
                               host.genes=t(host.gene.counts),
                               mtx=microbial.gene.counts,
                               mtx.features=mtx.features.to.test)
        )
    )

    # There should be 2 rows with no significant results
    expect_equal(nrow(results), 6)
    expected.colnames <- c("baseMean", "log2FoldChange", "lfcSE",
                            "stat", "pvalue", "padj", "term")
    expect_equal(colnames(results), expected.colnames)
    expect_true(sum(results$padj < 0.05) == 0)
})
