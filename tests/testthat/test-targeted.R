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
    host.gene.counts <- data.frame(sapply(X=all.host.genes,
                                         FUN=function(x){
                                             rnbinom(n=1000, size=400000, prob=runif(1))
                                             }
                                         )
                                  )
    # simulate some mtx data
    microbial.gene <- "k00248_butyrate_gene_1"
    microbial.gene.counts <- data.frame(sapply(X=mtx.feature.names,
                                         FUN=function(x){
                                             rbinom(n=1000, size=1000, prob=0.05)
                                             }
                                         )
                                  )

    results <- suppressMessages(
               suppressWarnings(
                   go_targeted_diffex(targeted.genes=relevant.host.genes,
                                  host.genes=t(host.gene.counts),
                                  microbial.gene=microbial.gene,
                                  microbial.genes=microbial.gene.counts,
                                  verbose=F)
                   )
                )

    # There should be 2 rows with no significant results
    expect_equal(nrow(results), 2)
    expected.colnames <- c("baseMean", "log2FoldChange", "lfcSE",
                           "stat", "pvalue", "padj")
    expect_equal(colnames(results), expected.colnames)
    expect_true(sum(results$padj < 0.05) == 0)
})
