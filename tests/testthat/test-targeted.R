test_that("we can grab gene names from go terms", {

    expected.go.09 <- c("PIGV","ALG12")
    actual.go.09 <- suppressMessages(get_go_term_human_genes("GO:0000009"))

    expect_equal(actual_go_09, expected_go_09)
})
