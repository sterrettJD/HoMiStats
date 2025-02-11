test_that("Download works", {
    pull_GMMs_file("GMMs.txt")
    # file exists
    expect_true(file.exists("GMMs.txt"))

    # file isn't empty
    expect_equal(readLines(file("GMMs.txt","r"), n=1),
    "MF0001	arabinoxylan degradation")

    # func warns if the file is already there
    expect_warning(pull_GMMs_file("GMMs.txt"),
                   "File already exists.")
    if (Sys.info()[["sysname"]] != "Windows") {
        # file should be cleaned up
        # There are some weird file locks on Windows
        # when running tests in parallel
        # So suppress the warnings here
        file.remove("GMMs.txt")
    } else {
        suppressWarnings(file.remove("GMMs.txt"))
    }
})


test_that("parsing works", {


    if (Sys.info()[["sysname"]] != "Windows") {
        gmms.df <- suppressMessages(
            get_GMM_matrix("GMMs.txt", cleanup=TRUE)
        )
        # file should be cleaned up
        # There are some weird file locks on Windows
        # when running tests in parallel
        # So skip this on windows due to the locks

        expect_false(file.exists("GMMs.txt"))
    } else {
        gmms.df <- suppressMessages(
            get_GMM_matrix("GMMs.txt", cleanup=FALSE)
        )
    }
    # There should be 3 columns with these names
    expect_equal(colnames(gmms.df),
                 c("KEGG","Module","Module ID"))

    # There should be 581 rows
    expect_equal(nrow(gmms.df), 581)

    # Check that it also works with a predownloaded file
    pull_GMMs_file("predownloaded_GMMs.txt")
    predownloaded.gmms.df <- suppressMessages(
                                get_GMM_matrix("predownloaded_GMMs.txt")
                             )

    expect_equal(predownloaded.gmms.df, gmms.df)
})
