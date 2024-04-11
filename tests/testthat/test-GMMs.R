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

    file.remove("GMMs.txt")
})

