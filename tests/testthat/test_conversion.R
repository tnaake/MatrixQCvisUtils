## function biocrates
test_that("biocrates", {
    expect_is(biocrates, "function")
    expect_error(biocrates(file = ""), "File does not exist")

})

## function metaboscape
test_that("metaboscape", {
    expect_is(metaboscape, "function")
    expect_error(metaboscape(file = ""), "File does not exist")
})

## function maxquant
test_that("maxquant", {
    expect_is(maxquant, "function")
    expect_error(maxquant(file = "", intensity = "foo"), "should be one of")
    suppressWarnings(expect_error(maxquant(file = "", intensity = "iBAQ", 
        type = "txt"), "no lines available"))
    expect_error(maxquant(file = "", intensity = "iBAQ", type = "xlsx"), 
        "File does not exist")
    suppressWarnings(expect_error(maxquant(file = "", intensity = "LFQ", 
        type = "txt"), "no lines available"))
    expect_error(maxquant(file = "", intensity = "LFQ", type = "xlsx"), 
        "File does not exist")
    expect_error(maxquant(file = "", intensity = "LFQ", type = "foo"),
        "should be one of")
})

## function diann
test_that("diann", {
    expect_is(diann, "function")
    expect_error(diann(file = ""), "no lines available")
})

## function spectronaut
test_that("spectronaut", {
    expect_is(spectronaut, "function")
    expect_error(spectronaut(file = ""), "File does not exist")
})