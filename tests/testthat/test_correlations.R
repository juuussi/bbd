library(bbd)

context("Test correlations")

test_that("result is a data.frame", {
  expect_equal(class(correlations(data=iris, x.variables = colnames(iris)[1:3], parallel=FALSE)), "data.frame")
})

test_that("result has correct dimensions", {
  expect_equal(ncol(correlations(data=iris, x.variables = colnames(iris)[1:3], parallel=FALSE)), 8)
  expect_equal(nrow(correlations(data=iris, x.variables = colnames(iris)[1:3], parallel=FALSE)), 3)
})

test_that("NULL value input", {
  expect_error(correlations(NULL))
})
