# Entry point: R CMD check runs this file, which runs every
# tests/testthat/test-*.R file against the installed package.

library(testthat)
library(boundedHill)

test_check("boundedHill")
