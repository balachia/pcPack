context('Exponential utility functions')

test_that('Expected utility ratio is correct', {
    expect_equal(euratio.normal(0), 1/2)
    expect_equal(euratio.normal(log(10)), 10/11)
    expect_equal(euratio.normal(0, a=2, b=3), 9/5)

    expect_equal(euratio.logged(0), -log(2))
    expect_equal(euratio.logged(log(10)), log(10) - log(11))
    expect_equal(euratio.logged(0, a=2, b=3), log(9) - log(5))

    expect_equal(euratio(0), -log(2))
    expect_equal(euratio(0, log=TRUE), -log(2))
    expect_equal(euratio(0, log=FALSE), 1/2)
})

test_that('exparg is built correctly', {
    expect_equal(exparg.open(1, 0), 3/2)
    expect_equal(exparg.open(2, 1, sigma=2, b=2), 15)

    expect_equal(exparg.brid(1, 1, 5, 1, 2), 11/24)
    expect_equal(exparg.brid(1, 1, 5, 1, 2, b=2, sigma=2), 37/6)
})
