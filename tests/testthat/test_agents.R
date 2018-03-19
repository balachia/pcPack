context('Agents')

test_that('standard.agent Ops work', {
    ag11 <- make_standard_agent()
    ag11a <- make_standard_agent()
    ag12 <- make_standard_agent(a=1, b=2)
    ag22 <- make_standard_agent(a=2, b=2)
    ag_insert <- make_insert_agent(data.table::data.table())

    expect_true(ag11==ag11a)
    expect_false(ag11==ag12)
    expect_false(ag11==ag_insert)
})
