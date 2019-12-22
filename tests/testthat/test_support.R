context("Testing Support Functions")

test_that(desc="Testing Quantile Binning Function", {
  expect_equal(quant_bin_1d(1:30, nbin=3), rep(1:3,each=10))
})

test_that(desc="Testing Conversion to R-tree list", {
  expect_equal({iq_def <- iqbin(data=data.frame(x1=1:16,x2=1:16),
                               bin_cols=c("x1","x2"),
                               nbins=c(2,2), output="definition")
               rtree <- make_bin_list(bin_bounds=iq_def$bin_bounds,nbins=iq_def$nbins)
               rtree[[1]][[1]][[1]]},
               c(1.0, 4.5, 8.0))
})

test_that(desc="Testing R-tree index checker", {
  expect_equal(index_collapser(nbins=c(2,2,2),indeces=c(1,1,1)),1)
  expect_equal(index_collapser(nbins=c(2,2,2),indeces=c(1,1,2)),2)
  expect_equal(index_collapser(nbins=c(2,2,2),indeces=c(2,2,1)),7)
})

test_that(desc="Testing R-tree index lookup for new point", {
  expect_equal({
    iq_def <- iqbin(data=data.frame(x1=1:16,x2=1:16),
                    bin_cols=c("x1","x2"),
                    nbins=c(2,2), output="definition")
    bin_index_finder_nest(x=c(10,10), bin_def=iq_def, strict=TRUE)
  },3)
})

test_that(desc="Testing stack matrix build", {
  expect_equal(make_stack_matrix(3,2),{
    matrix(c(1,1,0,0,0,0,
             0,0,1,1,0,0,
             0,0,0,0,1,1),ncol=3)
  })
})

test_that(desc="Testing CV cohort balance", {
  expect_equal({
    tab1 <- table(make_cv_cohorts(dat=data.frame(1:17),5))
    attributes(tab1) <- NULL
    tab1
    },{
    c(4,4,3,3,3)
  })
})

test_that(desc="Testing Majority Vote Counter", {
  expect_equal(majority_vote(c("a","a","a","b","b","c")), "a")
})

test_that(desc="Testing Tuning Binlist Options", {
  expect_equal(make_nbins_list(c(2,3),3),
               list(c(2,2,2),c(3,2,2),c(3,3,2),c(3,3,3)))
})




