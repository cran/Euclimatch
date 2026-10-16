# Test of climatch_vec

test_that("climatch_vec computes correct format", {
i <- as.data.frame(matrix(runif(n=180, min=1, max=20), nrow=60)) # Fake source climate data
j <- as.data.frame(matrix(runif(n=300, min=10, max=40), nrow=100)) # Fake recipient data
variance <- c(600, 800, 450) # Fake global variance

cv <- climatch_vec(recipient = j, source = i, globvar = variance)

expect_type(cv, "double")
})
