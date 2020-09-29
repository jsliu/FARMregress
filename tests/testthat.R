library(testthat)
library(FARMregress)

test_check("FARMregress")

x <- fread("data/x.csv")
y <- fread("data/y.csv")[,y]
out <- farmRegress(x,y,type="classification")
