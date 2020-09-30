library(testthat)
library(data.table)
library(FARMregress)

test_check("FARMregress")

x <- fread("x.csv")
y <- read.csv("y.csv")$x
out <- farmRegress(x,y,type="regression")
