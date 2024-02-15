
library(dplyr)
library(tidyr)
library(rootSolve)

rm(list = ls())

# wczytanie danych

houses <- read.csv("path")
View(houses)
summary(houses)

# prepare variables
houses$pricelog = log(houses$price)
houses$sqlog = log(houses$sq)
houses$age = (2024 - houses$year)

y = houses$pricelog
x1 = houses$rooms
x2 = houses$sqlog
x3 = houses$floor
x4 = houses$airport
x5 = houses$age

# results of estimation 
uklad.rownan = function ( b) {
  b0 = b[1]
  b1 = b[2]
  b2 = b[3]
  b3 = b[4]
  b4 = b[5]
  b5 = b[6]
  r1 = sum (y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 )
  r2 = sum (x1 *(y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 ))
  r3 = sum (x2 *(y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 ))
  r4 = sum (x3 *(y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 ))
  r5 = sum (x4 *(y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 ))
  r6 = sum (x5 *(y - b0-b1*x1 -b2*x2 - b3*x3 - b4*x4 - b5*x5 ))
  return (c( r1, r2, r3, r4, r5, r6))
}
wynik = multiroot(f = uklad.rownan , start = rep(0, times =6))
wynik$root
b0.hat = wynik$root[1]
b1.hat = wynik$root[2]
b2.hat = wynik$root[3]
b3.hat = wynik$root[4]
b4.hat = wynik$root[5]
b5.hat = wynik$root[6]

# Var-Cov matrix
m = cbind (( y - b0.hat -b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5  ) ,
           (x1 *(y -b0.hat - b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5 )),
           (x2 *(y -b0.hat - b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5 )),
           (x3 *(y -b0.hat - b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5 )),
           (x4 *(y -b0.hat - b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5 )),
           (x5 *(y -b0.hat - b1.hat *x1 -b2.hat*x2 - b3.hat*x3 - b4.hat*x4 - b5.hat*x5 ))
           )

N = length(y)

D = 1/N * cbind(c(-N, -sum(x1), -sum(x2), -sum(x3), -sum(x4), -sum(x5)),
                c(sum(-x1), sum(-x1^2), sum(-x2*x1), sum(-x3*x1), sum(-x4*x1),sum(-x5*x1)),
                c(sum(-x2), -sum(x1*x2), sum(-x2^2),sum(-x3*x2), sum(-x4*x2),sum(-x5*x2)),
                c(sum(-x3), -sum(x1*x3), sum(x2*x3), sum(-x3^2), sum(-x4*x3),sum(-x5*x3)),
                c(sum(-x4), -sum(x1*x4),sum(x2*x4), sum(-x3*x4), sum(-x4^2),sum(-x5*x4)),
                c(sum(-x5), -sum(x1*x5),sum(x2*x5), sum(-x3*x5), sum(-x4*x5),sum(-x5^2)))

J = t(m)%*%m
Omega = solve(D)%*%J%*%t(solve(D))/N^2
std.err = sqrt(diag(Omega))
std.err
  

# comparison with mnk
model = lm(y ~ x1 + x2 + x3 + x4 + x5 )
( sum.mod = summary(model)$coefficients )
# Porownanie w jednej tablicy
porownanie = data.frame( wynik$root , std.err , sum.mod[,1], sum.mod[ ,2])
colnames (porownanie) = c("MM.b", "MM.std.err ", "MNK.b", " MNK.std.err ")

rownames(porownanie) = c("intercept", "rooms", "sqlog", "floor",
                      "airport", "age")

print ( porownanie )

# test of variables significant
Bety <- wynik$root
Bety
ztest0 <- (Bety[1] - 0)/std.err[1]
p_value0 <- 2 * pt(abs(ztest0), df = N - 1, lower.tail = FALSE)
ztest1 <- (Bety[2] - 0)/std.err[2]
p_value1 <- 2 * pt(abs(ztest1), df = N - 1, lower.tail = FALSE)
ztest2 <- (Bety[3] - 0)/std.err[3]
p_value2 <- 2 * pt(abs(ztest2), df = N - 1, lower.tail = FALSE)
ztest3 <- (Bety[4] - 0)/std.err[4]
p_value3 <- 2 * pt(abs(ztest3), df = N - 1, lower.tail = FALSE)
ztest4 <- (Bety[5] - 0)/std.err[5]
p_value4 <- 2 * pt(abs(ztest4), df = N - 1, lower.tail = FALSE)
ztest5 <- (Bety[6] - 0)/std.err[6]
p_value5 <- 2 * pt(abs(ztest5), df = N - 1, lower.tail = FALSE)


z.stats <- c()
for (i in 0:5) {
  z.stats <- c(z.stats, get(paste0("ztest",i )))
} 
z.stats


p_values <- c()
for (i in 0:5) {
  p_values <- c(p_values, get(paste0("p_value",i )))
} 
p_values

# results
wyniki = data.frame( wynik$root , std.err , z.stats, round(p_values, digits = 5))
colnames (wyniki) = c("Oszacowanie", "Błąd stand. ", "Statystyka z", "p-value")
rownames(wyniki) = c("intercept", "rooms", "sqlog", "floor",
                         "airport", "age")

print ( wyniki )

# save
write.table(wyniki, "wyniki.csv", sep = "\t", row.names = TRUE, fileEncoding = "UTF-8")

# additional test 
# h0: b3.hat = 1.0, h1: ~h0
ztest6 <- (Bety[3] - 1)/std.err[6]
p_value6 <- 2 * pt(abs(ztest5), df = N - 1, lower.tail = FALSE)
p_value6



