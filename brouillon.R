library(tidyverse)
hist(rbinom(1000, 10000, c(0.9)))
x = seq(1:380)
y = 1000+450*(x)^(1/2)
M = tibble(x,y)
n =rep(381:1000, each=1)
l = c(rep(9800, 10), rep(10000, c(length(n)-10)))
G = tibble(n,l) %>% rename(x = n, y = l)
H = rbind(M,G)

plot(H)

TU = seq(1:17)*50
YI = c(1000,3000, 5000, 6000, 6800, 7500, 7900, 8300, 8500, 9000, 9155, 9360, 9700, 9900, 10000, 10000, 10000)
TUYI = tibble(TU, YI)
plot(TUYI)

mean(YI)
mean(TU)
mean(YI)/ mean(TU)
