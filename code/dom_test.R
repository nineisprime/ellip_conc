## Experiments to explore the stochastic dominance question
##
## If F >= G, is F* >= G* ?

library(logcondens)

X = seq(from=0, to=1, by=0.1)
res = logConDens(X, smoothed=FALSE)

## Y dominates X
Y = sort(X)
eps = 0.2
Y[length(Y)] = Y[length(Y)] + eps
res_new = logConDens(Y, smoothed=FALSE)

## F_Y* (red line) should lie below F_X* (black line)
plot(sort(X), res$Fhat, type='l')
lines(sort(Y), res_new$Fhat, col='red')


