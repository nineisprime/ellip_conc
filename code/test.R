
library(Matrix)
source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")
source("helper.R")

options(digits=20)

n = 6000
p = 10

##X = rgamma(n, shape=k/2, scale=1/sqrt(k/2))
X = rgamma(n, shape=p/2, scale=2) 
Y = X^(1/2)
Y = sort(Y)


##Y = runif(n) + sqrt(p)

phi = fit_1d_density(Y, p, M=80000)


hat_density = exp(phi + (p-1)*log(Y))

plot(Y, hat_density, type="l")

ygap = Y[2:n] - Y[1:(n-1)]
sum(hat_density[2:n] * ygap)

##true_density = dgamma(X, shape=p, scale=1)

true_density = function(Y){
    tmp = (p - 1)*log(Y) - (1/2)*Y^2 - lgamma(p/2) + (p/2 - 1)*log(2)

    return(exp(tmp))
}

true_density_vec = true_density(Y)
sum(true_density_vec[2:n]*ygap)

lines(Y, true_density(Y), lty=3, col="red")

hell = compute_hellinger_ellip(Y, phi, p, M=300000, true_density, offset=sqrt(p),
    boundary=3)




## To simulate Gaussian norm
## k = p
##
## Y = sqrt(X)
##
## true_density = Y^(k-1) * exp(-(1/2)*Y^2) * 2 / (gamma(k/2)* 2^(k/2))


## To simulate Gaussian norm squared
## k = 2*p
##
## Y = X
##
## true_density = Y^(k/2 - 1) * exp(-(1/2)*Y) / (gamma(k/2)*2^(k/2))
