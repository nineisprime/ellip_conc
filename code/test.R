
library(Matrix)
source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")
source("helper.R")

options(digits=20)

n = 300
p = 1e12

##X = rgamma(n, shape=k/2, scale=1/sqrt(k/2))
X = rgamma(n, shape=p, scale=1)
Y = sqrt(X)
Y = sort(Y)


phi = fit_1d_density(Y, p, M=50000)

hat_density = exp(phi +
    (p-1)/sqrt(p) * (Y - sqrt(p)) -
    (p-1)/p * (Y - sqrt(p))^2 +
    (p-1)/p^(3/2) * (Y - sqrt(p))^3 -
    (p-1)/p^2 * (Y - sqrt(p))^4 +
    (p-1)/p^(5/2) * (Y - sqrt(p))^5 )

plot(Y, hat_density, type="l")


##true_density = dgamma(X, shape=p, scale=1)
true_density = function(Y){
    tmp = (2*p - 1)*log(Y) - Y^2 + log(2) -
        (p-1)*log(p-1) + (p-1) - (1/2)*log(2*pi*(p-1))

    return(exp(tmp))
}


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
