
library(Matrix)
source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")

options(digits=20)

n = 300
p = 50

##X = rgamma(n, shape=k/2, scale=1/sqrt(k/2))
X = rgamma(n, shape=p, scale=1)
Y = sqrt(X)
Y = sort(Y)

##A = c(2:10, 12:(n-1))
##phi = active_set_newton_largep(Y, A, p, diagnostic=TRUE)

I = c(1, 300)
I = sort(I)
#I = c(5, 300)
A = setdiff(1:n, I)

phi = active_set_newton_largep(Y, A, p, diagnostic=TRUE)
phi_old = phi



Y = sort(Y)
n = length(Y)
ygap = Y[2:n] - Y[1:(n-1)]
## create V matrix representing constraints
## V'phi <= 0
diffmat = bandSparse(n-1, n, k=c(0,1),
    diagonals=list( -1/ygap, 1/ygap ))
diffmat2 = bandSparse(n-2, n-1, k=c(0,1),
    diagonals=list(rep(-1,n-2), rep(1, n-2)))
sec_diffmat = diffmat2 %*% diffmat
sec_diffmat = rBind(c(-1/ygap[1], 1/ygap[1], rep(0,n-2)), sec_diffmat, rep(1, n))
V = sec_diffmat
## rows of B represent dual-vectors
## Each V_i is a row of V
## Each B_i is a row of B
##
## B_i'V_i = -1
## B_i'V_j = 0
B = -t(solve(V))
V = V[1:(n-1), ]
##B = B[1:(n-1), ]




num_int_res = numerical_integration_helper(Y, M=500000, eps)
evalpts = num_int_res$evalpts
bdpts = num_int_res$bdpts
y_rvec = num_int_res$y_rvec
y_lvec = num_int_res$y_lvec
gap_vec = num_int_res$gap_vec
ind_vec2 = num_int_res$ind_vec2
kepler_global_wts = num_int_res$kepler_global_wts
kepler_ls = num_int_res$kepler_ls
ixs_ls = num_int_res$ixs_ls

phi = expand_phi(phi, Y, A)

grad = compute_gradient_phi_largep(phi, y_lvec=y_lvec, y_rvec=y_rvec, ind_vec2=ind_vec2,
    gap_vec=gap_vec, kepler_ls=kepler_ls, evalpts=evalpts, ixs_ls=ixs_ls, p=p)
dual_vec = B %*% grad




phi = fit_1d_density(Y, p)

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

hell = compute_hellinger(Y, phi, p, M=300000, true_density, offset=sqrt(p),
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
