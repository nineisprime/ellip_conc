

library(Matrix)

source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")
source("hellinger.R")

if (!exists("myname"))
    myname = "A"

continue_run = FALSE
options(digit=20)
ntrials = 50
#p_ls = c(5, 1e3, 1e6, 1e9, 1e12, 1e15)
p_ls = c(1e15)

res = matrix(0, ntrials, length(p_ls))

n = 250

if (continue_run){
    load("test_tmp.R")
    source("fit_1d_density.R")
    source("active_set_newton.R")
    source("active_set_newton_largep.R")
    source("numerical_integration_helper.R")
    source("hellinger.R")
}

for (ip in 1:length(p_ls)){
    for (it in 1:ntrials) {
        print(c(ip, it))

        if (continue_run && res[it, ip] != 0)
            next
        
        p = p_ls[ip]

        n = 250
        X = rgamma(n, shape=p, scale=1)
        Y = sqrt(X)
        Y = sort(Y)

        phi = fit_1d_density(Y, p)

        true_density = x_2p_minus_1_e_xsquared(Y, p)

        hell = compute_hellinger(Y, phi, p, M=400000, true_density, offset=sqrt(p),
            boundary=2)

        res[it, ip] = hell

        save.image(paste0("test_", myname, ".RData"))
    }
}
