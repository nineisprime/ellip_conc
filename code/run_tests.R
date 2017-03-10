

library(Matrix)

source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")
source("hellinger.R")
source("density_zoo.R")

if (!exists("myname"))
    myname = "A"

continue_run = FALSE
options(digit=20)
ntrials = 50
p_ls = c(1, 1e2, 1e4, 1e6, 1e8, 1e10, 1e12)
#p_ls = c(1e15)

if (!continue_run){
    res = matrix(0, ntrials, length(p_ls))
}

n = 300

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

        phi = fit_1d_density(Y, p, M=200000)

        true_density = function(Y) {
            x_2p_minus_1_e_xsquared(Y, p)
        }

        hell = compute_hellinger_ellip(Y, phi, p, M=400000, true_density, offset=sqrt(p),
            boundary=2)

        res[it, ip] = hell

        save.image(paste0("test_", myname, ".RData"))
    }
}
