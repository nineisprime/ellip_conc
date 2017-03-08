

library(Matrix)

source("fit_1d_density.R")
source("active_set_newton.R")
source("active_set_newton_largep.R")
source("numerical_integration_helper.R")


continue_run = TRUE
options(digit=20)
ntrials = 50
p_ls = c(50, 1e3, 1e6, 1e9)

res = matrix(0, ntrials, length(p_ls))


if (continue_run){
    load("test_tmp.R")
    source("fit_1d_density.R")
    source("active_set_newton.R")
    source("active_set_newton_largep.R")
    source("numerical_integration_helper.R")

}

for (ip in 1:length(p_ls)){
    for (it in 1:ntrials) {
        print(c(ip, it))

        if (continue_run && res[it, ip] != 0)
            next
        
        p = p_ls[ip]

        n = 300
        X = rgamma(n, shape=p, scale=1)
        Y = sqrt(X)
        Y = sort(Y)

        phi = fit_1d_density(Y, p)

        true_density = function(Y){
            tmp = (2*p - 1)*log(Y) - Y^2 + log(2) -
                (p-1)*log(p-1) + (p-1) - (1/2)*log(2*pi*(p-1))

            return(exp(tmp))
        }

        hell = compute_hellinger(Y, phi, p, M=400000, true_density, offset=sqrt(p),
            boundary=2)

        res[it, ip] = hell

        save.image("test_tmp.R")
    }
}
