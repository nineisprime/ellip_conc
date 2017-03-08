

## script used to compare objective of different optimization arguments
##
## phi and phi2 needs to be expanded first

phi_x = expand_phi(phi_old, Y, A)

## precompute quantities needed for numerical integration
M = 500000
if (nI > 1){
    num_int_res = numerical_integration_helper(Y, M, eps)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    ind_vec2 = num_int_res$ind_vec2
    kepler_global_wts = num_int_res$kepler_global_wts
    kepler_ls = num_int_res$kepler_ls
    ixs_ls = num_int_res$ixs_ls
}


phi_lvec = rep2(ind_vec2, phi_x[1:(n-1)])
phi_rvec = rep2(ind_vec2, phi_x[2:n])

fnpts_cand =  exp(
    ((y_rvec-evalpts)*phi_lvec +
     (evalpts-y_lvec)*phi_rvec )/gap_vec +
    ((p-1)/sqrt(p)) * (evalpts - sqrt(p)) -
    ((p-1)/p) * (evalpts - sqrt(p))^2 * 1/2 +
    ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^3 * 1/3 -
    ((p-1)/p^2) * (evalpts - sqrt(p))^4 * 1/4 +
    ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^5 * 1/5 -
    ((p-1)/p^3) * (evalpts - sqrt(p))^6 * 1/6 +
    ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^7 * 1/7)

fn_int_cand = sum(fnpts_cand * kepler_global_wts) +
    exp(phi_x[1] +
        sqrt(p)*(min(Y) - sqrt(p)) -
        (min(Y) - sqrt(p))^2 * 1/2 +
        (1/sqrt(p)) * (min(Y) - sqrt(p))^3 * 1/3 -
        (1/p) * (min(Y) - sqrt(p))^4 * 1/4 +
        (1/p^(3/2)) * (min(Y) - sqrt(p))^5 * 1/5 -
        (1/p^2) * (min(Y) - sqrt(p))^6 * 1/6 +
        (1/p^(5/2)) * (min(Y) - sqrt(p))^7 * 1/7 - log(sqrt(p)))

obj = (1/n)*sum(phi_x) - fn_int_cand


fn_so_far = exp(phi_x[1] +
        sqrt(p)*(min(Y) - sqrt(p)) -
        (min(Y) - sqrt(p))^2 * 1/2 +
        (1/sqrt(p)) * (min(Y) - sqrt(p))^3 * 1/3 -
        (1/p) * (min(Y) - sqrt(p))^4 * 1/4 +
        (1/p^(3/2)) * (min(Y) - sqrt(p))^5 * 1/5 -
        (1/p^2) * (min(Y) - sqrt(p))^6 * 1/6 +
        (1/p^(5/2)) * (min(Y) - sqrt(p))^7 * 1/7 - log(sqrt(p)))
for (i in 1:4)    
    fn_so_far = fn_so_far +  sum(fnpts_cand[ixs_ls[[i]]] * kepler_ls[[i]])







#################
## Try new Phi
phi2 = phi_x - 0.0001 * B[1, ]

phi2_lvec = rep2(ind_vec2, phi2[1:(n-1)])
phi2_rvec = rep2(ind_vec2, phi2[2:n])

fnpts_cand =  exp(
    ((y_rvec-evalpts)*phi2_lvec +
     (evalpts-y_lvec)*phi2_rvec )/gap_vec +
    ((p-1)/sqrt(p)) * (evalpts - sqrt(p)) -
    ((p-1)/p) * (evalpts - sqrt(p))^2 * 1/2 +
    ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^3 * 1/3 -
    ((p-1)/p^2) * (evalpts - sqrt(p))^4 * 1/4 +
    ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^5 * 1/5 -
    ((p-1)/p^3) * (evalpts - sqrt(p))^6 * 1/6 +
    ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^7 * 1/7)

fn_int_cand = sum(fnpts_cand * kepler_global_wts) +
    exp(phi2[1] +
        sqrt(p)*(min(Y) - sqrt(p)) -
        (min(Y) - sqrt(p))^2 * 1/2 +
        (1/sqrt(p)) * (min(Y) - sqrt(p))^3 * 1/3 -
        (1/p) * (min(Y) - sqrt(p))^4 * 1/4 +
        (1/p^(3/2)) * (min(Y) - sqrt(p))^5 * 1/5 -
        (1/p^2) * (min(Y) - sqrt(p))^6 * 1/6 +
        (1/p^(5/2)) * (min(Y) - sqrt(p))^7 * 1/7 - log(sqrt(p)))

obj2 = (1/n)*sum(phi2) - fn_int_cand

##
tmp = V %*% phi2
which( tmp > 1e-7)

tmp = V %*% phi_x
which( tmp > 1e-7)
