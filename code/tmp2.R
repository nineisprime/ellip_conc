


nI = length(I)

YI_sorted = Y_sorted[I]

## precompute grad_ypart
grad_ypart = rep(0, nI)
grad_ypart[1] = I[1] 

for (i in 1:(nI-1)){
    if (nI == 1) break
    
    ixs = (I[i]+1):I[i+1]
    
    y_lvec = rep(YI_sorted[i], length(ixs))
    y_rvec = rep(YI_sorted[i+1], length(ixs))
    ygap = YI_sorted[i+1] - YI_sorted[i]
    
    grad_ypart[i] = grad_ypart[i] + sum( (y_rvec - Y_sorted[ixs])/ygap )
    grad_ypart[i+1] = grad_ypart[i+1] + sum( (Y_sorted[ixs] - y_lvec)/ygap )
    
}
grad_ypart = grad_ypart/n


if (nI > 1){
    num_int_res = numerical_integration_helper(YI_sorted, M, eps)
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


phi_lvec = rep2(ind_vec2, phi[1:(nI-1)])
phi_rvec = rep2(ind_vec2, phi[2:nI])

fnpts =  exp(
    ((y_rvec-evalpts)*phi_lvec +
     (evalpts-y_lvec)*phi_rvec )/gap_vec +
    ((p-1)/sqrt(p)) * (evalpts - sqrt(p)) -
    ((p-1)/p) * (evalpts - sqrt(p))^2 * 1/2 +
    ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^3 * 1/3 -
    ((p-1)/p^2) * (evalpts - sqrt(p))^4 * 1/4 +
    ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^5 * 1/5 -
    ((p-1)/p^3) * (evalpts - sqrt(p))^6 * 1/6 +
    ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^7 * 1/7)


fn_int = sum(fnpts * kepler_global_wts) +
    exp(phi[1] +
        sqrt(p)*(min(YI_sorted) - sqrt(p)) -
        (min(YI_sorted) - sqrt(p))^2 * 1/2 +
        (1/sqrt(p)) * (min(YI_sorted) - sqrt(p))^3 * 1/3 -
        (1/p) * (min(YI_sorted) - sqrt(p))^4 * 1/4 +
        (1/p^(3/2)) * (min(YI_sorted) - sqrt(p))^5 * 1/5 -
        (1/p^2) * (min(YI_sorted) - sqrt(p))^6 * 1/6 +
        (1/p^(5/2)) * (min(YI_sorted) - sqrt(p))^7 * 1/7 -
        (1/p^3) * (min(YI_sorted) - sqrt(p))^8 -
        log(sqrt(p)))


new_obj = sum(grad_ypart*phi) - fn_int
