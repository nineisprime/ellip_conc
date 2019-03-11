


## Z must be sorted
## 
compute_Phi <- function(phi, Z, Z2, p){

    ## Sanity check
    stopifnot(length(Z) == length(unique(Z)))

    
    Z2ord = order(Z2)
    Z2 = sort(Z2)
    
    if (min(Z) > 0){
        Z = c(0, Z)
        phi = c(phi[1], phi)
    }

    overflow = which(Z2 > max(Z))
    inbound = which(Z2 <= max(Z))
    Z2 = Z2[inbound]
    
    allpts = sort(c(Z, Z2))
    tmp1 = (allpts %in% Z) & (!duplicated(allpts))
    tmp2 = !tmp1
    ZI = which(tmp1)
    Z2I = which(tmp2)

    n1 = length(Z)

    ind_vec = ZI[2:n1] - ZI[1:(n1-1)]
    ind_vec[1] = ind_vec[1]+1

    z_lvec = rep2(ind_vec, Z[1:(n1-1)])
    z_rvec = rep2(ind_vec, Z[2:n1])

    gap_vec = z_rvec - z_lvec
    
    phi_lvec = rep2(ind_vec, phi[1:(n1-1)])
    phi_rvec = rep2(ind_vec, phi[2:n1])

    
    phi_all = (z_rvec - allpts)/gap_vec * phi_lvec + (allpts - z_lvec)/gap_vec * phi_rvec

    res = c(phi_all[Z2I], rep(-Inf, length(overflow)))
    res2 = rep(0, length(Z2))
    res2[Z2ord] = res
    
    return(res2)
}


##
##

compute_Phi_lc <- function(phi, Z, Z2){

    Z2ord = order(Z2)
    Z2 = sort(Z2)
    

    overflow2 = which(Z2 > max(Z))
    overflow1 = which(Z2 < min(Z))
    inbound = c(which(Z2 >= min(Z) & Z2 <= max(Z)))
    Z2 = Z2[inbound]
    
    
    allpts = sort(c(Z, Z2))
    ZI = which(allpts %in% Z)
    Z2I = which(allpts %in% Z2)

    n1 = length(Z)

    ind_vec = ZI[2:n1] - ZI[1:(n1-1)]
    ind_vec[1] = ind_vec[1]+1

    z_lvec = rep2(ind_vec, Z[1:(n1-1)])
    z_rvec = rep2(ind_vec, Z[2:n1])

    gap_vec = z_rvec - z_lvec
    
    phi_lvec = rep2(ind_vec, phi[1:(n1-1)])
    phi_rvec = rep2(ind_vec, phi[2:n1])

    
    phi_all = (z_rvec - allpts)/gap_vec * phi_lvec + (allpts - z_lvec)/gap_vec * phi_rvec

    res = c(rep(-Inf, length(overflow1)), phi_all[Z2I], rep(-Inf, length(overflow2)))
    res2 = rep(0, length(Z2))
    res2[Z2ord] = res
    
    return(res2)
}



## INPUT:
##   Y -- of length n, of samples
##  phi -- of same length as Y, output of elliptical log-concave MLE 
##   p -- positive integer
##  M -- number of grid points to use for numerical integration
##  true_density -- function that takes in a vector of numbers and "p" as input and outputs a
##                  a vector of density values
##
##  offset -- will work with (Y - offset) instead of Y
##  boundary -- to compute hellinger from [minY - boundary, maxY + boundary]
## 
##
##
## OUTPUT:
##   dist -- real number indicating hellinger distance
##

compute_hellinger_sslc <- function(Y, phi, p, M=400000, true_density, offset=1, boundary = 1) {

    Y = sort(Y)
    n = length(Y)

    ## Y_aug contains Y and the limits of integration
    Y_aug = c(max(Y[1] - boundary, 0), Y, Y[n] + boundary)
    phi_aug = c(phi[1], phi, -Inf)

    
    num_int_res = numerical_integration_helper(Y_aug, M)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    kepler_global_wts = num_int_res$kepler_global_wts
    ind_vec2 = num_int_res$ind_vec2


    phi_lvec = rep2(ind_vec2, phi_aug[1:(n+1)])
    phi_rvec = rep2(ind_vec2, phi_aug[2:(n+2)])
      
    fnhat_pts =  exp(
       ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
       (p-1)*log(evalpts) )
    

    f_pts = true_density(evalpts, p)

    out = sum(kepler_global_wts * (sqrt(fnhat_pts) - sqrt(f_pts))^2)

    return(out)
}


compute_KL_sslc <- function(Y, phi, p, M=400000, logh0, offset=1, boundary = 1) {

    Y = sort(Y)
    n = length(Y)

    ## Y_aug contains Y and the limits of integration
    Y_aug = c(max(Y[1] - boundary, 1e-6), Y)
    phi_aug = c(phi[1], phi)

    
    num_int_res = numerical_integration_helper(Y_aug, M)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    kepler_global_wts = num_int_res$kepler_global_wts
    ind_vec2 = num_int_res$ind_vec2


    phi_lvec = rep2(ind_vec2, phi_aug[1:n])
    phi_rvec = rep2(ind_vec2, phi_aug[2:(n+1)])
      
    fnhat_pts =  exp(
       ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
       (p-1)*log(evalpts) )

    logfnhat_pts = ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
       (p-1)*log(evalpts)
    

    logf_pts = logh0(evalpts, p)

    out = sum(kepler_global_wts * fnhat_pts * (logfnhat_pts - logf_pts))

    return(out)
}





        
## INPUT
##   phi -- value of estimated log-concave MLE on each of the Y
##     Y -- same length as "phi"
##   
##   true_density -- a function indicating what the true density is

compute_hellinger_lc <- function(Y, p, phi, M=400000, true_density, boundary = 1){

    Y = sort(Y)
    n = length(Y)

    ## Y_aug contains Y and the limits of integration
    Y_aug = c(max(Y[1] - boundary, 0), Y, Y[n] + boundary)
    phi_aug = c(-1e100, phi, -Inf)

    
    num_int_res = numerical_integration_helper(Y_aug, M)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    kepler_global_wts = num_int_res$kepler_global_wts
    ind_vec2 = num_int_res$ind_vec2


    phi_lvec = rep2(ind_vec2, phi_aug[1:(n+1)])
    phi_rvec = rep2(ind_vec2, phi_aug[2:(n+2)])
    
    
    fnhat_pts = exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec )


    f_pts = true_density(evalpts, p)

    out = sum(kepler_global_wts * (sqrt(fnhat_pts) - sqrt(f_pts))^2)

    return(out)
}




compute_KL_lc <- function(Y, p, phi, M=400000, logh0, boundary = 1){

    Y = sort(Y)
    n = length(Y)

    ## Y_aug contains Y and the limits of integration
    Y_aug = Y
    phi_aug = phi

    
    num_int_res = numerical_integration_helper(Y_aug, M)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    kepler_global_wts = num_int_res$kepler_global_wts
    ind_vec2 = num_int_res$ind_vec2


    phi_lvec = rep2(ind_vec2, phi_aug[1:(n-1)])
    phi_rvec = rep2(ind_vec2, phi_aug[2:n])
    
    
    fnhat_pts = exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec )

    logfnhat_pts = ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec

    logf_pts = logh0(evalpts, p)

    out = sum(kepler_global_wts * fnhat_pts * (logfnhat_pts - logf_pts))

    return(out)
}



 
## INPUT
##   xs    -- grid on which density is evaluated
##   hhats -- value of the density on "xs"
##   
##   true_density -- a function indicating what the true density is

compute_hellinger_ker <- function(xs, p, hhats, true_density){

    M = length(xs)
    gaps = xs[2:M] - xs[1:(M-1)]


    h0 = true_density(xs[2:M], p)
    hhats = hhats[2:M]
    hell = sum(gaps*(sqrt(h0) - sqrt(hhats))^2)
    return(hell)
}
    

compute_KL_ker <- function(xs, p, hhats, logh0){

    M = length(xs)
    gaps = xs[2:M] - xs[1:(M-1)]


    logh0_pts = logh0(xs[2:M], p)
    hhats = hhats[2:M]
    out = sum( gaps * hhats * (log(hhats) - logh0_pts))
    return(out)
}
    

compute_hellinger_ker <- function(xs, hhats, true_loglike){

    M = length(xs)
    gaps = xs[2:M] - xs[1:(M-1)]

    hhats = hhats[2:M]
    out = sum( gaps *
               (sqrt(hhats) - sqrt(exp(true_loglike[2:M])))^2)
    return(out)
}
    


