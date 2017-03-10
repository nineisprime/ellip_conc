

## INPUT:
##   Y -- of length n
##  phi -- of same length as Y
##   p -- positive integer
##  M -- number of grid points to use for numerical integration
##  true_density -- function that takes in a vector of numbers as input and outputs a
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

compute_hellinger <- function(Y, phi, p, M=50000, true_density, offset, boundary = 1) {

    Y = sort(Y)
    n = length(Y)
    eps = 1e-6

    ## Y_aug contains Y and the limits of integration
    Y_aug = c(Y[1] - boundary, Y, Y[n] + boundary)
    phi_aug = c(phi[1], phi, -Inf)

    
    num_int_res = numerical_integration_helper(Y_aug, M, eps)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    kepler_global_wts = num_int_res$kepler_global_wts
    ind_vec2 = num_int_res$ind_vec2


    phi_lvec = rep2(ind_vec2, phi_aug[1:(n+1)])
    phi_rvec = rep2(ind_vec2, phi_aug[2:(n+2)])
    
    
    if (p > 100){
        
        fnhat_pts = exp(
            ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
            ((p-1)/sqrt(p)) * (evalpts - offset) -
            ((p-1)/p) * (evalpts - offset)^2 * 1/2 +
            ((p-1)/p^(3/2)) * (evalpts - offset)^3 * 1/3 -
            ((p-1)/p^2) * (evalpts - offset)^4 * 1/4 +
            ((p-1)/p^(5/2)) * (evalpts - offset)^5 * 1/5)
        
    } else {
        
        fnhat_pts =  exp(
            ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
            (p-1)*log(evalpts/offset) )
    }

    f_pts = true_density(evalpts)

    out = sum(kepler_global_wts * (sqrt(fnhat_pts) - sqrt(f_pts))^2)

    return(out)
}
        
        
