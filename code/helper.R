

## various helper functions

## Approximating
## 
##      
## exp( ((Y_{i+1} - r) phi_l + (r - Y_i) phi_r)/gap_i + 
##               (p-1) * log( 1 + (r - mu)/mu )
##
## = exp( phi(r) ) (r/mu)^(p-1)
##
##
hx_eval <- function(phi_lvec, phi_rvec, evalpts, y_lvec, y_rvec, gap_vec, mu, p) {

    if (p < 100){
        out =  exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
        (p-1) * log(1 + (evalpts - mu)/mu))
    } else {
        out = exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
        ((p-1)/mu)*(
            (evalpts - mu) - (1/2)*(evalpts - mu)^2 / mu +
            (1/3)*(evalpts - mu)^3 / mu^2 -
            (1/4)*(evalpts - mu)^4 / mu^3 +
            (1/5)*(evalpts - mu)^5 / mu^4 -
            (1/6)*(evalpts - mu)^6 / mu^5))
    }
    
    ## out = exp(
    ##     ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
    ##     sqrt(p)*(
    ##         ((p-1)/p) * (evalpts - sqrt(p)) -
    ##         ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^2 * 1/2 +
    ##         ((p-1)/p^(4/2)) * (evalpts - sqrt(p))^3 * 1/3 -
    ##         ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^4 * 1/4 +
    ##         ((p-1)/p^(6/2)) * (evalpts - sqrt(p))^5 * 1/5 -
    ##         ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^6 * 1/6 +
    ##         ((p-1)/p^(8/2)) * (evalpts - sqrt(p))^7 * 1/7))
    
    return(out)
}


hx1_int <- function(phi1, minY, mu, p){


    if (p < 100){
        out = exp(phi1 + p*log(1 + (minY - mu)/mu) + log(mu/p))
    } else {
        out = exp(phi1 + (p/mu)*
                  ((minY - mu) -
                   (1/2)*(minY - mu)^2 / mu +
                   (1/3)*(minY - mu)^3 / mu^2 -
                   (1/4)*(minY - mu)^4 / mu^3 +
                   (1/5)*(minY - mu)^5 / mu^4 -
                   (1/6)*(minY - mu)^6 / mu^5)
                  + log(mu/p))
    }
    
    ## out = exp(
    ##     phi1 + 
    ##     sqrt(p) * (
    ##         (minY - sqrt(p)) -
    ##         p^(-1/2) * (minY - sqrt(p))^2 * 1/2 +
    ##         p^(-2/2) * (minY - sqrt(p))^3 * 1/3 -
    ##         p^(-3/2) * (minY - sqrt(p))^4 * 1/4 +
    ##         p^(-4/2) * (minY - sqrt(p))^5 * 1/5 -
    ##         p^(-5/2) * (minY - sqrt(p))^6 * 1/6 +
    ##         p^(-6/2) * (minY - sqrt(p))^7 * 1/7 ) -
    ##     log(sqrt(p)) ) 
    return(out)
}
    
