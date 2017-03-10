

## various helper functions

## Approximating
## 
## exp( (Y_{i+1} - r) phi(old)_l + (r - Y_i) phi(old)_r + 
##               (p-1) * log( 1 + s/sqrt(p))  )
##      

hx_eval <- function(phi_lvec, phi_rvec, evalpts, y_lvec, y_rvec, gap_vec, p) {
    out = exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
        sqrt(p)*(
            ((p-1)/p) * (evalpts - sqrt(p)) -
            ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^2 * 1/2 +
            ((p-1)/p^(4/2)) * (evalpts - sqrt(p))^3 * 1/3 -
            ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^4 * 1/4 +
            ((p-1)/p^(6/2)) * (evalpts - sqrt(p))^5 * 1/5 -
            ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^6 * 1/6 +
            ((p-1)/p^(8/2)) * (evalpts - sqrt(p))^7 * 1/7))
    
    return(out)
}


hx1_int <- function(phi1, minY, p){
    out = exp(
        phi1 + 
        sqrt(p) * (
            (minY - sqrt(p)) -
            p^(-1/2) * (minY - sqrt(p))^2 * 1/2 +
            p^(-2/2) * (minY - sqrt(p))^3 * 1/3 -
            p^(-3/2) * (minY - sqrt(p))^4 * 1/4 +
            p^(-4/2) * (minY - sqrt(p))^5 * 1/5 -
            p^(-5/2) * (minY - sqrt(p))^6 * 1/6 +
            p^(-6/2) * (minY - sqrt(p))^7 * 1/7 ) -
        log(sqrt(p)) ) 
    return(out)
}
    
