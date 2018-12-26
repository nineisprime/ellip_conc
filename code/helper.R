
## Sample N points uniformly at random
## from surface of a ball in p dimensions

runifBall <- function(N, p){
    X = matrix(rnorm(p*N), N, p)
    X = diag(1/apply(X, 1, norm, '2')) %*% X
    #X = diag(runif(N)^(1/p)) %*% X
}

## Sample N points uniformly at random from
## surface of a cube in p dimensions

runifSquare <- function(N, p){
    X = matrix(2*runif(p*N)-1, N, p)
    coord_max = apply(X, 1, function(x){max(abs(x))})
    X = diag(1/coord_max) %*% X
}


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
    
    return(out)
}
    
