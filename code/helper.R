
###
##
##
##

sqrtm <- function(Sigma){
    tmp = eigen(Sigma)
    return(tmp$vectors %*% diag(tmp$values^(1/2)) %*% t(tmp$vectors))
}

neg_sqrtm <- function(Sigma){
    tmp = eigen(Sigma)
    return(tmp$vectors %*% diag(tmp$values^(-1/2)) %*% t(tmp$vectors))
}


###########
## Uniform simulation

##
## Draws samples from the density
##
## 1/Z * r^{p-1} * I( r < p )
## 

rSSLCUnif <- function(n, p){
    M = 10*n;

    ## draw random points from a triangle with rejection sampling
    while (TRUE){
        init_pts = cbind(runif(M), runif(M))*sqrt(2)
        rtriangle = init_pts[ init_pts[, 1] > init_pts[, 2], 1]

        if (length(rtriangle) > n) break
    }
    
    ## samples from a density h(r) = r^(p-1) Id(r < 2^(1/p)) p/2
    rdens = rtriangle^(2/p)
    
    rdens = rdens*(p/2^(1/p))
    
    return(rdens[1:n])
}


## Sample N points uniformly at random
## from surface of a ball in p dimensions

runifBall <- function(N, p){
    X = matrix(rnorm(p*N), N, p)

    xnorms = apply(X, 1, norm, '2')

    for (ii in 1:N)
        X[ii, ] = X[ii, ]/xnorms[ii]
    
    ##X = diag(runif(N)^(1/p)) %*% X
    return(X)
}

## Sample N points uniformly at random from
## surface of a cube in p dimensions

runifSquare <- function(N, p){
    X = matrix(2*runif(p*N)-1, N, p)
    coord_max = apply(X, 1, function(x){max(abs(x))})

    for (ii in 1:N)
        X[ii, ] = X[ii, ]/coord_max[ii]
    return(X)
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
    
