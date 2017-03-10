

x_2p_minus_1_e_xsquared <- function(Y, p) {

    if (p > 100){
        tmp = (2*p - 1)*log(Y/sqrt(p-1)) + (1 - (Y/sqrt(p-1))^2)*(p-1) + log(2) -
            (1/2)*log(2*pi*(p-1)) + (1/2)*log(p-1)
    } else {
        tmp = (2*p - 1)*log(Y) - Y^2 + log(2) -
            (p-1)*log(p-1) + (p-1) - (1/2)*log(2*pi*(p-1))
    }

    return(exp(tmp))
}


