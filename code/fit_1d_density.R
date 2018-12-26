##
## main function: fit_1d_density
## helper functions:
## expand_phi, compute_gradient_phi
##
## workhorse:
##    active_set_newton
##    active_set_newton_largep


##
##
## INPUT:
##   Y -- (n vector)
##   p -- scalar
## OUTPUT:
##   phi -- (n vector) of phi-values
##

##
## returns phi : [0, inf) -> [-inf, inf) such that
## 
## int_0^inf e^phi(r) r^(p-1) dr = 1
## 

fit_1d_density <- function(Y, p, M=1000000) {

    print("Starting ...")
    
    EPS = 1e-7
    EPS2 = 1
    CONV_THRESH = 1e-10
    
    Y = sort(Y)
    n = length(Y)

    mu = mean(Y)
    sigma = sd(Y)
    
    Y = Y/sigma
    ## production-TODO: check for outlier and positivity
    mu = mean(Y)
    
    ygap = Y[2:n] - Y[1:(n-1)]

    
    ## precompute quantities needed for numerical integration
    num_int_res = numerical_integration_helper(Y, M)
    evalpts = num_int_res$evalpts
    bdpts = num_int_res$bdpts
    y_rvec = num_int_res$y_rvec
    y_lvec = num_int_res$y_lvec
    gap_vec = num_int_res$gap_vec
    ind_vec2 = num_int_res$ind_vec2
    kepler_global_wts = num_int_res$kepler_global_wts
    kepler_ls = num_int_res$kepler_ls
    ixs_ls = num_int_res$ixs_ls
         
    A = 2:(n-1)

    print("First run ...")

    ##ptm = proc.time()
    phi = active_set_newton_largep(Y, A, p, M=M, evalpts=evalpts, bdpts=bdpts)
    phi_n = expand_phi(phi, Y, A)

    grad = compute_gradient_phi_largep(phi_n, y_lvec=y_lvec,
                                       y_rvec=y_rvec, ind_vec2=ind_vec2,
                                       gap_vec=gap_vec, kepler_ls=kepler_ls,
                                       evalpts=evalpts, ixs_ls=ixs_ls,
                                       mu=mu, p=p)
    


    
    safety = 1
    while (TRUE){
        safety = safety + 1
        stopifnot(safety < 200)

        grad_bar = grad - mean(grad)
        tmp = rev(cumsum(rev( grad_bar[2:n] )))
        dual_vec = rev(cumsum(rev( -ygap*tmp )))
        dual_vec = dual_vec/sqrt(seq(n-1, 1, -1))
        
        I = setdiff(1:n, A)

        ## ASSERT
        ## by optimality, b_i' gradient = 0 for loose constraints i
        if (length(I) > 1){
            dual_vec_I = dual_vec[I[1:(length(I)-1)]]
            stopifnot(max(abs(dual_vec_I)) < EPS2)
        }

        violation_val = max(dual_vec[A])
        violation_ind = which(dual_vec == violation_val)
        if (violation_val < CONV_THRESH) break

        ## diagnostics
        print(paste("iteration:", safety, "adding knot:", violation_ind,
                    "violation:", violation_val))

        
        ## add a knot (remove an active point)
        Aprime = A[ !(A == violation_ind) ]
        Iprime = setdiff(1:n, Aprime)

        
        phi_prime = active_set_newton_largep(Y, Aprime, p, M=M,
                                             evalpts=evalpts, bdpts=bdpts,
                                             init_phi=phi_n[Iprime])

        YIprime = Y[Iprime]
        nIprime = length(Iprime)
        Vprime = createV(YIprime[2:nIprime] - YIprime[1:(nIprime-1)], nIprime)

        slack_prime = Vprime %*% phi_prime

        newix = which(Iprime == violation_ind)
        phi_tmp = rep(0, nIprime)
        phi_tmp[setdiff(1:nIprime, c(newix))] = phi
        if (newix == 1)
            phi_tmp[newix] = phi[1]
        else
            phi_tmp[newix] = ( (YIprime[newix+1] - YIprime[newix])*phi_tmp[newix-1] +
                                (YIprime[newix] - YIprime[newix-1])*phi_tmp[newix+1] )/
                               (YIprime[newix+1] - YIprime[newix-1]) 
        
        
        slack = Vprime %*% phi_tmp

        ## ASSERT
        ## for any i in A, V_i * phi_prime should be 0
        if (slack_prime[newix] > EPS){
            print("Secondary finish")
            return(phi_n - (p-1)*log(mu) - p*log(sigma))
        }
        
      
        ## remove knots (add tight constraints A)
        ## until phi_prime satisfy all constraints
        while ( any(slack_prime > EPS) ){

            ## LOOP invariant:
            ##           phi_prime
            ##           Iprime
            ##   slack,  slack_prime
            ##
            ## all on the Iprime space

            slack_violate = slack[which(slack_prime > EPS)]
            slack_prime_violate = slack_prime[which(slack_prime > EPS)]

            tmp = (-slack_prime_violate/(slack_violate - slack_prime_violate))
            t_coef = max(tmp)

            rm_knot_ix_Iprime = which(slack_prime == slack_prime_violate[which.max(tmp)])
            keep_ixs = which(slack_prime != slack_prime_violate[which.max(tmp)])

            rm_knot_ix = Iprime[rm_knot_ix_Iprime]
            print(paste("  remove knot: ", rm_knot_ix))
            I_new = c(Iprime[keep_ixs], n)
            A_new = setdiff(1:n, I_new)

            slack = t_coef*slack + (1-t_coef)*slack_prime
            slack = slack[keep_ixs]

            phi_prime = active_set_newton_largep(Y, A_new, p, M=M,
                                                 evalpts=evalpts, bdpts=bdpts,
                                                 init_phi=phi_n[I_new])
            Iprime = I_new
            nIprime = length(Iprime)
            YIprime = Y[Iprime]

            if (nIprime == 1) break
            
            Vprime = createV(YIprime[2:nIprime] - YIprime[1:(nIprime-1)], nIprime)
            slack_prime = Vprime %*% phi_prime
        }
        ## Loop outputs phi_prime, dimension Iprime

        phi = phi_prime
        A = setdiff(1:n, Iprime)
        phi_n = expand_phi(phi, Y, A)
        grad = compute_gradient_phi_largep(phi_n, y_lvec=y_lvec,
                                           y_rvec=y_rvec, ind_vec2=ind_vec2,
                                           gap_vec=gap_vec, kepler_ls=kepler_ls,
                                           evalpts=evalpts, ixs_ls=ixs_ls,
                                           mu=mu, p=p)
        
    }
    return(phi_n - (p-1)*log(mu) - p*log(sigma))
}

## INPUT:
##  phi -- (nI vector)
##  y_lvec, y_rvec (length(evalpts) vector)
##  
##
## OUTPUT:
##   grad -- (nI vector) of gradient
##

compute_gradient_phi <- function(phi, y_lvec, y_rvec, gap_vec, kepler_ls, ind_vec2,
                                 evalpts, ixs_ls, p){
        
    n = length(phi)

    phi_lvec = rep2(ind_vec2, phi[1:(n-1)])
    phi_rvec = rep2(ind_vec2, phi[2:n])
    
    fnpts =  exp(
        ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
        (p-1)*log(evalpts/sqrt(p)) )

    grad = rep(1/n, n)
    grad[1] = grad[1] - exp(phi[1] + p*log(y_lvec[1]/sqrt(p)) - (1/2)*log(p))

    for (i in 1:(n-1)){
        ixs = ixs_ls[[i]]

        v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
        v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

        grad[i] = grad[i] - sum(v1 * kepler_ls[[i]] * fnpts[ixs])
        grad[i+1] = grad[i+1] - sum(v2 * kepler_ls[[i]] * fnpts[ixs])
    }

    return(grad)
}


## Same as before but for p large
compute_gradient_phi_largep <- function(phi, y_lvec, y_rvec, gap_vec, kepler_ls, ind_vec2,
                                 evalpts, ixs_ls, mu, p){
  
    n = length(phi)

    phi_lvec = rep2(ind_vec2, phi[1:(n-1)])
    phi_rvec = rep2(ind_vec2, phi[2:n])
    
    fnpts = hx_eval(phi_lvec, phi_rvec, evalpts, y_lvec, y_rvec, gap_vec, mu, p)

    grad = rep(1/n, n)
    grad[1] = grad[1] - hx1_int(phi[1], y_lvec[1], mu, p)
    
    ## for (i in 1:(n-1)){
    ##     ixs = ixs_ls[[i]]

    ##     v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
    ##     v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

    ##     grad[i] = grad[i] - sum(v1 * kepler_ls[[i]] * fnpts[ixs])
    ##     grad[i+1] = grad[i+1] - sum(v2 * kepler_ls[[i]] * fnpts[ixs])
    ## }
    
    a_vec = rep(0, n-1)
    b_vec = rep(0, n-1)

    K = 300
    lengths = sapply(ixs_ls, FUN=length)
    sm_gaps = which(lengths <= K)
    big_gaps = which(lengths > K)
    
    ixs_sm_gaps =  lapply(ixs_ls[sm_gaps],
                            function(x){length(x) = K; return(x)})
    ixs = do.call("cbind", ixs_sm_gaps)
    ixs[is.na(ixs)] = 1

    kepler_sm_gaps = lapply(kepler_ls[sm_gaps],
                            function(x){length(x) = K; return(x)})
    kepler_mat = do.call("cbind", kepler_sm_gaps)
    kepler_mat[is.na(kepler_mat)] = 0
    
    v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
    v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

    v2_mat = matrix(v2, nrow(ixs), ncol(ixs))
    v1_mat = matrix(v1, nrow(ixs), ncol(ixs))

    fn_mat = matrix(fnpts[ixs], nrow(ixs), ncol(ixs))

    a_vec[sm_gaps] = apply(- v1_mat * kepler_mat * fn_mat, MARGIN=2, FUN=sum)
    b_vec[sm_gaps] = apply(- v2_mat * kepler_mat * fn_mat, MARGIN=2, FUN=sum)
    
    for (i in big_gaps){
        ixs = ixs_ls[[i]]

        v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
        v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

        a_vec[i] = - sum(v1 * kepler_ls[[i]] * fnpts[ixs])
        b_vec[i] = - sum(v2 * kepler_ls[[i]] * fnpts[ixs])
    }

    grad[1:(n-1)] = grad[1:(n-1)] + a_vec
    grad[2:n] = grad[2:n] + b_vec
    
    ## afn <- function(i){
    ##     ixs = ixs_ls[[i]]
    ##     v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]
    ##     return(- sum(v1 * kepler_ls[[i]] * fnpts[ixs]))
    ## }
    
    ## bfn <- function(i){
    ##     ixs = ixs_ls[[i]]
    ##     v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs]
    ##     return(- sum(v2 * kepler_ls[[i]] * fnpts[ixs]))
    ## }
        
    ## a_coefs = sapply(1:(n-1), afn)
    ## b_coefs = sapply(1:(n-1), bfn)

    ## grad[1:(n-1)] = grad[1:(n-1)] + a_coefs
    ## grad[2:n] = grad[2:n] + b_coefs
    
    return(grad)
}




## INPUT:
##   phi -- (nI vector)
##   A -- (n - nI vector) of active coordinates
##   Y -- (n vector) sorted
##
##  
##
## OUTPUT:
##  

## A must not contain "n"

expand_phi <- function(phi, Y, A){

    n = length(Y)
    nI = length(phi)

    I = setdiff(1:n, A)
    YI = Y[I]
    
    if (nI == 1){
        phi_out = rep(phi, n)
        return(phi_out)
    }

    phi_out = rep(0, n)
    phi_out[1:I[1]] = phi[1]
    
    ind_vec = I[2:nI] - I[1:(nI-1)]
    y_lvec = rep2(ind_vec, YI[1:(nI-1)])
    y_rvec = rep2(ind_vec, YI[2:nI])
    
    gap_vec = y_rvec - y_lvec

    phi_lvec = rep2(ind_vec, phi[1:(nI-1)])
    phi_rvec = rep2(ind_vec, phi[2:nI])

    ixs = (I[1] + 1):n
    phi_out[ixs] = (y_rvec - Y[ixs])/gap_vec * phi_lvec +
        (Y[ixs] - y_lvec)/gap_vec * phi_rvec

    
    stopifnot(sum(abs(phi_out[I] - phi)) < 1e-6)
    
    return(phi_out)
}

    
## create V matrix

createV <- function(gaps, nI) {

    if (nI == 2){
        V = matrix(c(-1, 1), 1, 2)
        return(V)
    }
    
    
    diffmat = bandSparse(nI-1, nI, k=c(0,1),
                         diagonals=list( -1/gaps, 1/gaps ))

    diffmat2 = bandSparse(nI-2, nI-1, k=c(0,1),
                          diagonals=list(rep(-1,nI-2), rep(1, nI-2)))

    sec_diffmat = diffmat2 %*% diffmat
    sec_diffmat = rBind(c(-1/gaps[1], 1/gaps[1], rep(0,nI-2)), sec_diffmat)
    V = sec_diffmat

    return(V)
}
