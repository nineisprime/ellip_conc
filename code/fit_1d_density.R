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


fit_1d_density <- function(Y, p, M=500000) {

    print("Starting ...")
    
    eps = 5e-5
    eps2 = 5e-4
    
    Y = sort(Y)
    n = length(Y)

    ygap = Y[2:n] - Y[1:(n-1)]

    ## create V matrix representing constraints
    ## V'phi <= 0
    diffmat = bandSparse(n-1, n, k=c(0,1),
        diagonals=list( -1/ygap, 1/ygap ))

    diffmat2 = bandSparse(n-2, n-1, k=c(0,1),
        diagonals=list(rep(-1,n-2), rep(1, n-2)))

    sec_diffmat = diffmat2 %*% diffmat
    sec_diffmat = rBind(c(-1/ygap[1], 1/ygap[1], rep(0,n-2)), sec_diffmat, rep(1, n))
    V = sec_diffmat

    ## rows of B represent dual-vectors
    ## Each V_i is a row of V
    ## Each B_i is a row of B
    ##
    ## B_i'V_i = -1
    ## B_i'V_j = 0
    
    B = -t(solve(V))
    
    V = V[1:(n-1), ]
    ##B = B[1:(n-1), ]
    
    
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

    ptm = proc.time()
    if (p < 200){
        phi = active_set_newton(Y, A, p, M=M)
    } else {
        phi = active_set_newton_largep(Y, A, p, M=M)

    }
    phi = expand_phi(phi, Y, A)

    print(proc.time() - ptm)
    ##ptm = proc.time()
    
    if (p < 200) {
        grad = compute_gradient_phi(phi, y_lvec=y_lvec, y_rvec=y_rvec, ind_vec2=ind_vec2,
            gap_vec=gap_vec, kepler_ls=kepler_ls, evalpts=evalpts, ixs_ls=ixs_ls, p=p)
    } else {
        grad = compute_gradient_phi_largep(phi, y_lvec=y_lvec, y_rvec=y_rvec, ind_vec2=ind_vec2,
            gap_vec=gap_vec, kepler_ls=kepler_ls, evalpts=evalpts, ixs_ls=ixs_ls, p=p)
    }

    ##print(proc.time() - ptm)
    
    slack = V %*% phi
    stopifnot(slack[A] > -eps2)

    
    safety = 1
    while (TRUE){
        safety = safety + 1
        stopifnot(safety < 200)

        dual_vec = B %*% grad

        I = setdiff(1:n, A)

        if (length(I) > 1){
            stopifnot(abs(dual_vec[I[1:(length(I)-1)]]) < eps)
        }
        
        violation_val = max(dual_vec[A])
        violation_ind = which(dual_vec == violation_val)

        if (violation_val < eps) break

        ## diagnostics
        print(paste("iteration: ", safety, "adding knot: ", violation_ind))

        ## add a knot (remove an active point)
        A = A[ !(A == violation_ind) ]
        I = setdiff(1:n, A)

        ##print(I)
        ##ptm = proc.time()
        if (p < 200){
            phi_prime = active_set_newton(Y, A, p, M=M, init_phi=phi[I])
        } else {
            phi_prime = active_set_newton_largep(Y, A, p, M=M, init_phi=phi[I])
        }
        ##print(proc.time() - ptm)
        print("Checkpoint")
        
        phi_prime = expand_phi(phi_prime, Y, A)

        slack_prime = V %*% phi_prime
        
        ##diagnostic
        ## for any i in A, V_i * phi_prime should be 0
        stopifnot(slack_prime[A] > -eps2)
        
        ## remove knots until phi_prime satisfy all constraints
        safety2 = 1
        while ( any(slack_prime > eps2) ){
            safety2 = safety2 + 1
            stopifnot(safety2 < 200)
            
            slack_violate = slack[slack_prime > eps2]
            slack_prime_violate = slack_prime[slack_prime > eps2]

            t_coef = max(-slack_prime_violate/(slack_violate - slack_prime_violate))
            
            phi = t_coef*phi + (1-t_coef)*phi_prime

            slack = V %*% phi
            
            A_new = which(slack > -eps2)

            ## diagnostics
            ## by removing knot, the active set should be larger
            stopifnot( all(A %in% A_new) )
            stopifnot(length(A_new) > length(A))
            
            rm_knot_ix = setdiff(A_new, A)
            print(paste("  remove knot: ", rm_knot_ix))

            A = A_new
            
            I = setdiff(1:n, A)
            if (p < 200){
                phi_prime = active_set_newton(Y, A, p, M=M, init_phi=phi[I])
            } else {
                phi_prime = active_set_newton_largep(Y, A, p, M=M, init_phi=phi[I])
            }
            
            phi_prime = expand_phi(phi_prime, Y, A)
            slack_prime = V %*% phi_prime
            
        }

    
        phi = phi_prime
        slack = slack_prime
        
        A = which(slack > -eps2)
        
        if (p < 200){
            grad = compute_gradient_phi(phi, y_lvec=y_lvec, y_rvec=y_rvec, ind_vec2=ind_vec2,
                gap_vec=gap_vec, kepler_ls=kepler_ls, evalpts=evalpts, ixs_ls=ixs_ls, p=p)
        } else {
          grad = compute_gradient_phi_largep(phi, y_lvec=y_lvec, y_rvec=y_rvec, ind_vec2=ind_vec2,
                gap_vec=gap_vec, kepler_ls=kepler_ls, evalpts=evalpts, ixs_ls=ixs_ls, p=p)
        }

    }
    return(phi)
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
                                 evalpts, ixs_ls, p){
  
    n = length(phi)

    phi_lvec = rep2(ind_vec2, phi[1:(n-1)])
    phi_rvec = rep2(ind_vec2, phi[2:n])
    
    fnpts = hx_eval(phi_lvec, phi_rvec, evalpts, y_lvec, y_rvec, gap_vec, p)

    grad = rep(1/n, n)
    grad[1] = grad[1] - hx1_int(phi[1], y_lvec[1], p)

    for (i in 1:(n-1)){
        ixs = ixs_ls[[i]]

        v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
        v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

        grad[i] = grad[i] - sum(v1 * kepler_ls[[i]] * fnpts[ixs])
        grad[i+1] = grad[i+1] - sum(v2 * kepler_ls[[i]] * fnpts[ixs])
    }

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

    
