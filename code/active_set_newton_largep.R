
## Main function:
## active_set_newton
##
## Helper functions:
## rep2,  make_kepler_wts,  numerical_integration_helper
##
##

## INPUT:
## Y -- (n--vector) of positive samples
## A -- (vector) of indices of active constraints
## M -- scalar, number of grid points to perform numerical integration
##
## OUTPUT:
## phi -- ( length(I)-vector ) of function values
##

active_set_newton_largep <- function(Y, A, p, M=800000, init_phi=NULL, diagnostic=FALSE){

    n = length(Y)
    I = setdiff(1:n, A)
    nI = length(I)

    Y_sorted = sort(Y)
    
    eps = 1e-5
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
    ##

    if (is.null(init_phi)){
        phi = - ((p-1)/sqrt(p)) * (YI_sorted - sqrt(p))
    } else {
        phi = init_phi
    }

    ## precompute quantities needed for numerical integration
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

    
    safety = 1
    stepsize = 1
    while (TRUE && safety < 2000) {
        safety = safety + 1
        stopifnot(safety < 1000)
        
        ## compute function evaluations
        ## fnpts[i] is function evaluated at evalpts[i]
        if (nI > 1){
            phi_lvec = rep2(ind_vec2, phi[1:(nI-1)])
            phi_rvec = rep2(ind_vec2, phi[2:nI])

            fnpts =  exp(
                ((y_rvec-evalpts)*phi_lvec + (evalpts-y_lvec)*phi_rvec )/gap_vec +
                ((p-1)/sqrt(p)) * (evalpts - sqrt(p)) -
                ((p-1)/p) * (evalpts - sqrt(p))^2 * 1/2 +
                ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^3 * 1/3 -
                ((p-1)/p^2) * (evalpts - sqrt(p))^4 * 1/4 +
                ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^5 * 1/5 -
                ((p-1)/p^3) * (evalpts - sqrt(p))^6 * 1/6 +
                ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^7 * 1/7)

            stopifnot( y_rvec - evalpts >= 0)
            stopifnot( evalpts - y_lvec >= 0)
            stopifnot( y_rvec - y_lvec == gap_vec)
            ##plot(fnpts)
            
            stopifnot(fnpts >= 0)
        }
        
        ## compute gradient and Hessian
        grad = grad_ypart
        Hdiag = rep(0, nI)
        Hdiag_off = rep(0, nI-1)

        fn_int = exp(phi[1] +
                sqrt(p)*(min(YI_sorted) - sqrt(p)) -
                (min(YI_sorted) - sqrt(p))^2 * 1/2 +
                (1/sqrt(p)) * (min(YI_sorted) - sqrt(p))^3 * 1/3 -
            (1/p) * (min(YI_sorted) - sqrt(p))^4 * 1/4 +
            (1/p^(3/2)) * (min(YI_sorted) - sqrt(p))^5 * 1/5 -
            (1/p^2) * (min(YI_sorted) - sqrt(p))^6 * 1/6 +
            (1/p^(5/2)) * (min(YI_sorted) - sqrt(p))^7 * 1/7 -
            log(sqrt(p)))
        
        grad[1] = grad[1] - fn_int
                
                            
        Hdiag[1] = Hdiag[1] - fn_int
        

        
        for (i in 1:(nI-1)){
            if (nI == 1) break
            
            ixs = ixs_ls[[i]]

            v2 = (evalpts[ixs] - y_lvec[ixs])/gap_vec[ixs] 
            v1 = (y_rvec[ixs] - evalpts[ixs])/gap_vec[ixs]

            fn_int = fn_int + sum(kepler_ls[[i]] * fnpts[ixs])

            grad[i] = grad[i] - sum(v1 * kepler_ls[[i]] * fnpts[ixs])
            grad[i+1] = grad[i+1] - sum(v2 * kepler_ls[[i]] * fnpts[ixs])

            Hdiag[i] = Hdiag[i] - sum(v1^2 * kepler_ls[[i]] * fnpts[ixs])
            Hdiag[i+1] = Hdiag[i+1] - sum(v2^2 * kepler_ls[[i]] * fnpts[ixs])

            Hdiag_off[i] = Hdiag_off[i] - sum(v1*v2 * kepler_ls[[i]] * fnpts[ixs])
        }


        if (nI > 1){
            H = bandSparse(nI, k=-1:1, diag=list(Hdiag_off, Hdiag, Hdiag_off))
        } else{
            H = Hdiag[1]
        }

        
        ## test for convergence
        if (sqrt(sum(grad^2)) < 1e-7)
            break

        
        ##
        phi_step = -solve(H, grad, sparse=TRUE)
        ##phi_step[phi_step > 300] = (phi_step[phi_step > 300])^(1/3)
        ##phi_step[phi_step < -300] = -(-phi_step[phi_step < -300])^(1/3)


        ## diagnostics
        if (diagnostic){
            print(safety)
            print(fn_int)
            if (nI > 1)
                plot(fnpts[seq(1, length(fnpts), by=100)])
            
            print("phi_step")
            print(phi_step)
            print("phi")
            print(phi)
            print("grad")
            print(grad)
            print("stepsize")
            print(stepsize)
            print("")
        }
        
        
        
        ## debug
        if (any(abs(phi_step) > 500)){
            phi_step = 100*grad
        }
        
        old_obj = sum(grad_ypart*phi) - fn_int
        new_obj = -Inf

        ## backtracking line search
        sigma = 0.02
        step_beta = 0.9
        stopifnot( sum(phi_step*grad) >= 0 )
        armijo_safety = 1
        while (new_obj - old_obj < sigma * stepsize * sum(phi_step*grad) - 1e-7){
            if (new_obj > -Inf)
                stepsize = step_beta * stepsize

            phi_cand = phi + stepsize*phi_step

            ## recompute new_obj
            if (nI > 1){
                phi_cand_lvec = rep2(ind_vec2, phi_cand[1:(nI-1)])
                phi_cand_rvec = rep2(ind_vec2, phi_cand[2:nI])
                
                fnpts_cand =  exp(
                    ((y_rvec-evalpts)*phi_cand_lvec +
                     (evalpts-y_lvec)*phi_cand_rvec )/gap_vec +
                    ((p-1)/sqrt(p)) * (evalpts - sqrt(p)) -
                    ((p-1)/p) * (evalpts - sqrt(p))^2 * 1/2 +
                    ((p-1)/p^(3/2)) * (evalpts - sqrt(p))^3 * 1/3 -
                    ((p-1)/p^2) * (evalpts - sqrt(p))^4 * 1/4 +
                    ((p-1)/p^(5/2)) * (evalpts - sqrt(p))^5 * 1/5 -
                    ((p-1)/p^3) * (evalpts - sqrt(p))^6 * 1/6 +
                    ((p-1)/p^(7/2)) * (evalpts - sqrt(p))^7 * 1/7)

                fn_int_cand = sum(fnpts_cand * kepler_global_wts) +
                    exp(phi_cand[1] +
                        sqrt(p)*(min(YI_sorted) - sqrt(p)) -
                        (min(YI_sorted) - sqrt(p))^2 * 1/2 +
                        (1/sqrt(p)) * (min(YI_sorted) - sqrt(p))^3 * 1/3 -
                        (1/p) * (min(YI_sorted) - sqrt(p))^4 * 1/4 +
                        (1/p^(3/2)) * (min(YI_sorted) - sqrt(p))^5 * 1/5 -
                        (1/p^2) * (min(YI_sorted) - sqrt(p))^6 * 1/6 +
                        (1/p^(5/2)) * (min(YI_sorted) - sqrt(p))^7 * 1/7 -
                        log(sqrt(p)))
                        
                
            } else {
                fn_int_cand =  exp(phi_cand[1] +
                    sqrt(p)*(min(YI_sorted) - sqrt(p)) -
                    (min(YI_sorted) - sqrt(p))^2 * 1/2 +
                    (1/sqrt(p)) * (min(YI_sorted) - sqrt(p))^3 * 1/3 -
                    (1/p) * (min(YI_sorted) - sqrt(p))^4 * 1/4 +
                    (1/p^(3/2)) * (min(YI_sorted) - sqrt(p))^5 * 1/5 -
                    (1/p^2) * (min(YI_sorted) - sqrt(p))^6 * 1/6 +
                    (1/p^(5/2)) * (min(YI_sorted) - sqrt(p))^7 * 1/7 -
                    log(sqrt(p)))
            }
            new_obj = sum(grad_ypart*phi_cand) - fn_int_cand

            armijo_safety = armijo_safety + 1
            ##if (armijo_safety > 100)
            ##    break
            stopifnot(armijo_safety < 100)
            
        }

            
        ## update phi
        phi = phi_cand

        ## alternative convergence criterion
        ##if (new_obj - old_obj < 1e-10)
        ##   break
    }

    return(phi)
}

