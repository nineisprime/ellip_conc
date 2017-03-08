


## INPUT:
##   Y -- (nI--vector) must be sorted in increasing order
##   M -- scalar, num of boundary grid points to use
##
## OUTPUT:
##   evalpts -- vector, of length (2*length(bdpts) - 1)
##
##   bdpts -- (M + nI -2 vector) of boundary points, could be shorter
##            all nI Y points are boundary, as well as M grid points seq(minY, maxY)
##
##   ind_vec2 -- (nI-1 vector) of how many evalpts is in (Y_i, Y_{i+1}]
##                except first interval, which is [Y_1, Y_2]
##
##   y_rvec -- same length as evalpts. y_rvec[i] is the right-most Y value associated with
##             evalpts[i]
##   y_lvec --  similar to y_rvec
##   gap_vec -- similar to y_rvec
##                 "(y_rvec - evalpts)/gap_vec * phi_lvec + (evalpts - y_lvec)/gap_vec * phi_rvec"
##                 are the estimated "phi" value for any points between Ymin and Ymax
##
##   kepler_global_wts -- (2(M+nI)-1 vector)
##   kepler_ls -- list of kepler weights for each [Y_i, Y_{i+1}]
##   ixs_ls -- list of indices of evalpts in each [Y_i, Y_{i+1}]
##
## sum(ind_vec2) == length(evalpts)


## evalpts[ixs_ls[i]] contains all the evaluation points in [Y_i, Y_{i+1}]




numerical_integration_helper <- function(Y, M, eps){

    Y = sort(Y)
    nI = length(Y)
    
    grid_width = (max(Y) - min(Y))/(M-1)
    grid = seq(min(Y), max(Y), length.out=M)

    evalpts = c(grid, Y)
    evalpts = unique(evalpts)
    evalpts = sort(evalpts)
    
    nn = length(evalpts)
    midpts = (evalpts[2:nn] + evalpts[1:(nn-1)])/2
    bdpts = evalpts
    evalpts = sort(c(midpts, evalpts))

    ## ind_vec[i] is the number of boundary-evaluation points in [Y_1, Y_i]
    ## both right and left inclusive
    ind_vec = which(bdpts %in% Y)[2:nI]
    ##ind_vec = ceiling( (Y[2:nI] - min(Y))/grid_width - 1) + 2

    stopifnot(ind_vec[nI-1] == length(bdpts))
    
    ## for i=1, ind_vec2[1] is the number of eval points in [Y_1, Y_2]
    ## for i > 1, ind_vec2[i] is number of evaluation points in the (Y_{i-1}, Y_i]
    ## interval, right-inclusive
    ## there are (nI-1) such intervals, ind_vec2 has length (nI-1)
    if (nI > 2){
        ind_vec2 = c(ind_vec[1], ind_vec[2:(nI-1)] - ind_vec[1:(nI-2)] ) 
    } else {
        ind_vec2 = c(ind_vec[1]) 
    }

    ind_vec2 = ind_vec2*2           #count midpoints
    ind_vec2[1] = ind_vec2[1] - 1   #over counted for first cell

    
    ## y_rvec repeats each Y_i by ind_vec2[i] times
    y_rvec = rep2(ind_vec2, Y[2:nI])
    y_lvec = rep2(ind_vec2, Y[1:(nI-1)])
    gap_vec = rep2(ind_vec2, Y[2:nI] - Y[1:(nI-1)])


    ## precompute kepler_global_wts
    ## assumption: more than 2 boundary points globally
    kepler_global_wts = make_kepler_wts(bdpts)
    stopifnot(length(kepler_global_wts) == length(evalpts))

    
    ## precompute kepler_ls and ixs
    kepler_ls = list(0)
    ixs_ls = list(0)
    ## ixs_ls[[i]] is a list of indices of evalpts that lie inside (Y_i, Y_{i+1}]
    
    
    ix1 = 1
    bd_ix1 = 1
    
    for (i in 1:(nI-1)){
        if (i==1) {
            ixs = ix1:(ix1 + ind_vec2[i]-1)
            ix1 = ind_vec2[i] 
        } else {
            ixs = ix1:(ix1 + ind_vec2[i])
            ix1 = ix1 + ind_vec2[i]
        }
        ixs_ls[[i]] = ixs


        if (i==1){
            bd_ixs = bd_ix1:(bd_ix1 + (ind_vec2[i]-1)/2)
            bd_ix1 = bd_ix1 + (ind_vec2[i]-1)/2
        } else {
            bd_ixs = bd_ix1:(bd_ix1 + ind_vec2[i]/2)
            bd_ix1 = bd_ix1 + ind_vec2[i]/2
        }

        stopifnot(length(bd_ixs) == (length(ixs)+1)/2)
        stopifnot(max(bd_ixs) <= length(bdpts))
        
        cur_bdpts = bdpts[bd_ixs]
        
        kepler_wts = make_kepler_wts(cur_bdpts)
        kepler_ls[[i]] = kepler_wts

        
        stopifnot(kepler_wts >= 0)
        stopifnot(length(kepler_wts)==length(ixs))
    }

    return(list(evalpts = evalpts,
                bdpts = bdpts,
                ind_vec2 = ind_vec2,
                y_lvec = y_lvec,
                y_rvec = y_rvec,
                gap_vec = gap_vec,
                kepler_global_wts = kepler_global_wts,
                kepler_ls = kepler_ls,
                ixs_ls = ixs_ls))
}







## INPUT:
## bdpts -- (m--vector) of boundary points
##
## OUTPUT:
## kepler_wts -- (2m-1 vector) of kepler weights
##
## 1/6 * (delta1, 4*delta1, delta1+delta2, 4*delta2, delta2+delta3 ...)
##
## For Kepler numerical integration, which approximates
## int_a^b f(x)dx as (b-a)/6 * (f(a) + 4f( (a+b)/2 ) + f(b) )

make_kepler_wts <- function(bdpts) {
    m = length(bdpts)
    stopifnot(m > 1)
    
    kepler_gaps = bdpts[2:m] - bdpts[1:(m-1)]
    
    if (m == 2){
        kepler_wts = c(kepler_gaps, 4*kepler_gaps, kepler_gaps)
    } else {
        
        vec1 = 4*kepler_gaps
        vec2 = kepler_gaps[1: (length(kepler_gaps)-1) ] +
            kepler_gaps[2: (length(kepler_gaps)) ]
        
        tmp = rbind(vec1, c(vec2,0))
        dim(tmp) = NULL

        kepler_wts = c(kepler_gaps[1], tmp[1:(length(tmp)-1)],
            kepler_gaps[length(kepler_gaps)])/6
    }

    stopifnot(length(kepler_wts) == (2*m - 1))
    
    return(kepler_wts)
}






## Helper functions
##
## INPUT:
## ix_vec -- (n-vector) of non-negative integer
## val_vec -- (n-vector) of values
##
## OUTPUT:
## a vector of the form:
## [val_vec[1], val_vec[1], ...   #repeated ix_vec[1] times
##  val_vec[2], val_vec[2], ...   #repeated ix_vec[2] times
## ...


rep2 <- function(ix_vec, val_vec){

    stopifnot(length(ix_vec) == length(val_vec))
    
    if (length(ix_vec) == 1){
        return(rep(val_vec, ix_vec))
    }
    
    tmp = rbind(ix_vec, val_vec)
    out_vec = sapply( split(tmp, rep(1:ncol(tmp), each=nrow(tmp))), 
        FUN=function(x){rep(x[2], x[1])})
    out_vec = as.array(unlist(out_vec))

    return(out_vec)

}

