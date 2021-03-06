
THRESH = 1e-6

p = 1000

##k_ls = seq(5, 400, 2)
r0_ls = seq(4/sqrt(p), 1, 0.003)

envels = rep(0, length(k_ls))

a_ls = rep(0, length(k_ls))

prev_mult = 1

for (ii in 1:length(r0_ls)) {

    r0 = r0_ls[ii]
    loc = sqrt(p) + r0

    hi = prev_mult
    lo = 0

    while (TRUE){

        a = sqrt(p*(1 - ((hi - lo)/2 + lo)))

        ##print(paste(hi, lo))
        
        x1 = pgamma(loc, shape=p, rate=a)
        x2 = pgamma(loc, shape=p+2, rate=a)

        if (x1 == 0 && x2 == 0)
            Fval = (a^2/((p+1)*p))
        else {
            Fval = (x1/x2)* ( a^2/((p+1)*p))
        }
        
        #print(val)
        #print(1/p)
        
        if (Fval/(1/p) - 1 > THRESH){
            lo = (hi - lo)/2 + lo
        } else if (Fval/(1/p) - 1 < -THRESH){
            hi = (hi - lo)/2 + lo
        } else {
            prev_mult = (hi - lo)/2 + lo
            break
        }
        
    }


    
    envel = dgamma(loc, shape=p, rate=a)/pgamma(loc, shape=p, rate=a) 

    envels[ii] = envel

}

pdf(paste0(p, ".pdf"))
plot(r0_ls, log(envels), type="l")
dev.off()

pdf(paste0(p, "b.pdf"))
plot(r0_ls, 1/envels, type="l")
dev.off()
## if val > 1/p, then loc should be more to the left
## if val < 1/p, then loc should be more to the right

