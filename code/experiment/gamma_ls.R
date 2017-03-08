
library(ggplot2)
library(reshape2)

THRESH = 1e-7

p_ls = c(1000, 2000, 4000, 8000, 16000, 32000, 64000)
r0_ls = seq(0.045, 1, 0.002)

envels = matrix(0, length(r0_ls), length(p_ls))

## a_ls = rep(0, length(k_ls))


for (jj in 1:length(p_ls)){
    p = p_ls[jj]

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

            
            ##print(paste(a/sqrt(p), x1/x2))
            
            if (x1 == 0 && x2 == 0)
                Fval = (a^2/((p+1)*p))
            else {
                Fval = (x1/x2)* ( a^2/((p+1)*p))
            }
            
            
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
        envels[ii, jj] = envel

    }

}



##whole plot
df = data.frame(cbind(r0_ls, envels))
colnames(df) = c("r0", "1000", "2000", "4000", "8000", "16000", "32000", "64000")
melted = melt(df, id.vars="r0")
ggplot() + geom_line(data=melted, aes(x=r0, y=value, colour=variable))



## 1/x
df = data.frame(cbind(r0_ls, 1/envels))
colnames(df) = c("r0", "1000", "2000", "4000", "8000", "16000", "32000", "64000")
melted = melt(df, id.vars="r0")
ggplot() + geom_line(data=melted, aes(x=r0, y=value, colour=variable))



## log(x)
df = data.frame(cbind(r0_ls, log(envels)))
colnames(df) = c("r0", "1000", "2000", "4000", "8000", "16000", "32000", "64000")
melted = melt(df, id.vars="r0")
ggplot() + geom_line(data=melted, aes(x=r0, y=value, colour=variable))



## 1/x on 0.05 to 0.1
subs = 1:50
df = data.frame(cbind(r0_ls[subs], 1/envels[subs,]))
colnames(df) = c("r0", "1000", "2000", "4000", "8000", "16000", "32000", "64000")
melted = melt(df, id.vars="r0")
ggplot() + geom_line(data=melted, aes(x=r0, y=value, colour=variable))


## log(x) on 0.64 to 1
subs = 300:length(r0_ls)
df = data.frame(cbind(r0_ls[subs], log(envels[subs,])))
colnames(df) = c("r0", "1000", "2000", "4000", "8000", "16000", "32000", "64000")
melted = melt(df, id.vars="r0")
ggplot() + geom_line(data=melted, aes(x=r0, y=value, colour=variable))


