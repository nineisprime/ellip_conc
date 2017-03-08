

c_ls = seq(30, 0, -0.1)

rat_ls = rep(0, length(c_ls))

for (ii in 1:length(c_ls)){
    c = c_ls[ii]

    rat = dnorm(x=0, mean=c)/pnorm(q=0, mean=c)

    rat_ls[ii] = rat
}
    
