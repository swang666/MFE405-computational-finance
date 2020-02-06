
t = 2
dt = 0.01
N = t/dt
set.seed(12345)
v0 = 0.06
alpha = 0.45
beta = -5.105
gamma = 0.25
S0 = 20
r = 0.05
rhos = c(-0.75, 0, 0.75)

generate_bivariate = function (rand1, rand2, rho) {
  
  y = rho * rand1 + sqrt(1 - rho*rho) * rand2
  
  return(y)
}


calc_price = function(rho){
  payoffs = rep(0, 10000)
  for (j in 1:10000){
    rand1 = rnorm(N)
    rand2 = rnorm(N)
    z1 = rand1
    z2 = generate_bivariate(rand1, rand2, rho)
    
    vt = rep(v0, N+1)
    St = rep(S0, N+1)
    for (i in 2:(N+1)){
      vt[i] = vt[i-1] +(alpha + beta*max(vt[i-1],0))* dt + gamma*sqrt(max(vt[i-1],0))*sqrt(dt)*z1[i-1]
      St[i] = St[i-1] + r*St[i-1]*dt + St[i-1] * sqrt(max(vt[i-1],0))* sqrt(dt)*z2[i-1]
    }
    A = mean(St)
    payoff = max(St[N+1] - A, 0)
    payoffs[j] = payoff
  }
  return(exp(-r*t)*mean(payoffs))
}

price1 = calc_price(rhos[1])
price2 = calc_price(rhos[2])
price3 = calc_price(rhos[3])