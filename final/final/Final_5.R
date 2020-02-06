set.seed(12345)
r = 0.05
q =0
rf = 0.04
sigma1 = 0.1
sigma2 = 0.15
gamma = -0.04
lambda = 1.5
K = 60
t = 1
dt = 0.01
N = t/dt
rho = -0.25
S0 = 6000
E0 = 0.0096

generate_bivariate = function (rand1, rand2, rho) {
  y = rho * rand1 + sqrt(1 - rho*rho) * rand2
  return(y)
}

out = rep(0, 10000)
for (j in 1 : 10000){
  rand1 = rnorm(N)
  rand2 = rnorm(N)
  z1 = rand1
  z2 = generate_bivariate(rand1, rand2, rho)
  
  St = rep(S0, N+1)
  Et = rep(E0, N+1)
  
  exp1 = rexp(1, lambda)
  for (i in 2:(N+1)){
    St[i] = St[i-1]+St[i-1]*(r-q)*dt + St[i-1]*sigma1 *sqrt(dt)*z1[i-1]
    Et[i] = Et[i-1] + (r-rf)*Et[i-1]*dt + sigma2 * Et[i-1]*sqrt(dt)* z2[i-1]
    curr = dt*(i-1)
    if (curr > exp1){
      exp1 = exp1 + rexp(1, lambda)
      St[i] = St[i] * (1 + gamma)
    }
  }
  payoff = max(St[N+1]*Et[N+1]-K, 0)
  out[j] = payoff
}
price = exp(-r*t)*mean(out)