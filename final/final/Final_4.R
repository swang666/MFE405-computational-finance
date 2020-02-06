set.seed(12345)
r0 = 0.05
alpha = 0.36
beta = -5.86
sigma = 0.36
gamma = 2
t = 0
cap_t = 0.5
S = 1
K = 9800
face = 10000

dt = 0.01
N = S/dt

out = rep(0,10000)
for (k in 1: 10000){
  rand = rnorm(N)
  rt = rep(r0, N+1)
  for (i in 2:(N+1)){
    rt[i] = rt[i-1] + (alpha + beta*rt[i-1])*dt + sigma*rt[i-1]^gamma*sqrt(dt)*rand[i-1]
  }
  rTS0 = rt[N/2+1]
  rTS = rep(rTS0, N/2+1)
  P = rep(0, 100)
  for (i in 1: 100){
    for (j in 2:(N/2+1)){
      rTS[j] = rTS[j-1] + (alpha + beta*rTS[j-1])*dt + sigma*rTS[j-1]^gamma*sqrt(dt)*rnorm(1)
    }
    disrate = -sum(rTS[2:(N/2+1)])*dt
    PTS = face * exp(disrate)
    P[i] = PTS
  }
  disrate = -sum(rt[2:(N/2+1)])*dt
  payoff = exp(disrate)*max(K- mean(P),0) 
  out[k] = payoff
  
}
price = mean(out)


