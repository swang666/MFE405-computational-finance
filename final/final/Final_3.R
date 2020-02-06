set.seed(12345)
t = 5
dt = 0.01
N = t/dt
S0 = 100
r = 0.05
sigma = 0.35
K = 100

dts = seq(0, t, dt)
L = 50*exp(0.138629 * dts)
U = 200 - L


total = 0
ltotal = 0
payoffs = rep(0, 10000)
for (j in 1:10000){
  rand1 = rnorm(N)
  St = rep(S0, N+1)
  for (i in 2:(N+1)){
    St[i] = St[i-1] + St[i-1]*r*dt + St[i-1] * sigma* sqrt(dt) * rand1[i-1]
    
    if (St[i] <= L[i]){
      payoff = exp(-r*(i-1)*dt)*(K - St[i])
      payoffs[j] = payoff
      ltotal = ltotal + 1
      total = total + 1
      break
    }else if(St[i] >= U[i]){
      payoff = exp(-r*(i-1)*dt)*(St[i]- K)
      payoffs[j] = payoff
      total = total + 1
      break
    }
  }
}
price = mean(payoffs)
prob = ltotal/ total
