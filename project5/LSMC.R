Monomial = function (x, k){
  if (k == 1)
    return(x^0)
  if (k == 2)
    return(x)
  if (k == 3)
    return(x^2)
  if (k == 4)
    return(x^3)
}

Hermit = function (x, k){
  if (k == 1)
    return(x^0)
  if (k == 2)
    return(2*x)
  if (k == 3)
    return(4*x^2-2)
  if (k == 4)
    return(8*x^3-12*x)
}

Lagur = function (x, k){
  if (k == 1)
    return(exp(-x/2))
  if (k == 2)
    return(exp(-x/2)*(1-x))
  if (k == 3)
    return(exp(-x/2)*(1- 2*x+x^2/2))
  if (k == 4)
    return(exp(-x/2)*(1 - 3*x + 3*x^2/2 -x^3/6))
}

LSMC = function (N, t, K, S0, r, sigma, num, k, FUN){
  paths = matrix(rep(0, len = num *(N+1)), nrow = num)
  rand = matrix(rnorm(num*N), nrow = num)
  CF = matrix(rep(0, len = num *(N)), nrow = num)
  dt = t/N
  paths[, 1] = S0
  zeros = rep(0 , len = num)
  zeros_row = rep(0, len = N)
  for (i in 1:N){
    paths[, i+1] = paths[, i] *exp((r-sigma^2/2)*dt + sigma*rand[,i]*sqrt(dt))
  }
  CF[, N] = apply(cbind(K - paths[,N+1], zeros), 1, max)
  for (i in (N-1) : 1){
    temp = (K - paths[,i+1] > 0)
    x = paths[, i+1] * temp
    dis_factor = exp(seq(-r*dt,-(N-i)*r*dt, -r*dt))
    if (i == N-1){
      y = CF[, (i+1):N] * dis_factor * temp
    }else{
      y = CF[, (i+1):N] %*% dis_factor * temp
    }
    pos = seq(1, num, 1) * temp
    y_clean = y[pos]
    x_clean = x[pos]
    A = matrix(rep(0, len = k*k), nrow = k)
    B = matrix(rep(0, len =k), nrow = k)
    X = c()
    for (a in 1:k){
      for (b in 1:k){
        A[a,b] = sum(FUN(x_clean, a) * FUN(x_clean, b))
      }
      B[a,1] = sum(FUN(x_clean,a) * y_clean)
      X = cbind(X, FUN(x, a))
    }
    reg = solve(A,B)
    continuation = X %*% reg * temp
    newCF = apply(cbind(K - paths[,i+1], zeros), 1, max)

    for (j in 1: num){
      if (newCF[j] > continuation[j]){
        CF[j, ] = zeros_row
        CF[j,i] = newCF[j]
      }
    }
  }
  D = exp(seq(-r*dt,-N*r*dt, -r*dt))
  return(mean(CF%*% D))
}


N = 250
S0 = 36
r = 0.06
sigma = 0.2
K = 40
t = 1
num = 100000
k = 3
#out = LSMC(N, t, K, S0, r, sigma, num, k, Monomial)
#out2 = LSMC(N, t, K, S0, r, sigma, num, k, Hermit)
out3 = LSMC(N, t, K, S0, r, sigma, num, k, Lagur)
