set.seed(12345)
result = rep(0,10000)
for (i in 1:10000){
  rand = runif(1000)
  for (j in 1:1000){
    temp = sum(rand[1:j])
    if (temp > 1.1){
      result[i] = j
      break
    }
  }
}

out = 4.54 - result
out = out*(out >= 0)
answer = mean(out)
answer