"logrow" <-
  function(N0,b,d,K,bcontr=0.5,out=0) {
    c <- (K-N0)/N0
    r <- b-d
    t <- 0
    N <- N0
    while(N[t+1] < 0.998*K) {
      t <- t+1
      tmp <- K/(1+c*exp(-r*t))
      N <- c(N,tmp)
    }
    tend <- length(N)
    t <- seq(0,length(N)-1)
    plot(t,N,type="l",col="red",ylim=c(0,1.2*K))
    abline(h=K)
    time <- 0
    now <- 0
    simN <- N0
    size <- N0
    while (now < tend) {
      adj <- (b-d)*(size/K)
      bN <- b - bcontr*adj
      dN <- d + (1-bcontr)*adj
      events <- (bN+dN)*size
      step <- rexp(1,events)
      now <- now+step
      if (runif(1) <= bN/(bN+dN)) {
        size <- size+1
      }
      else {
        size <- size-1
      }
      simN <- c(simN,size)
      time <- c(time,now)
      if (size == 0) break
    }
    lines(time,simN,col="blue")
    final <- cbind(time,simN)
    if (out != 0) return(final)
  }




