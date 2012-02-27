`rickerplot` <-
    function(N,R,K,b)
{
    i <- 1
    a <- ((R-1)^(1/b))/K
    over <- 0
    tmp <- N
    repeat { 
        N <- (N*R)/(1+(a*N)^b)
        tmp <- c(tmp,N)
        i <- i+1
        if (N > 0.99*K) over <- over+1
        if (over > 10) break
    }
    steps <- i
                                        #  par(mfrow=c(1,2))
    plot(tmp,type="l",xlab="t",ylab="N")
    abline(h=K,col="red")
    x <- seq(from=0,to=max(tmp)+0.5,len=100)
    plot(x,(x*R)/(1+(a*x)^b),type="l",xlab=expression(N[t]),ylab=expression(N[t+1]))
    lines(x,x,type="l",lty=2)
    for (i in 1:(steps-1)) {
        now <- tmp[i]
        nex <- tmp[i+1]
        segments(now,now,now,nex,col="red")
        segments(now,nex,nex,nex,col="red")
        i <- i+1
    }
}
