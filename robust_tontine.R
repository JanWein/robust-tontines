# Reproducible checks for the robust heterogeneous tontine manuscript.
# Base R only. Run from this directory: Rscript robust_tontine.R 10000 results
# Exact event-driven simulation with constant individual death intensities.
args <- commandArgs(trailingOnly=TRUE)
M <- if(length(args)) as.integer(args[1]) else 10000L
out <- if(length(args)>1) args[2] else "results"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261002)

degrees <- function(k, epsilon=0, eta=Inf) {
  b <- if(epsilon==0) k else pmin(k,eta*k/(2*epsilon))
  B <- sum(b)
  pmax(0,pmin(b,B-b))
}
flux_col <- function(d,j) {
  L <- sum(d)
  if(L==0) return(d*0)
  end <- cumsum(d); start <- end-d
  s <- (start[j]+L/2)%%L; e <- s+d[j]
  v <- pmax(0,pmin(end,min(e,L))-pmax(start,s))
  if(e>L) v <- v+pmax(0,pmin(end,e-L)-start)
  v[j] <- 0
  v
}
flux <- function(d) matrix(vapply(seq_along(d),function(j) flux_col(d,j),numeric(length(d))),length(d))

# Exhaustive vertices of the uncertainty box for small n.
check_states <- function() {
  ans <- list(); l <- 0L
  for(n in c(1,2,3,5,10,20)) for(rep in seq_len(30)) {
    k <- exp(rnorm(n,0,3)); delta <- exp(rnorm(n,0,3)); eps <- 0.2
    b <- pmin(k,delta/(2*eps)); B <- sum(b); K <- sum(k)
    d <- pmax(0,pmin(b,B-b)); F <- flux(d)
    scale <- max(1,K)
    row_error <- max(abs(rowSums(F)-d))/scale
    sym_error <- max(abs(F-t(F)))/scale
    frontier <- K-B+max(0,2*max(b)-B)
    frontier_error <- abs(K-sum(d)-frontier)/scale
    worst_error <- NA_real_
    if(n<=5) {
      U <- as.matrix(expand.grid(rep(list(c(-eps,eps)),n)))
      G <- U%*%t(F)-sweep(U,2,d,"*")
      worst_error <- max(abs(apply(abs(G),2,max)-2*eps*d))/scale
    }
    stopifnot(row_error<1e-11,sym_error<1e-11,frontier_error<1e-11,
              min(F)>=0,all(diag(F)==0),all(d<=b+1e-10*scale))
    if(!is.na(worst_error)) stopifnot(worst_error<1e-11)
    l <- l+1L; ans[[l]] <- data.frame(n,rep,row_error,sym_error,frontier_error,worst_error)
  }
  do.call(rbind,ans)
}

# Each account is explicitly liquidated at T. This validates finite-horizon
# identities including terminal capital, not an infinite-horizon payout claim.
path <- function(pop,z,epsilon,eta,exit_rate,T=20,r=0.02,check=5) {
  n <- nrow(pop); mu <- pop$mu; entry <- pop$entry; capital <- pop$capital
  death <- entry+rexp(n,mu*z)
  exit <- if(exit_rate>0) entry+rexp(n,exit_rate) else rep(Inf,n)
  cp <- entry+check
  ev <- rbind(data.frame(t=entry,type=1L,id=seq_len(n)),
              data.frame(t=death,type=2L,id=seq_len(n)),
              data.frame(t=exit,type=3L,id=seq_len(n)),
              data.frame(t=cp,type=4L,id=seq_len(n)),
              data.frame(t=T,type=5L,id=0L))
  ev <- ev[ev$t<=T,]; ev <- ev[order(ev$t,ev$type),]
  A <- pv <- bound <- numeric(n); active <- rep(FALSE,n)
  rate <- rep(NA_real_,n); payment_alive <- rep(NA_real_,n)
  inputs <- outputs <- 0; last <- 0; balance <- txn <- 0
  for(e in seq_len(nrow(ev))) {
    t <- ev$t[e]; h <- t-last; idx <- which(active)
    if(h>0 && length(idx)) {
      A0 <- A[idx]; decay <- exp(-mu[idx]*h)
      pv[idx] <- pv[idx]+A0*exp(-r*last)*(-expm1(-(r+mu[idx])*h))
      # Integral of delta_i=eta*mu_i*A_i, along the observed path.
      if(is.finite(eta)) bound[idx] <- bound[idx]+eta*mu[idx]*A0*exp(-r*last)*
        (-expm1(-(r+mu[idx])*h))/(r+mu[idx])
      # Undiscounted aggregate balance: r*A-c = -mu*A.
      consumption <- (r+mu[idx])*A0*(-expm1(-mu[idx]*h))/mu[idx]
      interest <- r*A0*(-expm1(-mu[idx]*h))/mu[idx]
      outputs <- outputs+sum(consumption)-sum(interest)
      A[idx] <- A0*decay
    }
    j <- ev$id[e]; typ <- ev$type[e]
    if(typ==1L) {
      before <- A
      A[j] <- capital[j]; active[j] <- TRUE; inputs <- inputs+capital[j]
      txn <- max(txn,max(abs(A[-j]-before[-j]),0))
    } else if(typ==2L && active[j]) {
      k <- mu*A; d <- degrees(k,epsilon,eta)
      credit <- flux_col(d,j)/mu[j]
      H <- A[j]-sum(credit)
      stopifnot(H>=-1e-7)
      H <- max(0,H); pv[j] <- pv[j]+exp(-r*t)*H
      outputs <- outputs+H; A <- A+credit; A[j] <- 0; active[j] <- FALSE
    } else if(typ==3L && active[j]) {
      before <- A; outputs <- outputs+A[j]
      pv[j] <- pv[j]+exp(-r*t)*A[j]; A[j] <- 0; active[j] <- FALSE
      txn <- max(txn,max(abs(A[-j]-before[-j]),0))
    } else if(typ==4L) {
      # Alive payments include zeros after exit. Active-and-alive is distinct.
      payment_alive[j] <- if(death[j]>t) (r+mu[j])*A[j] else NA_real_
      if(active[j]) rate[j] <- (r+mu[j])*A[j]
    } else if(typ==5L) {
      pv <- pv+exp(-r*T)*A; outputs <- outputs+sum(A)
      A[] <- 0; active[] <- FALSE
    }
    balance <- max(balance,abs(sum(A)+outputs-inputs)/max(1,inputs))
    last <- t
  }
  # Report all individual benefits relative to own entry date.
  list(ratio=pv*exp(r*entry)/capital,
       bound=bound*exp(r*entry)/capital,
       rate=rate/((r+mu)*capital),
       alive_rate=payment_alive/((r+mu)*capital),balance=balance,txn=txn)
}

summ <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  c(mean=if(n) mean(x) else NA_real_,se=if(n>1) sd(x)/sqrt(n) else NA_real_,n=n)
}
states <- check_states()
write.csv(states,file.path(out,"state_checks.csv"),row.names=FALSE)
ages <- c(55,70,85)
mu <- 0.008*exp(0.09*(ages-65))
closed <- data.frame(age=ages,mu=mu,capital=c(10000,100000,1000000),entry=0)
open <- data.frame(age=rep(ages,4),mu=rep(mu,4),
                   capital=rep(c(10000,100000,1000000),4),entry=rep(c(0,2,4,6),each=3))
single <- data.frame(age=70,mu=mu[2],capital=100000,entry=0)
scenarios <- list(
  closed_nominal=list(pop=closed,z=rep(1,3),eps=0,eta=Inf,exit=0),
  closed_robust=list(pop=closed,z=rep(1,3),eps=.2,eta=.1,exit=0),
  open_exit_nominal=list(pop=open,z=rep(1,12),eps=0,eta=Inf,exit=.08),
  open_exit_robust_error=list(pop=open,z=rep(c(.8,.8,1.2),4),eps=.2,eta=.1,exit=.08),
  singleton=list(pop=single,z=1,eps=0,eta=Inf,exit=0))
individual <- aggregate <- list(); l <- 0L
for(name in names(scenarios)) {
  sc <- scenarios[[name]]; n <- nrow(sc$pop)
  ratios <- bounds <- rates <- alive_rates <- matrix(NA_real_,M,n)
  maxbal <- maxtxn <- 0
  for(m in seq_len(M)) {
    p <- path(sc$pop,sc$z,sc$eps,sc$eta,sc$exit)
    ratios[m,] <- p$ratio; bounds[m,] <- p$bound
    rates[m,] <- p$rate; alive_rates[m,] <- p$alive_rate
    maxbal <- max(maxbal,p$balance); maxtxn <- max(maxtxn,p$txn)
  }
  for(i in seq_len(n)) {
    s <- summ(ratios[,i]); q <- summ(rates[,i]); a <- summ(alive_rates[,i])
    l <- l+1L
    individual[[l]] <- data.frame(scenario=name,id=i,age=sc$pop$age[i],
      capital=sc$pop$capital[i],entry=sc$pop$entry[i],z=sc$z[i],
      epv_ratio=s[1],se=s[2],lower=s[1]-1.96*s[2],upper=s[1]+1.96*s[2],
      expected_bound=if(is.finite(sc$eta)) mean(bounds[,i]) else NA_real_,
      active_payment_ratio=q[1],active_payment_se=q[2],active_n=q[3],
      alive_payment_ratio=a[1],alive_payment_se=a[2],alive_n=a[3])
  }
  agg <- ratios%*%sc$pop$capital/sum(sc$pop$capital)
  # Entrant-date EPVs cannot be aggregated to a time-zero funding identity.
  # Convert them back to time-zero EPVs before using a capital-weighted total.
  agg <- ratios%*%(sc$pop$capital*exp(-.02*sc$pop$entry))/
    sum(sc$pop$capital*exp(-.02*sc$pop$entry))
  a <- summ(agg)
  aggregate[[name]] <- data.frame(scenario=name,M=M,epv_ratio=a[1],se=a[2],
    lower=a[1]-1.96*a[2],upper=a[1]+1.96*a[2],max_relative_balance=maxbal,
    max_other_account_change=maxtxn)
}
ind <- do.call(rbind,individual); agg <- do.call(rbind,aggregate)
write.csv(ind,file.path(out,"individual.csv"),row.names=FALSE)
write.csv(agg,file.path(out,"aggregate.csv"),row.names=FALSE)
writeLines(c(paste("Seed: 20261002; replications:",M),capture.output(sessionInfo())),
           file.path(out,"sessionInfo.txt"))
print(agg,row.names=FALSE)
print(ind[ind$scenario%in%c("closed_nominal","closed_robust","singleton"),
          c("scenario","id","epv_ratio","se","alive_payment_ratio","alive_payment_se")],row.names=FALSE)
cat("State-check maximum relative error:",max(states$row_error,states$sym_error,states$frontier_error),"\n")
cat("Singleton analytic conditional payment ratio at year 5:",exp(-single$mu*5),"\n")
