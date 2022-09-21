model{
	
	for(i in 1:nobs){
		r[i]~dbin(p[i],n[i])
		p[i] <- ilogit(mu[s[i]]+delta[s[i],t[i]])
		
		rhat[i] <- p[i]*n[i]
		dev[i] <- 2*(r[i]*(log(r[i])-log(rhat[i]))+(n[i]-r[i])*(log(n[i]-r[i])-log(n[i]-rhat[i])))
		
		std[i] <- sqrt(n[i]*p[i]*(1-p[i]))
		pres[i] <- (r[i]-rhat[i])/std[i]
		
		sign[i] <- (r[i]-rhat[i])/abs(r[i]-rhat[i])
		sqdev[i] <- sqrt(abs(dev[i]))
		res[i] <- sign[i]*sqdev[i]
		
		delta[s[i],t[i]] ~ dnorm(md[i],tau)
#		delta[s[i],t[i]] ~ dnorm(md[i],taud[i])
		md[i] <- d[t[i]]-d[b[s[i]]]

# multi-arm	trial correction
#		taud[i] <- ifelse(na[i]==1,1000000,tau*2*(na[i]-1)/na[i])
		
#FE		delta[s[i],t[i]] <- ifelse(t[i]==b[s[i]],0,(d[t[i]]-d[b[s[i]]]))
	}
	for(i in 1:ns){
		mu[i]~dnorm(0,.0001)
		resdev[i] <- sum(dev[i])
		AbsTrEf[i] <- ifelse(b[i]==1,mu[i],0)
	}
	totresdev <- sum(resdev[])
	
	pi0 <- sum(AbsTrEf[])/incb
	
	rk <- nt-rank(pi[])
	
	tau <- pow(sd,-2)
	sd~dunif(0.00001,2)
	
	d[1] <- 0
	for(k in 2:nt){
		d[k]~dnorm(0,.0001)
	}
	
	for(j in 1:nt){
		logit(pi[j]) <- pi0+d[j]
		
		Best[j] <- equals(rk[j],1)
		worst[j] <- equals(rk[j],nt)
		
		for(k in 1:nt){
			lor[j,k] <- d[j]-d[k]
			log(or[j,k]) <- lor[j,k]
			rr[j,k] <- pi[j]/pi[k]
		}
	}
}