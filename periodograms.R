###The following function is adapted from bgls.py by Mortier et al. 2015
bgls <- function(t, y, err, ofac=1,fmax=1,fmin=NA,tspan=NULL){
    t <- t-min(t)
    if(is.null(tspan)){
        tspan = max(t)-min(t)
    }
    if(tspan<100){
       ofac <- 10
    }
    step <- 1/(tspan*ofac)
#    fmax <- 1
    if(is.na(fmin)){
        fmin <- 1/(ofac*tspan)
    }
    f = seq(fmin,fmax,by=step)
    nout <- length(f)
    omegas <- 2*pi*f
    err2 <- err * err
    w = 1./err2
    W <- sum(w)
    bigY <- sum(w*y)  # Eq. (10)
    p <- c()
    constants <- rep(NA,length(omegas))
    exponents <- rep(NA,length(omegas))
    for(k in 1:length(omegas)){
        omega <- omegas[k]
        theta = 0.5*atan(sum(w*sin(2*omega*t))/sum(w*cos(2*omega*t)))
        x = omega*t - theta
        cosx = cos(x)
        sinx = sin(x)
	 
        wcosx = w*cosx
        wsinx = w*sinx

        C = sum(wcosx)
        S = sum(wsinx)
        YCh = sum(y*wcosx)
        YSh = sum(y*wsinx)
        CCh = sum(wcosx*cosx)
        SSh = sum(wsinx*sinx)
        if(CCh != 0 & SSh != 0){
            K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)
            L = (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
            M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2.*CCh*SSh)
            constants[k] <- 1./sqrt(CCh*SSh*abs(K))
        }else if(CCh == 0){
            K = (S*S - W*SSh)/(2.*SSh)
            L = (bigY*SSh - S*YSh)/(SSh)
            M = (YSh*YSh)/(2.*SSh)
            constants[k] <- 1./sqrt(SSh*abs(K))
        }else if(SSh == 0){
            K = (C*C - W*CCh)/(2.*CCh)
            L = (bigY*CCh - C*YCh)/(CCh)
            M = (YCh*YCh)/(2.*CCh)
            constants[k] <- 1./sqrt(CCh*abs(K))
        }
        if(K > 0){
            cat('K is positive. This should not happen.')
        }
        exponents[k] <- M - L^2/(4*K)
    }
    logp  <- log10(constants) + (exponents*log10(exp(1)))
#    p <- 10^logp
    logp <- (logp-max(logp))*log(10)  # normalize and with a base of natural constant
    return(list(P=1/f, logp=logp,p=exp(logp)))
}

###based on Zechmeister09.pdf or ZK09 and the lsp function in the 'lomb' library
gls <- function(t, y, err,ofac=1, norm="Cumming",fmax=1,fmin=NA,tspan=NULL){
    t <- t-min(t)
    if(is.null(tspan)){
        tspan = max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
#    f.max <- 1
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    f = seq(fmin,fmax,by=step)
    nout <- length(f)
#    nout = ofac * hifac * length(t)/2
#    xdif = max(t)-min(t)
#    f = 1./(xdif*ofac) + c(1:nout)/(ofac*xdif)
    omegas <- 2*pi*f
    err2 <- err * err
    W <- sum(1/err2)
    w <- 1/err2/W
    bigY <- sum(w*y)  # Eq. (7)
    p <- rep(NA,length(omegas))
    for(i in 1:length(omegas)){
        omega <- omegas[i]
        bigC <- sum(w*cos(omega*t))
        bigS <- sum(w*sin(omega*t))
        YY.hat <- sum(w*(y^2))
        YC.hat <- sum(w*y*cos(omega*t))
        YS.hat <- sum(w*y*sin(omega*t))
        CC.hat <- sum(w*cos(omega*t)^2)
        SS.hat <- sum(w*sin(omega*t)^2)
        CS.hat <- sum(w*sin(omega*t)*cos(omega*t))
        YY <- YY.hat-bigY*bigY
        YC <- YC.hat-bigY*bigC
        YS <- YS.hat -bigY*bigS
        CC <- CC.hat - bigC*bigC
        SS <- SS.hat - bigS*bigS
        CS <- CS.hat - bigC*bigS
        bigD <- CC*SS-CS^2
        p[i] <- (SS*YC^2+CC*YS^2-2*CS*YC*YS)/(YY*bigD)
    }
    N <- length(y)
 # An ad-hoc estimate of the number of independent frequencies (ZK_09 Eq. 24)
    #M <- (max(f)-min(f))*(max(t)-min(t))
    M <- 2*nout/ofac#ref lsp
    if(norm=='Scargle'){
        popvar <- 1#arbitrary value or input
        power <- p/popvar
    }
    if(norm=='HorneBaliunas'){
        power <- (N-1)*p/2
    }
    if(norm=='Cumming'){
        power <- ((N-3)/2)*p/(1-max(p))
    }
    PN <- power
    PN.max <- max(PN)
    peak.freq <- f[PN==PN.max]
    peak.per <- 1/f[PN==PN.max]
    FAP <- c(0.1,1e-2,1e-3)#significance level of FAP
    level <- powerLevel(FAP,M,N,norm)#power level
    pp <- M*prob(Pn=PN.max,N=N,norm=norm)
    if(pp>0.01){
        pp <- 1-(1-prob(Pn=PN.max,N=N,norm=norm))^M
    }
#    cat('level:',level,'\n')
    return(list(P=1/f, logp=power, pvalue=pp, sig.level=level))
}

###generalized lomb-scargle periodogram with trend component
glst <- function(t, y, err,ofac=1, norm="Cumming",fmax=1,fmin=NA,tspan=NULL){
    unit <- 365.24#to make the elements of the matrix in the function of 'solve' on the same order
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    f <- seq(fmin,fmax,by=step)*unit
    nout <- length(f)
    t <- (t-min(t))/unit
    omegas <- 2*pi*f
    err2 <- err * err
    W <- sum(1./err2)
    w <- 1/err2/W
    bigY <- sum(w*y)  # Eq. (7)
    p <- rep(NA,length(omegas))
    bigT <- sum(w*t)
    YY.hat <- sum(w*y^2)
    YT.hat <- sum(w*y*t)
    TT.hat <- sum(w*t^2)
    YY <- YY.hat - bigY*bigY
    YT <- YT.hat - bigY*bigT
    TT <- TT.hat - bigT*bigT
###optimized parameterse for the trend model
    d0 <- YT/TT
    c0 <- bigY-d0*bigT
    chi2.ref <- sum(W*w*(y-(c0+d0*t))^2)
    for(k in 1:length(f)){
        omega <- omegas[k]
        bigC <- sum(w*cos(omega*t))
        bigS <- sum(w*sin(omega*t))
        YY.hat <- sum(w*y^2)
        YC.hat <- sum(w*y*cos(omega*t))
        YS.hat <- sum(w*y*sin(omega*t))
        CC.hat <- sum(w*cos(omega*t)^2)
        SS.hat <- sum(w*sin(omega*t)^2)
        CS.hat <- sum(w*sin(omega*t)*cos(omega*t))
        ST.hat <- sum(w*sin(omega*t)*t)
        CT.hat <- sum(w*cos(omega*t)*t)
        YC <- YC.hat - bigY*bigC
        YS <- YS.hat - bigY*bigS
        CC <- CC.hat - bigC*bigC
        SS <- SS.hat - bigS*bigS
        ST <- ST.hat - bigS*bigT
        CS <- CS.hat - bigC*bigS
        bigD <- CC*SS-CS^2
        lin.mat <- matrix(c(CC.hat,CS.hat,bigC,CT.hat,CS.hat,SS.hat,bigS,ST.hat,bigC,bigS,1,bigT,CT.hat,ST.hat,bigT,TT.hat),byrow=TRUE,nrow=4)
        vec.rh <- c(YC.hat,YS.hat,bigY,YT.hat)
        pp <- solve(lin.mat,vec.rh,tol=1e-16)
        yp <- pp[1]*cos(omega*t)+pp[2]*sin(omega*t)+pp[3]+pp[4]*t
        chi2 <- sum(W*w*(y-yp)^2)
        p[k] <- (chi2.ref-chi2)/chi2.ref
    }
    N <- length(y)
    M <- 2*nout/ofac#ref lsp
    if(norm=='Scargle'){
        popvar <- 1#arbitrary value or input
        power <- p/popvar
    }
    if(norm=='HorneBaliunas'){
        power <- (N-1)*p/2
    }
    if(norm=='Cumming'){
        power <- ((N-3)/2)*p/(1-max(p))
    }
    PN <- power
    PN.max <- max(PN)
    peak.freq <- f[PN==PN.max]
    peak.per <- 1/f[PN==PN.max]
#    FAP <- c(0.317,0.046,3e-3)#significance level of FAP
    FAP <- c(0.1,0.01,0.001)
    level <- powerLevel(FAP,M,N,norm)#power level
    pp <- M*prob(Pn=PN.max,N=N,norm=norm)
    if(pp>0.01){
        pp <- 1-(1-prob(Pn=PN.max,N=N,norm=norm))^M
    }
    P <- unit/f
    ind.max <- which.max(power)
    return(list(P=unit/f, logp=power, pvalue=pp, sig.level=level,Popt=P[ind.max]))
}

#give a power, calcuate the the p value
prob <- function(Pn,N,norm='Cumming'){
    if(norm=="Scargle") return(exp(-Pn))
    if(norm=="HorneBaliunas") return((1-2*Pn/(N-1))^((N-3)/2))
    if(norm=="Cumming") return((1+2*Pn/(N-3))^(-(N-3)/2))
}
###Inverse of prob
probInv <- function(Prob,N,norm='Cumming'){
    if(norm=="Scargle") return(-log(Prob))
    if(norm=="HorneBaliunas") return((N-1)/2*(1-Prob^(2/(N-3))))
    if(norm=="Cumming") return((N-3)/2*(Prob^(-2/(N-3))-1))
}
####
powerLevel <- function(FAPlevel,M,N,norm='Cumming'){
    return(probInv(1-(1-FAPlevel)^(1/M),N,norm))
}
#####
tv.ls <- function(t, y, err,Dt,nbin,fmax=1,ofac=1,fmin=1/1000,tspan=NULL,per.type='bgls'){
    n <- nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    tmid <- (tstart+tend)/2
    rel.powers <- c()
    powers <- c()
    levels <- c()
    ndata <- c()
    for(j in 0:n){
        inds <- which(t>=min(t)+j*dt & t<min(t)+j*dt+Dt)
#        tmp <- glst(t=t[inds],y=y[inds],err=err[inds],fmax=fmax,ofac=ofac,fmin=1/Dt,tspan=Dt)
        if(per.type=='bgls'){
            tmp <- bgls(t=t[inds],y=y[inds],err=err[inds],fmax=fmax,ofac=ofac,fmin=1/Dt,tspan=Dt)
        }else if(per.type=='gls'){
            tmp <- gls(t=t[inds],y=y[inds],err=err[inds],fmax=fmax,ofac=ofac,fmin=1/Dt,tspan=Dt)
        }else if(per.type=='glst'){
            tmp <- glst(t=t[inds],y=y[inds],err=err[inds],fmax=fmax,ofac=ofac,fmin=1/Dt,tspan=Dt)
        }
        index <- sort(tmp$P,index.return=TRUE)$ix
        ##scaling the power
        dlevel <- tmp$sig.level[2]-tmp$sig.level[1]
#        rel.power <- (tmp$logp[index]-tmp$sig.level[1])/dlevel
#        rel.power <- (tmp$logp[index]-min(tmp$logp))#/length(inds)
        rel.power <- scale(tmp$logp[index])#/length(inds)
        powers <- cbind(powers,tmp$logp[index])
        rel.powers <- cbind(rel.powers,rel.power)
        levels <- cbind(levels,tmp$sig.level)
        ndata <- c(ndata,length(inds))
    }
    return(list(tmid=tmid,ps=tmp$P[index],powers=powers,rel.powers=rel.powers,levels=levels,ndata=ndata))
}
