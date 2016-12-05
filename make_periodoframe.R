####This file is an example for making PeriodoFrame: the computation part
source('prepare_data.R',local=TRUE)
############################################
#####PartIV: make PeriodoFrame
############################################
#######make BFP 
t1 <- proc.time()
bf <- BFP(tab[,1],tab[,2],tab[,3],Nma=Nma,NI=NI,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmax=fmax,tol=tol)
t2 <- proc.time()
dur1 <- format((t2-t1)[3],digit=3)
cat('BFP computation time:',dur1,'s\n\n')

#######make MLP 
t1 <- proc.time()
ml <- MLP(tab[,1],tab[,2],tab[,3],Nma=Nma,NI=NI,ofac=ofac,mar.type='part',model.type=model.type,fmax=fmax,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
t2 <- proc.time()
dur2 <- format((t2-t1)[3],digit=3)
cat('MLP computation time:',dur2,'s\n')
Nma <- ml$Nma
NI <- ml$NI
fname <- paste0('periodograms_Ndata_',Ndata,'modeltype',model.type,'_',star,'_NI',NI,'Nma',Nma,'_Nc',Nc,'_ofac',ofac,'_opt',opt.type,'_MLP',MLP.type,'_Nap',Nap)
fname <- paste0('results/',fname)

################################################
#####PartV: save the data and plot periodograms
################################################
###save data
obj.name <- paste0(fname,'.Robj')
cat(obj.name,'\n')
save(list = ls(all.names = TRUE),file=obj.name)

######plot
leg.pos <- 'topright'
pdf.name <- paste0(fname,'.pdf')
pdf(pdf.name,8,8)
size <- 1.0
par(mfrow=c(2,2),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
####BFP
ylim <- c(median(bf$logBF),max(bf$logBF))
plot(bf$P,bf$logBF,xlab='Period[d]',ylab='log(BF)',type='l',log='x',ylim=ylim)
if(np==0){
    if(length(bf$Popt)>1){
#        inds <- 1
        inds <- 1:2
    }else{
        inds <- 1
    }
    ps <- bf$Popt[inds]
    ks <- 2*rev(inds)
    if(length(inds)>1){
        text(x=bf$Popt[inds],y=bf$logBF.opt[inds],pos=c(2,4),labels=paste0(format(bf$Popt[inds],digit=3),'d'),col='red')
    }else{
        text(x=bf$Popt[inds],y=bf$logBF.opt[inds],pos=4,offset=-0.1,labels=paste0(format(bf$Popt[inds],digit=3),'d'),col='red')
    }
}
abline(v=ps,col='red',lty=3,lwd=ks)
abline(h=log(150),lty=2)
legend(leg.pos,legend=c('BFP',paste0(dur1,'s')),bty='n')

####MLP 
plot(ml$P,ml$logBF-max(ml$logBF),xlab='Period[d]',ylab='Relative probability',type='l',log='x')#expression('log(ML/M'*L[max]*')')
abline(v=ps,col='red',lty=3,lwd=ks)
legend(leg.pos,legend=c('MLP',paste0(dur2,'s')),bty='n')

#####BGLS
t1 <- proc.time()
bg <- bgls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmax=fmax)
t2 <- proc.time()
dur3 <- format((t2-t1)[3],digit=3)
plot(bg$P,bg$logp,xlab='Period[d]',ylab='Relative probability',type='l',log='x')
abline(v=ps,col='red',lty=3,lwd=ks)
legend(leg.pos,legend=c('BGLS',paste0(dur3,'s')),bty='n')

#####GLS
t1 <- proc.time()
g <- gls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmax=fmax)
t2 <- proc.time()
plot(g$P,g$logp,xlab='Period[d]',ylab='Power',type='l',log='x')
abline(h=g$sig.level[1],lty=2)
abline(h=g$sig.level[2],lty=3)
abline(h=g$sig.level[3],lty=4)
abline(v=ps,col='red',lty=3,lwd=ks)
legend(leg.pos,legend=c('GLS',paste0(dur3,'s')),bty='n')

####output pdf
cat('output pdf:\n')
cat(pdf.name,'\n')
dev.off()

