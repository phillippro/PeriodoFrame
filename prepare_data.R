source('periodoframe.R')#for PeriodoFrame
source('periodograms.R')#for other periodograms
############################################
#####PartI: set parameters
############################################
#star <- 'HD177565'
star <- 'HD41248'
np <- 0
NI <- 3
Nma <- 1
opt.type <- 'sl'#nl(nonlinear fitting all parameters; i.e. without using the formula ) or sl(semi-linear fitting)
model.type <- 'man'#manually setting the number of MA components and differential RVs
#model.type <- 'auto'#automatically determine the optimal noise model
tol <- 1e-12#numerical fitting precision tolerance 
Nap <- 0#number of aperture 
Nc <- 1#1

###To select the numbers of differential RVs and MA components, it is better not to include calibration data. 
if(model.type=='auto'){
    Nc <- 0
}
NI0 <- NI#the number of activity indices
####include all other noise proxies as indices
if(Nap>0){
    NI <- NI+Nap-1+Nc
    NI0 <- NI0+Nc
}
leg.pos <- 'topright'
MLP.type <- 'sub'#subtract the correlated noise from the data and then marginalize the non-correlated-noise parameters
#MLP.type <- 'assign'#fix the correlated noise parameters at their optimal values and then marginalize the other parameters
############################################
#####PartII: load data
############################################
f <- paste0('data/',star,'_TERRA_1AP1.dat')
tmp <- read.table(f)
if(NI>0){
    tab <- tmp[,c(1:3,4:(3+NI0-Nc))]
}else{
    tab <- tmp[,c(1:3)]
}
###load calibration data
if(Nc>0){
    tmp <- read.table(paste0('data/',star,'_3AP2-1_calibration.dat'))
    tab <- cbind(tab,tmp[,2])
}
####load differential RVs
if(model.type=='auto'){
    Naps <- c(3,6)
    for(Nap in Naps){
        tmp <- read.table(paste0('data/',star,'_TERRA_',Nap,'AP_dRVs.dat'))
        tab <- cbind(tab,tmp[,2:ncol(tmp)])
    }
}else if(Nap>0){
        tmp <- read.table(paste0('data/',star,'_TERRA_',Nap,'AP_dRVs.dat'))
        tab <- cbind(tab,tmp[,2:ncol(tmp)])
}
###sort the data in case the times are not sorted
ind.sort <- sort(tab[,1],index.return=TRUE)$ix
tab <- tab[ind.sort,]

############################################
#####PartIII: prepare data
############################################
########################
####set up
########################
#update <- FALSE
#opt.types <- c('full','part','wr','rw','rep')
ofac <- 1
fmax <- 1#Pmin=0.5d
mar.type <- 'part'
Ndata <- nrow(tab)
#model.type <- 'MA'#'MA' or 'auto'
cat('file:',f,';Ndata=',Ndata,'; Nma=',Nma,'; NI=',NI,'; Nap=',Nap,'; Nc=',Nc,'; opt.type=',opt.type,'; model.type=',model.type,'; ofac=',ofac,'; tol=',tol,'\n')

######rescale indices
Indices <- NULL
sortInd <- FALSE
if(ncol(tab)>3){
    Indices <- as.matrix(tab[,4:ncol(tab)])
    if(!is.matrix(Indices) & !is.data.frame(Indices)){
        Indices <- matrix(Indices,ncol=1)
    }
    for(j in 1:ncol(Indices)){
        Indices[,j] <- scale(Indices[,j])
    }
}
t <- tab[,1]
y <- tab[,2]
dy <- tab[,3]
data <- cbind(t,y,dy)
