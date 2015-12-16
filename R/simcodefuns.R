library(VGAM)
require(stats)


calcBaf<-function(str){
  out <- strsplit(as.character(str),",")[[1]]
  naN <- as.numeric(out[1])
  nbN <- as.numeric(out[2])
  bafN <- nbN/(nbN+naN)
  c(naN,nbN,bafN,nbN+naN)
}

binomP<-function(p,x,size) dbinom(x,size=size,prob=p,log=TRUE)

sampBeta<-function(vN,iterN,rep,val) {
 if(!is.null(val)) return(rep(val,iterN))
 else if(rep) return(rep(vN[2]/(vN[1]+vN[2]),iterN))
 else return (rbeta(iterN, vN[2]+1, vN[1]+1))
}

#[normal,tumour,mixture]
probF<-function(v,c,iter,rep,val){
   sampsN = sampBeta(v[1,],iter[1],rep[1],val[1]) 
   sampsT = sampBeta(v[2,],iter[2],rep[2],val[2]) 
   prob = outer(c,sampsT) + outer( (1-c),sampsN)
   Lmatr = apply(prob,c(1,2),binomP, v[3,2], v[3,4])
   LL = apply(Lmatr,1,sum)/iter
   LL
}

#[normal,tumour,mixture]

calcLL<-function(subdatk,c, thresh,iter, rep,val,multF = c(1,1,1)){
   inds = c(5,9,13)  #indices of counts for tumour, normal and mixture respectively
   v = array(dim = c(length(inds),4))
   for(k in 1:length(inds)) {
	v[k,] = calcBaf(subdatk[inds[k]])
       v[k,c(1,2,4)] = v[k,c(1,2,4)]*multF[k]
   }
   
   if(v[1,4]>=thresh[1] & v[2,4]>=thresh[2])
   {
     return(probF(v,c,iter,rep,val))
    #prob = c*vT[3] + (1-c)*vN[3]
    #LL = dbinom(x = vM[2], size = vM[1]+vM[2], prob = prob,log=TRUE)
   }else return(rep(0,length(c)) )
}

#e.g. set multF = c(1,0,1) in order to investigate having no tumour data
calcLLsum<-function(subdat,Cset,thresh, iter, rep,val, multF = c(1,1,1), maxlen = dim(subdat)[1]){
LL = apply(subdat[1:maxlen,],1,calcLL,Cset,thresh,iter,rep,val,multF=multF)
dimnames(LL)[[1]] = Cset
LLsum = apply(LL,1,sum)
###LLsum = LLsum - LLsum[maxi]
plot(Cset, LLsum, xlab="cellularity", ylab="log likelihood", type="b", pch=16)
LLsum
}

##
makeInitialCNDist<-function(max=4,rate=4,by = 1){
 x = seq(0,max,by = by)
 res = cbind(x/2,(dgamma(x,shape=2*rate+1,rate=rate)))
 res[1,2] = res[5,2]
 res[,2] = res[,2]/sum(res[,2])
 res
}


###R = [R_M, R_N]; S = [S_M,S_N]
binomPvec<-function(R_,S,paramsCN,mix,prop=1.0){
 CN = paramsCN[,1]
 pCN = paramsCN[,2]
 
 prob = (R_[2]/S[2])*((CN-1)*mix+1)
 #R2 = prob*S[2]
vec = dbinom(R_[1],size=S[1],prob=prob,log=TRUE) + log(pCN)
# vec = dbetabinom.ab(R_[1],size = S[1],shape1 = prop*(R2+1),shape2 = prop*(S[2]-R2+1),log=TRUE)+log(pCN)
 m = max(vec)
 pr = exp(vec-m)
 sp = sum(pr)
 m1 = m+log(sp)
 c(m1,exp(vec-m1))
}

#Pi = [Pi_M,Pi_N]
calcC<-function(v,CN){
 res =  (v[1]/v[2] -1)/(CN-1)
}

##thresh is the threshold
##order is R_M, R_N and S_M, S_N
#  Wenhan modified
calcSigLeft<-function(RS, Pi, thresh = 0,log=T,prop = 1/3,minLogP = -750, c = 0, RCN = -1){
 R_  = RS[1:2]
 S = RS[3:4]
 x = R_[1] 

 prob = R_[2]/S[2]

 #                                Nfb              Nfb * RCN 
 #     P = f(c, RCN) = (1 - c) ───────────  + c ────────────────────
 #                              Nfa + Nfb        Nfa + Nfb * RCN 
# if ( length(RCN) == 1 ||RCN != -1)
# {
#     prob = (1 - c) * prob + c * R_[2] *  RCN / (S[2] - R_[2] + R_[2]*RCN)
# }
 p = pbinom(x,size=S[1],prob=Pi,log=log)  ##P (X< =x) 
 d = dbinom(x,size=S[1],prob=Pi,log=log)  ##P (X< =x) 

 if(log) return(min(0,p+log(2))) ####need to fix this
 else return(min(1, 2*p-d))

}

calcSigRight<-function(RS, Pi,thresh = 0,log=T,prop = 1/3,minLogP = -29, c = 0, RCN = -1){
 R_ = RS[1:2]
 S = RS[3:4]
 x = R_[1]-1  

 prob = R_[2]/S[2]

#                                Nfb              Nfb * RCN
# prob = f(c, RCN) = (1 - c) ───────────  + c ────────────────────
#                              Nfa + Nfb        Nfa  + Nfb * RCN
#
# if (length(RCN) == 1 ||RCN != -1)
# {
#     prob = (1 - c) * prob + c * R_[2] * RCN/ ((S[2] - R_[2])  + R_[2] * RCN)
# }

 p = pbinom(x,size=S[1],prob=Pi,lower.tail=F,log=log)
 d = dbinom(x+1,size=S[1],prob=Pi,log=log)

 

 if(log) return(min(0,p+log(2))) ###need to think about this
 else return(min(1, 2*p - d))
}

calcSigAll1<-function(Mk,thresh=0,log=T,prop=1/3){
  calcSigAll(Mk$Pi,Mk$R,Mk$S, thresh=thresh,log=log,prop=prop)
}

##first calculates pbinom, then if pv<thresh (or log(pv)<thresh if logged) then uses betabinom
##order is [R_M, R_N] and [Pi_M, P_N]
# Wenhan modified
calcSigAll<-function(Pi, R, S, thresh = 0,log = T,prop = 1/3, tt = 1, c = 0, RCN = -1) {

  data=cbind(R,S) 

  pleft=c()
  pright=c()

  if (tt == 1)
  {
    leftinds  = (Pi[,1] < Pi[,2])
    rightinds = (Pi[,1] >= Pi[,2])
    left  = data[leftinds, ,drop=F]
    right = data[rightinds,,drop=F]
    p_l = Pi[ leftinds,,drop=F]
    p_r = Pi[rightinds,,drop=F]
    pleft  = rep(NA, length=dim(p_l)[1])
    pright = rep(NA, length=dim(p_r)[1])
    if(dim(left)[1] != 0)
    {
        for (xx in 1:dim(p_l)[1])
        {
            pleft[xx] = calcSigLeft(left[xx,], p_l[xx, 2],thresh = thresh,log=log,
                                    prop=prop, c = c, RCN = RCN)
        }
      #pleft =   apply(matrix(left,ncol=ncol(data)),1,calcSigLeft, thresh = thresh,log=log,prop=prop, c = c, RCN = RCN)
    }
    if(dim(right)[1] != 0)
    {
        for (xx in 1:dim(p_r)[1])
        {
            pright[xx] = calcSigRight(right[xx,], p_r[xx, 2], thresh = thresh, log=log,
                                    prop=prop, c = c, RCN = RCN)
        }
       #pright =   apply(matrix(right,ncol=ncol(data)),1,calcSigRight,thresh = thresh,log=log,prop=prop, c = c, RCN = RCN)
    }
    res = rep(NA,length(pleft)+length(pright))
    res[leftinds] = pleft
    res[rightinds] = pright
    list(left = pleft,right = pright,all= res)
  } else if (tt == 2)
  {  

    res = rep  (NA,length(R))


    data=array(R[,1], c(nrow(R), 2,2))


    dim(data)
    print (length(S[,1])) 
    print (length(R[,1])) 

    data[,1,2] = S[,1] - R[,1]
    data[,2,1] = R[,2]
    data[,2,2] = S[,2] - R[,2]



    res = apply(data, 1, fisher.test)
    #res = apply(cbind(R, S - R),1, chisq.test)

    pvSet=rep(0, length(res))

    for (ii in 1:length(res)){
        pvSet[ii]=res[[ii]]$p.value
    }



#    res = apply(cbind(R, S - R),1, chi)


    #pvSet = res

    list(left = NA,right = NA,all= pvSet)
    
  } else
  {
      print("Provide tt as the test type.")
  }

}

chi <- function (freqs)
{
    x <- matrix(freqs, ncol = 2)
    
    p = chisq.test(x)$p.value           # 0.4233
    return(p)
    #   ----- if using Monte Carlo -----
    # chisq.test(x, simulate.p.value = TRUE, B = 10000)$p.value
}

calcSigAllBack<-function(Pi,R,S, thresh= 0,log=T,prop=1/3) {
  leftinds = which(Pi[,1] < Pi[,2])
  rightinds = which(Pi[,1] >= Pi[,2])
  Rleft  = R[leftinds,]
  Rright = R[rightinds,]
  pleft  = apply(Rleft,1,calcSigLeft,S, thresh = thresh,log=log,prop=prop)
  pright = apply(Rright,1,calcSigRight,S,thresh = thresh,log=log,prop=prop)
  res = rep(NA,length(pleft)+length(pright))
  res[leftinds] = pleft
  res[rightinds] = pright
  list(left = pleft,right = pright,all= res)
}

##order is [Pi_M,Pi_N]
objF<-function(mix1,Ep,Pi,CN,S,priorweight = 1e7){
 expCNT = (Ep%*%CN)
 Pi_T = (Pi[,2] * expCNT)
 Pi_M = Pi_T * mix1 + Pi[,2] * (1-mix1)
# R_ = cbind(Pi[,1]*S[1],Pi_M*S[2])
# sig = calcSigAll(cbind(Pi[,1],Pi_M),R_,S,thresh=1,log=T,prop=1/3)
# res =  -sum(sig)

 res =  sum((abs(Pi_M - Pi[,1]))*S[1]) + abs(mix1)*priorweight
 res
}

##R is a matrix with col order [R_M, R_N]
##S is a vector with S_M, S_N
##paramsCN is a matrix with [CN,p_CN]
### alpha is pseudo counts
calcEmiss<-function(Mk,paramsCN,mix,noiter = 1, alpha = 1e3,br=1,priorweight=1e7){
R = Mk$R
S = Mk$S
Pi = Mk$Pi
normInd = which(paramsCN[,1]==1)
prevLL = -Inf
for(i in 1:noiter){
 res = t(apply(R,1,binomPvec,S,paramsCN,mix))
 dimnames(res)[[2]] = c("log_sumprobs",paramsCN[,1])
 logsums =  res[,1]
 Ep = res[,-1]
  opt = optim(mix,objF,gr=NULL,Ep,Pi,paramsCN[,1],S,priorweight,method="Brent",lower=0,upper=1)
  mix = opt$par
 if(alpha < 1e8){
    s = apply(Ep,2,sum) + alpha * paramsCN[,2]
    s = s/sum(s) 
   print(round(s*100))
   paramsCN[,2] = s
 }
 LL = sum(logsums)
 print(paste("mix",mix))

 print(paste("logL",sum(logsums)))
 if(LL < prevLL){
   warning('LL should decrease')
 }
 if(LL-prevLL<br) break
 prevLL = LL
}
list(paramsCN = paramsCN, LL = LL,logsums = logsums,mix=mix,Ep=Ep,loc = Mk$loc)
}

combineCounts<-function(vec,weights){
  res =  round(vec[1:length(weights)]%*% weights)
  res1 = c(res,vec[-(1:length(weights))])
  res1
}


## if synthetic mix, then calculates as proportion of first and last
extractChroms<-function(datTCGA1,nme="P60",thresh =100,len = 1, syntheticMix = TRUE){
 indRef = which(names(datTCGA1)=="Ref")
 restot = list()
 nme1 = nme
 if(syntheticMix){
   nme1 = list() 
   for(k in 1:length(nme)) nme1[[k]] = c(nme[1],nme[length(nme)])
 }
 ind =list()
 allinds = c(indRef)
 for(k in 1:length(nme)){
  ind[[k]] = which(names(datTCGA1) %in% nme1[[k]])
  allinds = unique(c(allinds,ind[[k]]))
 }
 nainds = which(is.na(apply(datTCGA1[,allinds,drop=F],1,sum)) |  datTCGA1[,indRef]<=thresh)
 datTCGA = datTCGA1[-nainds,] 
 locs = datTCGA[,1:3]
 if(syntheticMix) R = apply(datTCGA[,c(ind[[1]],indRef)],c(1,2),as.numeric)
 for(k in 1:length(nme)){
   if(!syntheticMix) R = apply(datTCGA[,c(ind[[k]],indRef)],c(1,2),as.numeric)
   S = apply(R,2,sum)
   if(dim(R)[2]>2){
     mix = as.numeric(sub('P','',nme[k]))/100
     weights = (1/S[1:2]) * c(mix,1-mix)
     weights = weights/sum(weights)
     R1 = t(apply(R,1,combineCounts,weights))
     dimnames(R1)[[2]][1] = nme[k] 
     R = R1
     S = apply(R,2,sum)
   }

   Pi = cbind(R[,1]/S[1],R[,2]/S[2])
   dimnames(Pi)= dimnames(R)
   restot[[k]] = list(R = R,S =S,Pi = Pi, locs = locs)
}
names(restot) = nme
resall = list()
arms = levels(as.factor(locs[,1]))
for(j in 1:length(arms)){
 arminds = which(locs[,1]==arms[j])
 res = list()
 for(k in 1:length(nme)){
   R = restot[[k]]$R[arminds,]
   Pi = restot[[k]]$Pi[arminds,]
   S =restot[[k]]$S
   res[[k]] = joinNeighbours(list(R = R,S =S,Pi = Pi, locs = locs[arminds,]),len=len)
 }
 names(res) = nme
 resall[[j]] = res
}
if(len>1){  ## remake the total if chroms have been joined
   for(k in 1:length(nme)){
      R = resall[[1]][[k]]$R
      Pi = resall[[1]][[k]]$Pi
      locs = resall[[1]][[k]]$locs
      S = resall[[1]][[k]]$S
      for(j in 2:length(arms)){
	R = rbind(R,  resall[[2]][[k]]$R)
	Pi = rbind(Pi,  resall[[2]][[k]]$Pi)
	locs = rbind(locs,  resall[[2]][[k]]$locs)
      }
      restot[[k]] = list(R=R,S=S,Pi=Pi,locs=locs)
   }
}
resall[[j+1]] = restot
names(resall) = c(arms,"all")
resall
}

joinNeighbours<-function(Mk,len=10){
  if(len==1) return(Mk)
     totl = dim(Mk$locs)[[1]]
     st = seq(1,totl,by = len)
     en = st +len
     excl = which(en>totl)
     if(length(excl>0)){
       st = st[-excl] 
       en = en[-excl]
     }
     R = Mk$R[st,]
     for(k in 1:(len-1)){
       R = R + Mk$R[st+k,,drop=F]
     }
     S = Mk$S
     locs = Mk$locs[st,,drop=F]
     locs[,3] =Mk$locs[en,3]
     Pi = cbind(R[,1,drop=F]/S[1],R[,2,drop=F]/S[2])
     dimnames(Pi)= dimnames(R)
     res = list(R = R,S =S,Pi = Pi, locs = locs)
       names(res) = names(Mk)
 res
}

addChroms<-function(Mall){
 arms = names(Mall)
 nme = names(Mall[[1]])
 res = list()
 for(k in 1:length(nme)){
  nmek = dimnames(Mall[[1]][[k]]$R)[[2]]
  nmek1 = dimnames(Mall[[1]][[k]]$locs)[[2]]
  R = array(dim= c(length(arms),2),dimnames = list(arms,nmek))
  Pi = array(dim= c(length(arms),2),dimnames = list(arms,nmek))
  locs = array(dim= c(length(arms),3),dimnames = list(arms,nmek1))
  for(j in 1:length(arms)){
    R[j,] = apply(Mall[[j]][[k]]$R,2,sum)
    locs[j,1] = arms[j]
    locs[j,2:3] = c(1,Inf)
  }
  S = apply(R,2,sum)
  names(S) = nmek
  Pi = cbind(R[,1]/S[1],R[,2]/S[2])
  dimnames(Pi)= dimnames(R)
  res[[k]] = list(R = R,S =S,Pi = Pi, locs = locs)
 }
 names(res) = nme
 res
}



## just function to get number of segments per chrom
getLengths<-function(Mall){
  nme = names(Mall)
 res = array(dim = c(length(nme),2))
 for(k in 1:length(Mall)){
   len = dim(Mall[[k]]$P0$R)[1]
   res[k,] = c(nme[k],len)
 }
 res[order(as.numeric(res[,2])),]
}


##Note, M must be a list of different proportions (list can have one element though)
calcEmissAll<-function(M,Cseq,paramsCN,noiter=1,plot=F,alpha = 10, br = 1,priorweight = 1e7){
 LL = array(dim = c(length(M),length(Cseq)))
 mix = array(dim = c(length(M),length(Cseq)))
 maxl = rep(0,length(M))
 minl = rep(0,length(M))
 dimnames(LL) = list(names(M),Cseq)
 EImax = list()
 bestind = rep(0,length(M))

 for(k in 1:length(M)){
  Mk = M[[k]]
  EIs = list()
  for(i in 1:length(Cseq)){
    print(Cseq[i])
    Ei  = calcEmiss(Mk, paramsCN,Cseq[i],noiter = noiter,alpha = alpha,br = br,priorweight=priorweight)
    mix[k,i] = Ei$mix
    LL[k,i] = Ei$LL
    EIs[[i]] = Ei
  }
  maxl[k] = max(LL[k,])
  
  bestind[k] = which(LL[k,]==maxl[k])
  EImax[[k]] = EIs[[bestind[k]]]
  minl[k] = min(LL[k,])
  LL[k,] = LL[k,]
  if(plot) plot(Cseq,LL[k,])
 }
 list(LL=LL,maxl=maxl,minl = minl,ymin = min(minl), ymax = max(maxl),EImax = EImax,bestind=bestind,mix = mix)
}


##log if the values are logp
p2chisq<-function(p,log=F) qchisq(p,df =1,lower.tail=F,log.p = log)

## calcs cumulative significance statistic
cumulSig<-function(pvs,title, plot=T, ggplotEnable=T, log = TRUE, log1 = FALSE, ymax=NA)
{
    pv = if(log) exp(pvs$all) else pvs$all
    p = pchisq(sum(p2chisq(pvs$all,log = log)),length(pvs$all),lower.tail=F,log = log1)
    if(plot)
    {
        print("Start plotting ...")
        if(ggplotEnable == T)
        {
            p <- ggplotQQ(pv, conf = c(0.025,0.975), title=title, ymax=ymax)
        } else
        {
            logQQ(pv,conf = c(0.025,0.975), title=paste(title,"pv=",sprintf("%5.2g",p)))
        }
    }
    p
}

merge<-function(tab,len=100){
rep = ceiling(dim(tab)[1]/len)
res = matrix(0,nrow = rep,ncol = dim(tab)[2])
for(k in 1:rep){
  start = (k-1)*len +1
  end = min(start +len,dim(tab)[1])
  res1 = apply(tab[start:end,,drop=F],2,sum)
    tmp=tab[start:end,,drop=F]
  fluctuation[k] <<- sd(tmp$NRef/(tmp$NAlt+tmp$NRef))  # lexically scoped
  res[k,] = res1
}

dimnames(res)[[2]] = dimnames(tab)[[2]]
colnames(res) = colnames(tab)
rownames(res) = rownames(tab)[seq(1,length(rownames(tab)), len)]

names(fluctuation) <<- rownames(tab)[seq(1,length(rownames(tab)), len)]
data.frame(res)
}



mergeByDistTest2<-function(tab, cnStates, genotypes, dist=10000)
{

mySum2 <- function (alleleDMatrix, cnStates, geno)
{
    sum_NRef = alleleDMatrix$nB[1] 
    sum_NAlt = alleleDMatrix$nSum[1]
    sum_TRef = alleleDMatrix$tB[1]
    sum_TAlt = alleleDMatrix$tSum[1]

#baf     = alleleDMatrix$NRef / (alleleDMatrix$NRef + alleleDMatrix$NAlt)
    stopAdd = F
    if(length(cnStates) >= 2)
    {
        for(i in 2:length(alleleDMatrix$nB))
        {
            current = c(geno[cnStates[i] * 2 -1],   geno[cnStates[i] * 2])
            prev    = c(geno[cnStates[i-1] * 2 -1], geno[cnStates[i-1] * 2])
            if(identical(current,rev(prev)))
            {
                switchPos <<- switchPos +1                      # lexical scopying switchPos
            }
            if(switchPos == 0 & !stopAdd)
            {
                sum_NRef = alleleDMatrix$nB[i] + sum_NRef
                sum_NAlt = alleleDMatrix$nSum[i] + sum_NAlt
                sum_TRef = alleleDMatrix$tB[i] + sum_TRef
                sum_TAlt = alleleDMatrix$tSum[i] + sum_TAlt

            } else
            {
                stopAdd = T
              #  sum_NRef = alleleDMatrix$nSum[i] - alleleDMatrix$nB[i] + sum_NRef
              #  sum_NAlt = alleleDMatrix$nSum[i] + sum_NAlt
              #  sum_TRef = alleleDMatrix$tSum[i] - alleleDMatrix$tB[i] + sum_TRef
              #  sum_TAlt = alleleDMatrix$tSum[i] + sum_TAlt
            }
        }
    }
    c(sum_NRef, sum_NAlt, sum_TRef, sum_TAlt)

}  # ( function mySum2)  ends
    pos <- as.numeric(rownames(tab))
    sep = c(0)
    length(sep) <- dim(tab)[1]
    idx = 1
    sep[idx] = 1
    idx = idx +1
    for(k in seq(2,length(pos),1))
    {
        if (pos[k] - pos[sep[idx -1]] > dist)
        {
            sep[idx] = k
            idx = idx +1
        }
    }

    if(idx < length(pos))
    {
        sep[idx] = length(pos)
        idx = idx +1
    }

    res = matrix(0, nrow = idx -2, ncol = dim(tab)[2])
    rownames(res) <-  rep("",idx -2)

    #print (c(idx-1,sep[idx-1], length(pos), pos[length(pos)]))
    numDist = c()
    for(k in 2:(idx-1)){
        start = sep[k -1]
        end = sep[k]

numDist = c(numDist, nrow(tab[start:end,,drop=F]))

        if(nrow(tab[start:end,,drop=F]) > 5)
        {
        switchPos = 0
        res1 = mySum2(tab[start:end,,drop=F], cnStates[start:end], genotypes)
        breakpointMark <<- c(breakpointMark, switchPos)       # lexical scoping breakpointMark
        
        #res1 = apply(tab[start:end,,drop=F],2,sum)
        res[k-1,] = res1
        } else
        {
            res[k-1,] = rep(NA, ncol(tab[start:end,,drop=F]))
        }

        rownames(res)[k-1] = toString(pos[start])
        #print(c(start,sep[start],pos[sep[start]]))
    }
    pdf("numDist.pdf")
    hist(numDist)
    dev.off()

    colnames(res) = colnames(tab)
    data.frame(res)

}


mergeByDistTest<-function(tab,dist=10000)
{
# stops at band switching
mySum <- function (alleleDMatrix)
{
    sum_NRef = 0
    sum_NAlt = 0
    baf     = alleleDMatrix$NRef / (alleleDMatrix$NRef + alleleDMatrix$NAlt)
    prevBaf = NA

    stopAdd = F
    for(i in 1:length(alleleDMatrix$NRef))
    {
        if(i!=1)
        {
            if(!(is.na(baf[i]) | is.na(prevBaf)) &&abs(prevBaf - baf[i]) > 0.5)
            {
                switchPos <<- switchPos +1                      # lexical scopying switchPos
            } else
            {
                if(!stopAdd)
                {
                    sum_NRef = alleleDMatrix$NRef[i] + sum_NRef
                    sum_NAlt = alleleDMatrix$NAlt[i] + sum_NAlt
                }
            }
        }
        prevBaf = baf[i]
    }
    c(sum_NRef, sum_NAlt)

}



    pos <- as.numeric(rownames(tab))
    sep = c(0)
    length(sep) <- dim(tab)[1]
    idx = 1
    sep[idx] = 1
    idx = idx +1
    for(k in seq(2,length(pos),1))
    {
        if (pos[k] - pos[sep[idx -1]] > dist)
        {
            sep[idx] = k
            idx = idx +1
        }
    }

    if(idx < length(pos))
    {
        sep[idx] = length(pos)
        idx = idx +1
    }

    res = matrix(0, nrow = idx -2, ncol = dim(tab)[2])
    rownames(res) <-  rep("",idx -2)

    #print (c(idx-1,sep[idx-1], length(pos), pos[length(pos)]))
    numDist = c()
    for(k in 2:(idx-1)){
        start = sep[k -1]
        end = sep[k]
        numDist = c(numDist, nrow(tab[start:end,,drop=F]))
        if(nrow(tab[start:end,,drop=F]) > 5)
        {
            switchPos =0
            res1 = mySum(tab[start:end,,drop=F])
            
            breakpointMark <<- c(breakpointMark, switchPos)       # lexical scoping breakpointMark
            #res1 = apply(tab[start:end,,drop=F],2,sum)
            res[k-1,] = res1
        } else
        {
            res[k-1,] = rep(NA, ncol(tab[start:end,,drop=F]))
        }

        rownames(res)[k-1] = toString(pos[start])
        #print(c(start,sep[start],pos[sep[start]]))
    }
    pdf("numDist.pdf")
    hist(numDist)
    dev.off()

    colnames(res) = colnames(tab)
    data.frame(res)

}

mergeByDist<-function(tab,genotypes, dist=10000)
{

    pos <- as.numeric(rownames(tab))
    sep = c(0)
    length(sep) <- dim(tab)[1]
    idx = 1
    sep[idx] = 1
    idx = idx +1
    for(k in seq(2,length(pos),1))
    {
        if (pos[k] - pos[sep[idx -1]] > dist)
        {
            sep[idx] = k
            idx = idx +1
        }
    }

    if(idx < length(pos))
    {
        sep[idx] = length(pos)
        idx = idx +1
    }

    res = matrix(0, nrow = idx -2, ncol = dim(tab)[2])
    rownames(res) <-  rep("",idx -2)

    #print (c(idx-1,sep[idx-1], length(pos), pos[length(pos)]))
    numDist = c()
    for(k in 2:(idx-1)){
        start = sep[k -1]
        end = sep[k]

        res1 = apply(tab[start:end,,drop=F],2,sum)
        res[k-1,] = res1
        rownames(res)[k-1] = toString(pos[start])
        #print(c(start,sep[start],pos[sep[start]]))
    }
    colnames(res) = colnames(tab)
    data.frame(res)

}







multiArmSig<-function(Mall, arms, mixes, thresh,log,log1 = FALSE, plot = FALSE,prop=1/3){
 respvs = array(dim = c(length(arms),length(mixes)))
 pvsall = list()
 for(i in 1:length(arms)){
   Mi = Mall[[which(names(Mall)==arms[i])]]
   pvsalli = list()

   for(j in 1:length(mixes)){
        Mij = Mi[[which(names(Mi)==mixes[j])]]
	pvsij = calcSigAll1(Mij, thresh,log,prop=prop)
        pvsalli[[j]] = pvsij
        title = paste(arms[i], mixes[j])    
        pv = cumulSig(pvsij,title, plot=plot, log =log,log1 = log1)
    	respvs[i,j] = pv
        if(plot & (j<length(mixes) | i<length(arms))) dev.new()
        print(paste(title,"pv=",sprintf("%5.2g",pv)))
   }
   names(pvsalli) = mixes
   pvsall[[i]] = pvsalli
 }
 names(pvsall) = arms
 list(pvsall=pvsall,respvs = respvs)
}

plotMatrix<-function(EI,loc,cols = 1:(dim(EI$Ep)[2])){
mat = EI$Ep
loc = EI$loc
x = loc[,2]
ymin = min(mat)
ymax = max(mat)
xmax = x[length(x)]
xmax = xmax + (xmax - x[1])*0.5
plot(x,mat[,1],col='white',type="p",ylim = c(ymin,ymax),xlim = c(x[1],xmax),pch=i,cex=0.1)
colours = cols
for(k in 1:length(cols)){
  i = cols[k]
 colours[k] = i
 if(i==1) colours[k] =  'pink' 
 else if(i==3) colours[k]=  'grey' 
 lines(x,mat[,i],col=colours[k],type="p",pch = i,cex = 0.1)
}
print(cols)
print(colours)
nme = dimnames(mat)[[2]][cols]
legend("bottomright",nme,col = colours,pch = cols)
}
