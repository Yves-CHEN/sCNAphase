require(stats)


qb<-function(np,conf){
  Obs = rep(0,np)
  for(s in 1:np)
              {
                Obs[s]=qbeta(conf, s, np-s+1)
              }
 Obs
}

collapse<-function(pvs,title=expression(paste("-log"[10],"QQ"))){
dimp = dim(pvs)
nosnp = dimp[[1]]
nophe = dimp[[2]]
pv = matrix(0,nrow = nosnp*nophe,ncol=5)
for(i in 1:nophe){
   pv[((i-1)*nosnp+1):(i*nosnp) ,1] = pvs[,i]
   pv[((i-1)*nosnp+1):(i*nosnp) ,2] = rep(i,nosnp)
   
}
pv = pv[!is.na(pv[,1]),]
pv = pv[order(pv[,1]),]
np = dim(pv)[[1]]
pv[,3] = (1:np)/(np+1)
pv[,4] = qb(np,0.95)
pv[,5] = qb(np,0.05)
ymax =       ymax=-log10(min(pv))		 
  plot(-log10(pv[,c(3,1)]),  xlim=c(0,log10(np+1)), ylim=c(0,ymax), xlab=expression(paste("-log"[10]," (Expected p-value)")), ylab=expression(paste("-log"[10]," (Observed p-value)")), main=title, pch=20, cex=0.5,col=0)
 lines(-log10(pv[,3]), -log10(pv[,3]), col="grey", lw=2)
 lines(-log10(pv[,4]), -log10(pv[,3]), col="grey", lw=2)
 lines(-log10(pv[,5]), -log10(pv[,3]), col="grey", lw=2)
for(i in 1:nophe){
    lines(-log10(pv[pv[,2]==i,c(3,1)]),  xlim=c(0,log10(np+1)), ylim=c(0,ymax), xlab=expression(paste("-log"[10]," (Expected p-value)")), ylab=expression(paste("-log"[10]," (Observed p-value)")), main=title, pch=i, cex=0.5,col=i,type="p")
 
}
legend("bottomright", legend =dimnames(pvs)[[2]],cex=0.5,col=(1:nophe),lty=1,pch=1:nophe)

}

logQQ<-function(pvs, conf=NA, yEqx=TRUE, cols=NA, pointcols=NA,colours=1, title=expression(paste("-log"[10],"QQ")), ymax=NA)
  {
    pv=pvs[!is.na(pvs)]
    if(sum(pv==0)>0)
    {
      print("Warning: Ommitting pvalue=0")
      pv=pv[pv!=0]
    }
    #print (paste("remove ones : ", length(pv[pv ==1])))
    #pv = pv[pv!=1]
    pv=sort(pv)
    np=length(pv)

    if(is.na(pointcols))
      {
        pointcols=rep(1,np)
      }
    if(is.na(sum(cols)))
      {
        if(yEqx)
          cols=1:(length(conf)+1)+1
        else
          cols=1:length(conf)+2
      }

    expec=rep(NA, np)
    for(s in 1:np)
      expec[s]=s/(np+1)

    plotpv=rep(FALSE, np)
    plotpv[1:min(1000,length(pv))]=TRUE
    plotpv[floor(10^(runif(n=4000,min=0,max=log10(length(pv)))))]=TRUE

    expec=expec[plotpv]
    pv=pv[plotpv]
    if(is.na(ymax))
      ymax=-log10(min(pv))		

   # png(outf)
    plot(cbind(-log10(expec), -log10(pv)),  xlim=c(0,log10(np+1)), ylim=c(0,ymax), xlab=expression(paste("-log"[10]," (Expected p-value)")), ylab=expression(paste("-log"[10]," (Observed p-value)")), main=title, pch=20, cex=0.5,col=1)


    if(yEqx)
      abline(0,1, lw=2, col=cols[1])
    
    if(!is.na(sum(conf)))
      {
        Obs=rep(NA, np)
        for(c in 1:length(conf))
          {
            for(s in 1:np)
              {
                Obs[s]=qbeta(conf[c], s, np-s+1)
              }
            Obs=Obs[plotpv]
            #expec=expec[plotpv]
            lines(-log10(expec), -log10(Obs), col=cols[c+as.numeric(yEqx)], lw=2)
          }
      }
    #dev.off()
  }





plotMatrix<-function(mat,loc,cols = NULL,title="",cex=0.5,chrom=NULL){
if(is.null(cols)) cols = 1:(dim(mat)[2])
if(!is.null(chrom)){
 inds = which(loc[,1] %in% chrom)
 mat = mat[inds,,drop=F]
 loc = loc[inds,,drop=F]
 
}
print(dim(mat))
arms = levels(as.factor(loc[,1]))
arms1 = as.numeric(unique((as.numeric(sub('q','',sub('p','',arms))))))
chr = loc[,1]
chr = as.numeric(sub('q','',sub('p','',chr)))
x = as.numeric(loc[,2])
x = (x-min(x))/((max(x)-min(x))*1.2)+0.1
for(i in 1:length(arms1)){
 chrind = which(chr==arms1[i])
 x[chrind] = arms1[i] + x[chrind]
}
print(head(mat))
ord = order(x)

x = x[ord]
mat = mat[ord,,drop=F]

ymin = min(mat[,cols])
ymax = max(mat[,cols])
print(ymin)
print(ymax)
xmax = x[length(x)]
xmax = xmax + (xmax - x[1])*0.1
colours = cols
for(k in 1:length(cols)){
  i = cols[k]
 colours[k] = i
# if(i==1) colours[k] =  'pink' 
# else if(i==3) colours[k]=  'grey' 
# else if(i==8) colours[k]='orange'
 if(k==1) plot(x,mat[,i],col=colours[k],type="p",ylim = c(ymin,ymax),xlim = c(x[1],xmax),pch=i,cex=cex,main=title)
 lines(x,mat[,i],col=colours[k],type="p",pch = i,cex = cex)
}
print(cols)
print(colours)
nme = dimnames(mat)[[2]][cols]
legend("bottomright",nme,col = colours,pch = cols)
}


#plotMatrix(-log10(table[,4,drop=F]),table[2:3],cols = 1,chr=10)


#index =9
#infile="file.txt";
#table = read.table(infile, sep = "\t", header=T);
#pdf(file=paste(infile,"pdf",sep='.'));
#logQQ(table[,index])
#dev.off();
#chsq = qchisq(table[,index],1,lower.tail=F)
#lambda = max(median(chsq),0.455)/0.455
#print(paste("lambda is ",lambda))
#table[,index] = pchisq(chsq/lambda,1,lower.tail=F)
#print(paste("inflation_factor is ", lambda))
#names(table)[9] = paste("pvalue", lambda, sep="_")
#write.table(table, file = paste(infile,"qq.txt", sep="_"),quote=F,col.names=T, row.names=F, sep="\t");
