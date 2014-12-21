##pv heatmap fuction
.heatm<-function(pvs, logfun, opts =list(title="",indexMatch = NULL, Colv=TRUE),symbreaks = FALSE, range = c(0, 20), file ="", xlab="phenotypes",ylab="snp", key=T){
  if(!is.list(pvs)) pvs = list(pvs)
  #par(mfrow = c(1,length(pvs)))
  title=opts$title
  if(symbreaks) col = c(rev(colorRampPalette(brewer.pal(9,"Reds"),bias = 0.5)(250)), colorRampPalette(brewer.pal(9,"Blues"), bias = 0.5)(250))
  else col=colorRampPalette(brewer.pal(9,"Blues"))(250)
  Colv=opts$mPhen.Colv
  
  if(symbreaks){
   min1 =-1
   max1 = 1
   min11 = -10
   max11 = 10
  }else{
    min11 = range[1]
    min1 =0
    max1 =1
    max11 = range[2]
  }
  for(k in 1:length(pvs)){
    pvm = pvs[[k]]
    if(is.null(pvm)){
      return(NULL)
    }
    if(!is.matrix(pvm)) pvm = as.matrix(pvm)
    heatMatrix66 = logfun(t(pvm))
     min1 = max(min11,min(min1, min(heatMatrix66,na.rm=T),na.rm=T))
     max1 = min(max11,max(max1, max(heatMatrix66,na.rm=T),na.rm=T))
  }
  min1 = min1 - 0.2*abs(min1)
  max1 = max1 + 0.2*abs(max1)
  min1 = floor(10*min1)/10
  max1 = ceiling(10*max1)/10
  if(symbreaks){
     max1 = max(-min1,max1,na.rm = T)
     min1 = -max1
  }  

  if(is.na(min1) || is.na(max1) || min1 == -Inf || max1 ==Inf || min1==Inf || max1==-Inf) breaks = NULL
  else breaks = sort(seq(min1,max1,length.out = (1+length(col))))
#  print(breaks)
#   breaks = NULL
  indexMatch = opts$mPhen.indexMatch
  dendro = Colv
  for(k in 1:length(pvs)){
    main = paste(title,names(pvs)[[k]], if(symbreaks) "beta" else "pv")
    pvm = pvs[[k]]
    if(is.null(pvm)){
      return(NULL)
    }
    if(!is.matrix(pvm)) pvm = as.matrix(pvm)
    heatMatrix66 = logfun(t(pvm))
    d = dim(heatMatrix66)
    if(d[1]==1) heatMatrix66 = rbind(heatMatrix66,rep(0,d[2]))
    d = dim(heatMatrix66)
    if(d[2] ==1) heatMatrix66 = cbind(heatMatrix66,rep(0,d[1]))  #1e-10*runif(d[1]))
    d = dim(heatMatrix66)
    dn = dimnames(heatMatrix66)
    cols66<-rep("white",ncol(heatMatrix66))
    if(!is.null(indexMatch)){
      indexChange<-grep(indexMatch,colnames(heatMatrix66))
      cols66[indexChange]<-"red"
    }
    cexCol = 0.5*min(1,0.2 + 1/log10(d[1]))
    cexRow =  0.5*min(1,0.2 + 1/log10(d[2]))
  
  
#    if(abs(max(heatMatrix66,na.rm=T) - min(heatMatrix66,na.rm=T))<1e-5) return(NULL)
    keysize = 0
    if(key == F) keysize = 0

par(oma=c(3,8,2,4))
     heatmap.2(t(heatMatrix66),Rowv=NA, main = title, xlab = "loci",key=F, keysize =0,
                                   na.rm=TRUE,breaks = breaks,
                                   Colv=NA, col = col,dendrogram = "none",
                                   density.info="none",trace="none", margins = c(3,3),
                                   symbreaks = symbreaks)

          p3 <- 
    tryCatch(
             
          heatmap.2(t(heatMatrix66),Rowv=NA, main = title, xlab = "loci", 
                    na.rm=TRUE,breaks = breaks,
                    Colv=NA, col = col,dendrogram = "none",
                    density.info="none",trace="none",
                    symbreaks = symbreaks)
      ,error = function(e) NULL)
    
    if(is.null(p3)){
p3 <- tryCatch(heatmap.2(t(heatMatrix66),Rowv=NA,
                    na.rm=TRUE,breaks =breaks,
                    Colv=NA, col = col,dendrogram = "none",
                    density.info="none",trace="none",xlab=xlab,ylab=ylab,key=FALSE,margins=c(8,8),
                    RowSideColors = c(cols66),scale="none",cexRow = cexRow,cexCol=cexCol,main = main,
                    symbreaks = symbreaks), error = function(e) NULL)
     }
     if(is.null(p3)){
warning(paste('unable to make heatmap',file[1], main[1]))
##print(heatMatrix66)
     }else{
       if(symbreaks) legend(x="topleft",legend="beta")
       else legend(x="topleft",legend="-log(p)")
       if(k==1)dendro = p3$colDendrogram
    }
  }
   invisible(dendro)
}
