
readVCF  <-  function (prefix, chr)
{

    addup   <- function(v) {c(v[1] + v[2], v[3] + v[4])}
    chr = paste("chr", chr, sep = "");
    aFile=paste(prefix, chr, "vcf", sep=".")
    #xx = read.table(aFile, comment.char = "#", nrows=140, stringsAsFactors=FALSE)
    xx = read.table(aFile, comment.char = "#", stringsAsFactors=FALSE)
    # remove INDELS and records without DP4 field
    sel = grep("INDEL", xx[[8]], invert=T)
    sel = intersect(grep("DP4", xx[[8]]), sel)
    xx = xx[sel,]
    # extract DP4 field from column 8
    m=regexec("DP4=[[:digit:]]+,[[:digit:]]+,[[:digit:]]+,[[:digit:]]+;", xx[[8]])
    tmp=regmatches( xx[[8]], m)
    #remove pos with no DP4 field

    digits=lapply(tmp, strsplit, split="[,=;]")
    nrow = dim(xx)[1]
    #res = as.data.frame(matrix(nrow=nrow, ncol=5))  # pos   ref alt   freq1 freq2
    res = (matrix(nrow=nrow, ncol=5))  # pos   ref alt   freq1 freq2
    for (i in 1:nrow)
    {
        res[i,4:5] =  addup(as.integer(digits[[i]][[1]][2:(2+3)]))
    }
    res = as.data.frame(res, stringsAsFactors = F)
    res[,1] = ( xx[[2]])
    res[,2] = ( xx[[4]])
    res[,3] = ( xx[[5]])
    rownames(res) = unlist( xx[[2]])
    colnames(res) = c("pos","ref", "alt", "NRef", "NAlt")
    res
}




# -------------------------------------------------------------------
#  tt specifies type of test: 2 for chiSq, 1 for betabinom
getPValues <- function(tB, nB, tSum, nSum, thresh, c = 0, RCN = -1, tt=1)
{
    R  <- cbind(tB[tSum > 1], nB[tSum > 1])
    S  <- cbind(tSum[tSum > 1], nSum[tSum > 1])
    res            = rep(NA,length(tSum))

    if(length(RCN) == 1 || RCN == -1)
    {
        Pi <- cbind(R[,1] / S[,1], R[,2] / S[,2])  # when tSum =0,  the /0 problem arises.
    }
    else
    {
        #                              (1 - c) * m_i_n + c * m_i_n * x
        #     P = f(c, RCN) = ──────────────────────────────────────────────────────
        #                      (1-c) * d_i_n + c * m_i_n * x + c *(d_i_n- m_i_n) * y
        x = RCN[1]
        y = RCN[2]

        m_i_n  = R[,2]
        d_i_n  = S[,2]
        D = (1-c) * d_i_n + c * m_i_n * x + c *(d_i_n- m_i_n) * y
        P  = ( (1-c) * m_i_n + c * m_i_n * x ) / as.double(D)
        Pi <- cbind(R[,1] / S[,1], P)  # when tSum =0,  the /0 problem arises.
    }
    log       <- F
    # pvs <- calcSigAll(Pi, R, S, thresh = -1, log = log, prop = 1, tt = tt)
    pvs <- calcSigAll(Pi, R, S, thresh = thresh, log = log, prop = 1, tt = tt, c = c, RCN = RCN)
    res[tSum > 1]  = pvs$all
    res[tSum <= 1] = 0
    list(all=res)
}

# -------------------------------------------------------------------------
# range for the min, max PValue to be maped to the color range of heatmap.
plotHeatm <- function (pvs, label, range=c(0, 20))
{
    xx = matrix (0, nrow =2, ncol = length(pvs$all))
    colnames(xx) = label
    xx [1,] = pvs$all
    xx [2,] = pvs$all
    logfun<-function(x) -log10(x)
    .heatm(xx[,1:200], logfun, range = range)
}


# -------------------------------------------------------------------
# This plots a single strip of vertical bars in different colors.
plotBars <- function (pvs, coord, level, width=1, thick = 1,colScale=250, ct = 1, ceiling = 1)
{
    if (ct==1)
    {
        col=colorRampPalette(brewer.pal(9,"Blues"))(colScale)
    } else {
        col=colorRampPalette(brewer.pal(9,"YlOrRd"))(colScale)
    }
    for (num in 1:length(pvs))
    {
        aPV = pvs[num]
        if(!is.na(aPV) && aPV < ceiling)
        {
        cidx = as.integer(colScale * aPV)
        lines(rep(coord[num],2),seq(level,width+level),col=col[cidx], lwd=thick)
        }
    }
}


average <- function (coord, rSquare)
{
    res =rep(0,length(coord))
    theNames =rep(0,length(coord))
    sum = 0
    pre = -1
    idx = 1
    count = 1
    for (num in 1:length(res))
    {
        if(pre != coord[num])
        {
            res[idx] = sum /count
            theNames[idx] = pre

            idx = idx +1
            count = 1
            pre = coord[num]
            sum =  rSquare[num] 

        }
        else
        {
            count = 1 + count

            sum = sum + rSquare[num]
        }
    }
    res[idx] = sum / count

    theNames[idx] = pre
    names(res) <- theNames
    print(names(res)[1:10])
res[2:idx]
}



normalize <- function(data)
{
        min = min(data)
    max = max(data)
        as.double(data - min) / (max-min) 
}




# ---------------------------------------------------------------------------------------------------
# Draw a heatmap in which color is determined by the scale of PValue (pvSet).
#   pvSet is a matrix m sets of pvs.
#   Each elements of pvs[m] is the pvalue at SNPi or Regioni which has a coordinate.
#
#   tt determines the way to draw the heatmap:
#    1.  Draw with the bars along with actual coodinates of the SNPi.     (default)
#    2.  Draw with the bars along with the order of the occurance of SNPi. 
#
plotHeatmap  <- function (coord, pvSet, tt=1, font=5, barWidth=1, barThick=1, cScale=250, xlab="", ylab="", main="")
{
    if (tt == 1)
    {
    endAt   = tail(coord, 1) + 40
    startAt = head(coord, 1) - 40
    # For just opening the canvas, otherwise line function won't work. The plot won't show in the graph.
    plot(endAt,dim(pvSet)[1] +2, ylim=c(1,dim(pvSet)[1] +1), xlim=c(startAt ,endAt),  xlab=xlab, ylab=ylab, yaxt='n', ann=FALSE,  main=theTitle)

    for (ii in 1:dim(pvSet)[1])
    {
        text(-40 , dim(pvSet)[1] +1 - ii +0.5, labels = rownames(pvSet)[ii], font=5, offset=1.5, adj=1)
        plotBars(pvSet[ii,], coord,  dim(pvSet)[1] +1 - ii)
    }
    } else {
        opts <- list(title = main, indexMatch = NULL, Colv=TRUE)
        #.heatm(pvSet,  logfun = function(x){ x}, range = c(0,20), opts = opts, key=F )
        .heatm(pvSet,  logfun = function(x){ x}, range = c(0,20), opts = opts )
    }
}

# --------------------------------------------------------------------------------------------------
#   Read in the frequency data includes the A,B allele frequences in tumor and control.
#   Parameters :
#     dir     is the directory where the frequency data are.
#     merge   is to specify how to merge the frequences at SNPi into regions. Don't merge by default.
#     rm.zero is for remove zero frequency records. Remove zero by default.
#
loadFreq  <- function(dir, chrID, merge, doflip=T, rm.zero = T, ana = "p.tumor", nPrefix, tPrefix)
{
    #  how to merge regions
    #preprocess = merge
    forceRead  = T
    countDataF = paste(ana, paste(chrID, "dat", sep="."), sep=".")
    print(countDataF)
    if ( !forceRead & file.exists(countDataF)) 
    {
        print(sprintf("loading chr%s from object", chrID) )
        load(countDataF)
    } else
    {
        print(sprintf("loading chr%s from file", chrID) )
        tAFreq = readVCF(tPrefix, chrID)
        nAFreq = readVCF(nPrefix, chrID)

        aFile=paste(nPrefix, paste("chr", chrID, sep=""), "haps", sep=".")
        rules  <- read.table(aFile,  stringsAsFactors = F)
        rules=rules[((rules[[6]] + rules[[7]]) ==1),]
        rules=(rules[, c(3, 6)])
        colnames(rules) = c("pos", "flip1")
        save(tAFreq, nAFreq, rules, file=countDataF)
    }
    
    # align the rows
    rules [with(rules,  order(pos)), ]
    tAFreq[with(tAFreq, order(pos)), ]
    nAFreq[with(nAFreq, order(pos)), ]

    multiAlt = grepl("[^.]*,", nAFreq$alt)

    # the loci in normal with multiple alternative alleles.
    multiAlts = c("loci"            =length(nAFreq$NAlt), 
                  "multiAlt"       =length(nAFreq$NAlt[multiAlt]),
                  "without ref"     = length(nAFreq$NRef[multiAlt & (nAFreq$NRef==0)])) 
    #cat("\t# overall SNP loci:\tmulti alt:\twithout ref: \n")
    #cat(sprintf("\t\t%d\t\t%d\t\t%d\n", length(nAFreq$NAlt),  length(nAFreq$NAlt[multiAlt]), length(nAFreq$NRef[multiAlt & (nAFreq$NRef==0)])))


    # align ref haps with normal snps
    rules  = rules [rules$pos  %in% nAFreq$pos, ]
    nAFreq = nAFreq[nAFreq$pos %in% rules$pos, ]

    if(nrow(nAFreq) == 0)
    {
        print("[error] Could not phase SNPs. Check if *.haps file is empty.")
    }

      # align tumor snps. This gives 
    #tSnps = data.frame( NA, dim(nAFreq), rownames=rownames(nAFreq), colnames(nAFreq))
    tSnps <- nAFreq 
    tSnps$NRef = NA
    tSnps$NAlt = NA
    tSnps$pos  = nAFreq$pos
    tSnps[tSnps$pos %in% tAFreq$pos, ] = tAFreq[tAFreq$pos %in% tSnps$pos, ]

    tAFreq = tSnps
    missingIdx = is.na(tSnps$NAlt)
    naInTumor = c("loci"      = length(tAFreq$NAlt), 
                  "missings"  = length(tAFreq$NAlt[is.na(tSnps$NAlt)])) 
    tAFreq$NAlt[missingIdx] = 1
    tAFreq$NRef[missingIdx] = 1

        # generated phased allelic depth
    if(doflip == T)
    {
        flipSites = rules$flip1 == 1
        nalt  = nAFreq$NAlt [ flipSites]
        nref  = nAFreq$NRef [ flipSites]
        nAFreq$NAlt [flipSites ] = nref  
        nAFreq$NRef [flipSites ] = nalt  

        nalt  = tAFreq$NAlt [flipSites]
        nref  = tAFreq$NRef [flipSites]
        tAFreq$NAlt [flipSites] = nref  
        tAFreq$NRef [flipSites] = nalt  
    }



    

    # heterogenous alts (normal.alt != tumor.alt)
    #  1.   alt_normal   alt_tumor    
    #           T              A
    #     This is inconsist
    #
    #  2.   alt_normal   alt_tumor
    #           T              .
    #     Consist, but the num of alt in tumor has to be zero

    heteroAlt=c( "loci"        = length(tAFreq$NAlt),
                 "heteroAlts"  =  0)
    #inconsis = ((tAFreq$alt != nAFreq$alt) & (tAFreq$NAlt > 4))

    inconsis = ((tAFreq$alt != nAFreq$alt) & (tAFreq$alt != "."))
    #inconsis = (tAFreq$alt != nAFreq$alt)
    #tAFreq$NAlt[tAFreq$alt == "." & tAFreq$NAlt < 5] = 0
    if(length(which(inconsis, T)))
    {
        numOfInconsis = length(which(inconsis, T))
        tAFreq = tAFreq[!inconsis,]
        nAFreq = nAFreq[!inconsis,]
        heteroAlt = c("loci"        = length(tAFreq$NAlt), 
                      "heteroAlts"  = numOfInconsis) 

    }

    if ( nrow(tAFreq) != nrow(nAFreq) )  stop(sprintf("Frequnencies from two samples are not well aligned : %d!=%d", nrow(tAFreq), nrow(nAFreq)) )


#    if (preprocess == "mergeSite") 
#    {
#        nAFreq <- merge(nAFreq[,c("NRef","NAlt")], len = 20)
#        tAFreq <- merge(tAFreq[,c("NRef","NAlt")], len = 20)
#    } else if (preprocess == "mergeDist") 
#    {
#        nAFreq <- mergeByDist(nAFreq[,c("NRef","NAlt")], dist=10000)
#        tAFreq <- mergeByDist(tAFreq[,c("NRef","NAlt")], dist=10000)
#    }

  
    nB     <- nAFreq$NRef
    nSum   <- nAFreq$NRef + nAFreq$NAlt
    tB     <- tAFreq$NRef
    tSum   <- tAFreq$NRef + tAFreq$NAlt

    minCount = tSum > 2
    nB    = nB   [minCount]
    nSum  = nSum [minCount]
    tB    = tB   [minCount]
    tSum  = tSum [minCount]

    names(tSum) = rownames(tAFreq)[minCount]
    names(nSum) = rownames(nAFreq)[minCount]
    names(tB)   = rownames(tAFreq)[minCount]
    names(nB)   = rownames(nAFreq)[minCount]



#    if (rm.zero == T)
#    {
#        nonZero = !((tSum - tB) == 0 | tB  == 0 | (cSum - cB) == 0  | cB == 0 )
#        tSum    = tSum[nonZero]
#        nSum    = nSum[nonZero]
#        tB      = tB  [nonZero]
#        nB      = nB  [nonZero]
#    }
#
    list(tSum=tSum, nSum=nSum, tB=tB, nB=nB, multiAlts = multiAlts, naInTumor = naInTumor, heteroAlt = heteroAlt) 

}

maxll <- function(tab,genotypes, tc)
{
    mll = 0
    state = 0
    res  <- .C( "maxll_improv",
       as.integer(as.vector(t(tab))),
       as.integer(dim(tab)[1]),
       as.integer(genotypes),
       as.integer(length(genotypes)/2),
       as.double(tc),
       state = as.integer(state),
       mll   = as.numeric(mll),
       PACKAGE = "cnProfile"
       )
       
  #  for(k in 1:(length(genotypes)/2))
  #  {
  #      sum = 0
  #      x=genotypes[k*2-1]
  #      y=genotypes[k*2]
  #      for(n in 1:dim(tab)[1])
  #      {
  #         each = tab[n,]
  #         prob= each$nB * x / ((each$nSum -each$nB) * y + each$nB * x)
  #         if(prob ==0)
  #         print(prob)
  #         #sum = sum + each$tB * log(prob) + (each$tSum - each$tB) * log(1-prob)
  #         sum = sum + log(pbinom(each$tB, each$tSum, prob))  
  #      }

  #      if(sum > mll)
  #      {
  #          mll = sum
  #          state = k
  #      }
  #  }
    res$mll
}

mergeByDistML<-function(tab,genotypes, tc=1, dist=10000)
{
    fluctuation = c()  
    fullLabel   = c()  
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
    depths = matrix(0, nrow = idx -2, ncol = dim(tab)[2])
    rownames(depths) <-  rep("",idx -2)
    numDist = c()
    inc=0
    for(k in 2:(idx-1)){
        start = sep[k -1]
        end = sep[k]

        if(dim(tab[start:end,,drop=F])[1] > 3)
        {
        inc = inc +1
        res1 = apply(tab[start:end,,drop=F],2,sum)
        fluctuation[inc] <- maxll(tab[start:end,,drop=F], tc=tc, genotypes) # lexically scoped
        depths[inc,] = res1
        rownames(depths)[inc] = toString(pos[start])
        fullLabel[(inc *2 -1):(inc * 2)]   <- rownames(tab) [c(start,end)] 
        }
    }

    names(fluctuation) <- rownames(depths[1:inc,])
    colnames(depths) = colnames(tab)
    list("fluctuation" = fluctuation, 
         "fullLabel"   = fullLabel,
         "depths"       = data.frame(depths[1:inc,]))

}

# calculating the Maximum log-likelihood for different genotypes. Then use ML to determine if there is band switching.
mergeUsingML <- function(tab, genotypes, tc=1, len=20)
{
    fluctuation = c()  
    fullLabel   = c()  
    rep = ceiling(dim(tab)[1]/len)
    res = matrix(0,nrow = rep,ncol = dim(tab)[2])
    for(k in 1:rep){
        start = (k-1)*len +1
        end = min(start +len,dim(tab)[1])
        res1 = apply(tab[start:end,,drop=F],2,sum)
        fluctuation[k] <- maxll(tab[start:end,,drop=F], genotypes, tc) 
        fullLabel[(k *2 -1):(k * 2)]   <- rownames(tab) [c(start,end)] 
        res[k,] = res1
    }

    dimnames(res)[[2]] = dimnames(tab)[[2]]
    colnames(res) = colnames(tab)
    rownames(res) = rownames(tab)[seq(1,length(rownames(tab)), len)]

    names(fluctuation) <- rownames(tab)[seq(1,length(rownames(tab)), len)]
    depths = data.frame(res)
    list("fluctuation" = fluctuation, 
         "fullLabel"   = fullLabel,
         "depths"       = depths)
}
filtering <- function (data, thresh)
{
    fluctuation = data[["fluctuation"]]
    fullLabel   = data[["fullLabel"]]
    depths      = data[["depths"]]
indexing<-function(idx) {as.numeric(names(idx))}

    selection =  fluctuation > thresh 
    #selection =  fluctuation > thresh & indexing(fluctuation) > 88000000 & indexing(fluctuation) < 89000000
    dim(fullLabel) = c(2,length(fullLabel)/2)
    fullLabel = t(data.frame((fullLabel),row.names=c("start", "end")))

    baf=depths$tB/depths$tSum
    hlRegion =  baf < 0.6 & baf > 0.4
    col = rep(1, length(fluctuation))
    col[hlRegion] = 2
    plot(indexing(fluctuation), fluctuation,pch=4,col=col, main = "Log likelihood at regions", 
         xlab="position (bp)", ylim=c(-300,0))
    #dev.new()
    plot(indexing(fluctuation), depths$tB/depths$tSum, pch=4, col=col, ylim=c(0,1), main="BAF before filtering")
    print("before filtering")

    col = rep(2, length(fluctuation))
    col[selection] = 1

    depths      = depths     [selection, ]
    fullLabel   = fullLabel  [selection, ]
    fluctuation = fluctuation[selection ]
    print(length(fluctuation))

    plot(indexing(fluctuation), depths$tB/depths$tSum, pch=4, col=col[selection], ylim=c(0,1),
         main="BAF after filtering", xlab="postion (bp)")
    plot(indexing(fluctuation), depths$nB/depths$nSum, pch=4, col=col[selection], main = "normal BAF", xlab="postion (bp)", ylim=c(0,1))
    #plot(indexing(fluctuation), depths$tSum, pch=4, col=col[selection], main = "normal RD", xlab="postion (bp)")
    #plot(indexing(fluctuation), depths$nSum, pch=4, col=col[selection], main = "tumor  RD", xlab="postion (bp)")
    list("fluctuation" = fluctuation, 
         "fullLabel"   = fullLabel,
         "depths"      = depths )
}

inferStates  <-  function (anaName, preprocess, chrID, runHMM = T, method = "Baum_Welch", 
                           nB,nSum, tB, tSum,
                           genotypes, maxiter=30)
{

    indexing<-function(idx) {as.numeric(names(idx))}

    dataFile = paste("res", anaName, preprocess,"chr", chrID,"dat", sep=".")
    if ( !runHMM & file.exists(dataFile)) 
    {
        load(dataFile)
    } else
    {
        method   =  "Baum_Welch"      #  "viterbi_train" or  "Baum_Welch" by default
        #DOA_range = c( 1.5, 2, 2.5)
        #DOA_range = c(0.9999999, 1.50000001)
        #DOA_range = c(1.4000001, 2.0000001, 2.50000001)
        DOA_range = c(0.800001, 1.4000001, 2.0000001, 2.50000001)

        allRes = list()

        for(start in c(1:(length(DOA_range)-1)))
        {
            print("notice")
            print(c(DOA_range[start], DOA_range[start +1]))
            if(!exists("maxRes"))
            {
                maxRes = segHMM2(rbind(nB,nSum,tB,tSum), 
                         length(genotypes)/2,
                         genotypes, 
                         maxiter = 100, 
                         DOA_range = c(DOA_range[start], DOA_range[start +1]),
                         method  = method)
                allRes[[start]] = maxRes
            } else{
                res= segHMM2(rbind(nB,nSum,tB,tSum), 
                         length(genotypes)/2,
                         genotypes, 
                         maxiter = 100, 
                         DOA_range = c(DOA_range[start], DOA_range[start +1]),
                         method  = method)
                allRes[[start]] = res

                print(maxRes$log.lik)
            # Do the inference for different DOA_ranges.
            # If the change is bigger than 1000, accept the better range.
                if(res$log.lik - maxRes$log.lik  > 1000) maxRes = res

                print(maxRes$log.lik)

            }


        }

        res=list()
        idx=1
        #for (cc in seq(0.001,0.999, by=1/as.double(numInits)))
      #  cc=0.95
      #  {
      #      tmp  <- segHMM2(rbind(nB,nSum,tB,tSum), length(RCNs$RCN), genotypes, tc = cc, maxiter=30, method=method)
      #      res[[idx]] = c(tmp$tc, tmp$log.lik,cc)

      #      if(tmp$log.lik > maxRes$log.lik)
      #      {
      #          maxRes = tmp
      #      }
      #      idx = idx +1
      #      print("============================================= \n");
      #  }
        #tmp=lapply(res, func <- function(aRes) {print(c(aRes$tc, aRes$log.lik))})
        tmp=lapply(res, func <- function(aRes) {print(aRes)})

        label = indexing(nB)
        save(nB,nSum, tB, tSum, label, maxRes, allRes, file=dataFile)
    }
    maxRes
}
plotPred  <- function(fluctuation, nB,nSum, tB, tSum, genotypes, maxRes, method = "Baum_Welch", len, chrID)
{

    indexing<-function(idx) {as.numeric(names(idx))}

    label    = indexing(fluctuation)

    roundTo <- function(x, k) format(round(x, k), nsmall=k)
    title = sprintf("copy_number.%s.merge%s.chr%s", method, len, chrID)
    showBound <- function(boundaries) 
    {

        for (ii in 1:length(boundaries))
        {
            b = boundaries[ii]
            lines(c(b, b), c(1, 7), col=3,  lty=4, lwd=3)
            text(b , 1 , labels = ii, font=5,  pos = 3, family="sans")
        }

    }


    ################################################################################################
    #    Slide 1:     This shows relative  read depth given the degree of amplification and tumor
    #                 cellularity.
    ################################################################################################


    sum_di_t = sum(tSum)
    sum_di_n = sum(nSum)
    rdTitle  = paste("Relative read depth", "tc", roundTo(maxRes$tc, 3), "DOA", roundTo(maxRes$DOA,3), sep="--")
    plot(label, tSum, main=rdTitle, xlab="loci", ylab="")
    plot(label, nSum, main=rdTitle, xlab="loci", ylab="")
    plot(label, tSum/nSum * (sum_di_n * maxRes$DOA/sum_di_t), main=rdTitle, xlab="loci", ylab="", ylim = c(0,4))

    print("rcov : ")
    print((sum_di_n * maxRes$DOA * maxRes$tc + (1 - maxRes$tc) * sum_di_n)/sum_di_t)

    ratio = sum_di_n * maxRes$DOA * maxRes$tc / (sum_di_n * maxRes$DOA * maxRes$tc + (1 - maxRes$tc) * sum_di_n)
    print(paste("ratio: ", ratio))
    print(paste("reads: ", ratio * sum_di_t))
    print(ratio*sum_di_t / (sum_di_n * maxRes$DOA * maxRes$tc))

    
    lines(c(min(label),max(label)),c(1/2, 1/2), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1,     1), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1.5, 1.5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(2,     2), col=3,  lty=4)
    lines(c(min(label),max(label)),c(2.5, 2.5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(3,     3), col=3,  lty=4)
    lines(c(min(label),max(label)),c(3.5, 3.5), col=3,  lty=4)
 

    ################################################################################################
    #    Slide 2:     This shows the maternal allelic freqency plot colored by copy
    #          numbers.
    ################################################################################################

    copyNum = genotypes[maxRes$hidden.states * 2] + 
       genotypes[maxRes$hidden.states * 2 -1]
    plot(label, tB/tSum , col=copyNum, pch=4, xlab= "loci", ylab="Maternal Allelic Frequency", main=title)
    lines(c(min(label),max(label)),c(1/2, 1/2), col=3,  lty=4)
    lines(c(min(label),max(label)),c(2/3, 2/3), col=3,  lty=4)
    lines(c(min(label),max(label)),c(3/4, 3/4), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/3, 1/3), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/4, 1/4), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/5, 1/5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(4/5, 4/5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/6, 1/6), col=3,  lty=4)
    lines(c(min(label),max(label)),c(5/6, 5/6), col=3,  lty=4)
    ################################################################################################
    #    Slide 3:      This shows the maternal allelic freqency plot colored by copy
    #          numbers.
    ################################################################################################


    plot(label, genotypes[maxRes$hidden.states * 2 -1] / copyNum, col = copyNum,
         pch=4, xlab = "pos", ylab = "genotype states",  main=title,  xaxt = "n")

    lines(c(min(label),max(label)),c(1/2, 1/2), col=3,  lty=4)
    lines(c(min(label),max(label)),c(2/3, 2/3), col=3,  lty=4)
    lines(c(min(label),max(label)),c(3/4, 3/4), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/3, 1/3), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/4, 1/4), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/5, 1/5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(4/5, 4/5), col=3,  lty=4)
    lines(c(min(label),max(label)),c(1/6, 1/6), col=3,  lty=4)
    lines(c(min(label),max(label)),c(5/6, 5/6), col=3,  lty=4)

    legend("topleft", inset=.05, title="copy number", legend=c(1:7), fill=c(1:7))


    #if(exists("breakPoints"))   showBound(breakPoints)
    ################################################################################################
    #     Slide 4:     This shows the cn states in comparision with climat.
    ################################################################################################

    plot(label, copyNum,
         pch=4, xlab = "pos", ylab = "CN states",  main=title,  xaxt = "n")

    #if(exists("breakPoints"))   showBound(breakPoints)
    #axis(side = 1, at = seq(1, length(label), by = length(label)/15), 
    #     labels = label[seq(1, length(label), by = length(label)/15)])
    #maxRes = segHMM2(rbind(nB,nSum,tB,tSum), length(RCNs$RCN), genotypes, tc = 0.001, maxiter=30,fix_tc=T)
    breakPoints = c()
    preState = -1
    for(i in 1:length(maxRes$hidden.states))
    {

        state = copyNum[i]
        if(state != preState)
        {
            breakPoints[length(breakPoints) +1] =i
            preState = state
        }

        lines(fullLabel[i, ], c(state,state),  col=3,  lty=1,lwd=5)
    }

    length(breakPoints)
    print(breakPoints[1:2])
    print(label[breakPoints[1:2]])
    #axis(side = 1, at = label[breakPoints], labels = round((as.integer(label[breakPoints])/1e+6), digit=1))


    #************************************************************************
    #   loading climat
    #************************************************************************

 ###    climat=read.table("./climat/G15511_HCC1143_1_sorted_chr8_alleleDepth.results", skip=12, header=T)
 ###    climat$StartPos
 ###    climat$EndPos
 ###    climat$CN
 ###    for(i in 1:length(climat$StartPos))
 ###    {
 ###        lines(c(climat$StartPos[i], climat$EndPos[i]), 
 ###              c(climat$CN[i]+0.1,climat$CN[i]+0.1), 
 ###              col=2,  lty=1,lwd=5)
 ###    }



    genoLabel = c("A", "AB", "ABB", "ABBB", "ABBBB" )
    #genoLabel = c("M", "MP", "MPP", "MPPP", "MPPPP" )
    for(ii in 1:length(genoLabel))
    {

        text(-40 , ii , labels = genoLabel[ii], font=5, offset=10.5 , adj=1)
    }
    legend("topleft", inset=.05, title="Methods", c("Climat","mine"), fill=c(2,3))

    #axis(side = 1, at = seq(1, length(label), by = as.integer(length(label)/15)), labels = label[seq(1, length(label), by = as.integer(length(label)/15))])
}



loadBamRegions  <- function (bamF, chrID, start,len, prefix, winSize = 10000, qc = 2)
{
    forceRead = F
    datFile   = sprintf("%s.chr%s.readDepth.win.%d.dat", prefix, chrID, winSize)

    if ( !forceRead & file.exists(datFile)) 
    {
            print(sprintf("loading chr%s from object", chrID) )
            load(datFile)
    } else
    {

        rDepths = rep(0, length(start))
        rFailure = F
        all <- .C("countReads", 
                  as.character(bamF), 
                  as.character(chrID), 
                  as.integer(start),
                  as.integer(len), 
                  as.integer(length(start)),
                  as.integer(qc), 
                  rDepths = as.integer(rDepths),
                  rFailure = as.logical(rFailure))

        if (all$rFailure) return(NA)
        names(all$rDepths) = start
        save(all, file=datFile)
    }
    all$rDepths
}



#============================================================================
#  ---------          FUNCTION      plotMAF         -------------------
# This is for plot of Major allele frequence with or without the LD profile.
#    LDdata is a list of two vectors for coords and r-square.
#
plotMAF <- function(pos, baf, prange=1, LDdata="", interval=1,  title="test", ceiling=1, ttp = 2)
{

    beginAt = head(pos,1)
    endAt   = tail(pos,1)
    if (length(interval) ==2 )
    {
        beginAt = interval[1]
        endAt   = interval[2]
    }
    range=beginAt < pos & pos < endAt
    rbind(pos[range], baf[range])
    if (ttp == 1)
    {
        plot(pos[range],baf[range],pch=20,col=3, xaxs ="i", ylim=c(0,1), xaxp=c(beginAt,endAt,20), las =2, ylab="MAF", xlab="")
        #plot(pos[range],baf[range],pch=20,col=3, xaxs ="i", ylim=c(0,1))
        if( is.list(LDdata) )
        {
            rspos=LDdata[[1]]
            r=LDdata[[2]]
            range=beginAt < rspos & rspos < endAt
            plotBars(r[range], rspos[range], level=0, ct=2, ceiling=ceiling)
        }
    }

    if (ttp  == 2)
    {
        colScale = 250
        col=colorRampPalette(brewer.pal(9,"YlOrRd"))(colScale)

        if(  is.list(LDdata) )
        {
            rspos = LDdata[[1]]
            r     = LDdata[[2]]
            rspos = rspos[rspos%in%pos]
            r     = r[rspos%in%pos]
            pos   = pos[pos%in%rspos]
            baf   = baf[pos%in%rspos]
            range2=beginAt < rspos & rspos < endAt

            if(is.vector(prange))
            {
                print(prange)

                plot(pos[range], baf[range], pch=20, col=col[as.integer(r[range] * colScale)],xlim=prange, ylim=c(0,1), xaxp=c(beginAt,endAt,20), las =2, ylab="MAF", xlab="", main=title)
            }else
            {
                plot(pos[range], baf[range], pch=20, col=col[as.integer(r[range] * colScale)], ylim=c(0,1), xaxp=c(beginAt,endAt,20), las =2, ylab="MAF", xlab="", main=title)
            }
        } else {
            print (length(pos[range]))
            print (length(baf[range]))
            print (length(r[range]))
            plot(pos[range], baf[range], pch=20,  ylim=c(0,1), xaxp=c(beginAt,endAt,20), las =2, ylab="MAF", xlab="", main=title)
        }
    }
}



#
#  CN(A)  CN(B)  RCN    Genotype
#  1        3      3     ABBB
#  1        2      2     ABB
#  1        1      1     AB
#  2        1      1/2   AAB
#  3        1      1/3   AAAB
#
# in which RCN is defined as 
#             CN(B)
#   RCN = ──────────────
#             CN(A)
#   
#   
# The posterier probabbility of RCN = x is give as :
#   
#   
#                      P(D | RCN = x) * P(RCN = x)
#   P(RCN = x | D) = ───────────────────────────────
#                      sum [ P(D | RCN) * P(RCN) ]
#   
#
#
#   P(D | RCN = x) = P(Nfb, Nfa, Tfb, Tfa | RCN) 
#                  => Nfb ~ Bino(Tfb + Tfa, P = f(c, RCN) )
# in which c is the cellular of tumor, aka the percentage of tumor in the mix.
#                                Nfb              Nfb * RCN 
#     P = f(c, RCN) = (1 - c) ───────────  + c ────────────────────
#                              Nfa + Nfb        Nfa + Nfb * RCN 
# The prior distribution of RCN is treated as a uniform distribution.
#    In this case,   P(RCN = x) = 1/5           


#  {M} {MMMP} {MMP} {MP} {MPP} {MPPP} {P}
# res  <- segHMM(baf, length(RCNs$RCN), c(0.999, 3/4, 2/3, 1/2, 1/3, 1/4, 0.001))

#  Scaling is to bring down the emission probility. 
#    This cause problems when it is being too big, 
#    because the transition probility needs to be summed on the probility (not logProb).
#    The sum will be 0, and an error occurs when divide 0.
#    Fix_tc determine if tc need to be estimated. If T, tc will be used as an initial guess.

segHMM2  <- function(allelicDepth, k, genotypes, tc = 0.5, fix_tc = F, diag.prob = 0.999999, maxiter = 10, eps = 0.0000001, DOA = 1.5, print.info = FALSE, method="Baum_Welch", DOA_range=c(1, 1.5) )
{
    numobs = dim(allelicDepth)[2]
    gamma     <- matrix((1 - diag.prob) / (k - 1), k, k)
    diag(gamma) <- diag.prob

    prob = rep(0.0, numobs * k)
        res <-
            .C(method,
               as.integer(numobs),
               as.integer(k),
               depths    = as.integer(as.vector(allelicDepth)),
               genotypes = as.integer(genotypes),
               tc        = as.double(tc),
               gamma     = as.double(gamma),
               pi        = as.double(rep(-log(k), k)),
               num.iter  = as.integer(maxiter),
               as.double(eps),
               log.lik = double(1),
               filtered.cond.probs = double(k * numobs),
               hidden.states = integer(numobs),
               prob          = as.double(prob),
               as.logical(fix_tc),
               DOA=as.double(DOA),
               DOA_range=as.double(DOA_range),
               as.logical(print.info),
               PACKAGE = "cnProfile"
               )
        res$hidden.states <- res$hidden.states + 1
        res$gamma <- matrix(res$gamma, nr = k)
        res
}







