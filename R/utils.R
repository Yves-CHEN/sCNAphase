

readVCF  <-  function (prefix, chr)
{

    addup   <- function(v) {c(v[1] + v[2], v[3] + v[4])}
    chr = paste("chr", chr, sep = "");
    aFile=paste(prefix, chr, "vcf", sep=".")
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
#     merge   is to specify how to merge the frequences at SNPi into regions. Don't merge by default.
#     rm.zero is for remove zero frequency records. Remove zero by default.
#
loadFreq  <- function(chrID, merge, doflip=T, rm.zero = T, ana = "p.tumor", nPrefix, tPrefix,  forceRead  = F, allelicMapability = F)
{
    #  how to merge regions
    #preprocess = merge
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
    tAFreq$NAlt[missingIdx] = 0.5
    tAFreq$NRef[missingIdx] = 0.5

    # 1. correct the different number of read mapped to ref than alt.

    if(allelicMapability == T)
    {
        allelicMapability =  0.953
        nAFreq$NAlt  = (nAFreq$NAlt  / allelicMapability)
        tAFreq$NAlt  = (tAFreq$NAlt  / allelicMapability)
    }

    # 2. generated phased allelic depth
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
    

    # 3. heterogenous alts (normal.alt != tumor.alt)
    #  1.   alt_normal   alt_tumor    
    #           T              A
    #     This is inconsist
    #
    #  2.   alt_normal   alt_tumor
    #           T              G
    #     Consist, but the num of alt in tumor has to be zero

    heteroAlt=c( "loci"        = length(tAFreq$NAlt),
                 "heteroAlts"  =  0)
    #inconsis = ((tAFreq$alt != nAFreq$alt) & (tAFreq$NAlt > 4))

    inconsis = ((tAFreq$alt != nAFreq$alt) & (tAFreq$alt != "."))
    if(length(which(inconsis, T)))
    {
        numOfInconsis = length(which(inconsis, T))
     #   tAFreq = tAFreq[!inconsis,]
     #   nAFreq = nAFreq[!inconsis,]
        heteroAlt = c("loci"        = length(tAFreq$NAlt), 
                      "heteroAlts"  = numOfInconsis,
                      "inconsis"    = inconsis) 
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
    names(tSum) = rownames(tAFreq)
    names(nSum) = rownames(nAFreq)
    names(tB)   = rownames(tAFreq)
    names(nB)   = rownames(nAFreq)

    minCount = nSum > 5
    nB    = nB   [minCount]
    nSum  = nSum [minCount]
    tB    = tB   [minCount]
    tSum  = tSum [minCount]
    
     # 4. mix into 0.1% tumor into the sample, this way all tB > 0.
    mixit=F
    if (mixit == T)
    {
        tB = tB + 0.001 * nB
        tSum = tSum + 0.001 * nSum
    }
#
    list(tSum=tSum, nSum=nSum, tB=tB, nB=nB, naInTumor = naInTumor, heteroAlt = heteroAlt) 

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
       PACKAGE = "sCNAphase"
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

# calculating the Maximum log-likelihood for different genotypes. Then use ML to determine if there is a state switching.
mergeUsingML <- function(tab, genotypes, tc=0.99, len=20, maxSeg=1000000)
{
    sumUp  <- function(seg) { 

        #print(tab[seg[1]:seg[2],,drop=F])
        res1 = apply(tab[seg[1]:seg[2],,drop=F],2,sum)
        as.integer(res1)
    }
    consistency  <- function(seg, doCheck=T)
    {
        #print(str(tab[seg[1]:seg[2],,drop=F]))

        #print(str(transform(tab[seg[1]:seg[2],,drop=F], nB=as.integer(nB), nSum=as.integer(nSum),
                                                #tB=as.integer(tB), tSum=as.integer(tSum)) ))
        #print(dim(as.integer()))
        if(doCheck)

        {
            res1 = maxll(transform(tab[seg[1]:seg[2],,drop=F], nB=as.integer(nB),
                        nSum=as.integer(nSum),
                        tB=as.integer(tB),
                        tSum=as.integer(tSum)),
                        genotypes,tc)
        }else
        {
            res1 = 0
        }


        #res1 = maxll(as.integer(tab[seg[1]:seg[2],,drop=F]), genotypes, tc)
        res1
    }
    rep = ceiling(dim(tab)[1]/len)
    res = matrix(0,nrow = rep,ncol = dim(tab)[2])

    startAt = (c(1:rep)-1)*len +1
    endAt = startAt +len
    endAt[length(endAt)] = dim(tab)[1]
    # filter sparse snp region
    sel = as.numeric(rownames(tab)[endAt]) - as.numeric(rownames(tab)[startAt]) < maxSeg
    #tab = tab[sel,,drop=F]
    startAt = startAt[sel]
    endAt  =  endAt[sel]

    res=apply(cbind(startAt,endAt),1, sumUp)
    res=t(res)
    fluctuation=apply(cbind(startAt,endAt),1, consistency, doCheck = T)
    fullLabel = cbind(rownames(tab)[startAt], rownames(tab)[endAt])
    colnames(fullLabel) = c("start", "end")

    colnames(res) = colnames(tab)
    rownames(res) = rownames(tab)[startAt]

    names(fluctuation) <- rownames(tab)[startAt]
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
    #dim(fullLabel) = c(2,length(fullLabel)/2)
    #fullLabel = t(data.frame((fullLabel),row.names=c("start", "end")))

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
                           genotypes, maxiter=5, fullLabel = fullLabel, DOA_range = c(1.0, 1.4),
                           underate = 1, tpmType = 1, reRun = T, mlen = 40)
{

    dataFile = paste("res", anaName, preprocess,"chr", chrID,"dat", sep=".")

    if ( !runHMM & file.exists(dataFile)) 
    {
        load(dataFile)
    } else
    {
        allRes = list()
        for(start in c(1:(length(DOA_range)-1)))
        {
            print(c(DOA_range[start], DOA_range[start +1]))
            if(!exists("maxRes"))
            {
                maxRes = segHMM2(rbind(nB,nSum,tB,tSum), 
                         length(genotypes)/2,
                         genotypes, 
                         maxiter = maxiter, 
                         DOA_range = c(DOA_range[start], DOA_range[start +1]),
                         method  = method, underate=underate, tpmType = tpmType)
                allRes[[start]] = maxRes
            } else{
                res= segHMM2(rbind(nB,nSum,tB,tSum), 
                         length(genotypes)/2,
                         genotypes, 
                         maxiter = maxiter, 
                         DOA_range = c(DOA_range[start], DOA_range[start +1]),
                         method  = method, underate=underate, tpmType = tpmType)
                allRes[[start]] = res
            # Do the inference for different DOA_ranges.
            # If the change is bigger than 1000, accept the better range.
                if(res$log.lik - maxRes$log.lik  > 1) maxRes = res
            }
        }

    }
    allRes
}

runForPloidy <-  function (anaName, preprocess, chrID, runHMM = T, method = "Baum_Welch", 
                           nB,nSum, tB, tSum, maxRes,
                           genotypes, maxiter=5, fullLabel = fullLabel, DOA_range = c(1.0, 1.4),
                           underate = 1, tpmType = 1, reRun = T, mlen = 40, unmergedDepth)
{

        mtB    = tB
        mnB    = nB
        mtSum  = tSum
        mnSum  = nSum

        #frequencyFile = paste(anaName, chrID, ".gw.dat", sep='.')
        #load(frequencyFile)
        startingPos = which(names(unmergedDepth$tB) %in% names(mtB))
        #startingPos = startingPos[1:(length(startingPos)-1)]
        genotypes   = maxRes$genotypes
        merge  <- function(aa, mlen=2){max(tapply(aa, as.integer(((1:length(aa)) +1) / mlen), sum))}
        print(genotypes)
        maxCN       = max(merge(genotypes, mlen=2))
        CN          = genotypes[maxRes$hidden.states*2-1] +  genotypes[maxRes$hidden.states*2]
        tab         = data.frame(cbind(unmergedDepth$nB,
                                       unmergedDepth$nSum,
                                       unmergedDepth$tB,
                                       unmergedDepth$tSum))
        pos         = as.numeric(names(mtB))
        ml = sapply(startingPos,
                    function(starting, tab, genotypes, tc, mlen)  {
                        ending = min(starting+mlen, nrow(tab))
                        maxll(tab[starting:ending, , drop=F], 
                              genotypes = genotypes, tc = tc) },
                    tab = tab, genotypes=genotypes, tc = maxRes$tc, mlen = mlen)

        #ml = c(ml,0)
        clean = 0.10
        # if the maxCN =12, the loci with cn > 6 will be defined as significantly AMP and filtered. 
        cnRange = c(1, maxCN/2)  
        selIncluded = rep(T, length(ml))
        for(eachGeno  in unique(maxRes$hidden.states))
        {
            aGenotype = maxRes$hidden.states == eachGeno
            selIncluded[aGenotype] = ml[aGenotype] > quantile(ml[aGenotype], prob=clean, type=1)
        }

        depthClean = 0.01
        outlayerBound = quantile(mtSum/mnSum, prob=1 - depthClean, type =1)

        bafClean = 0.01
        bafThresh = quantile(abs(mnB/mnSum - 0.5), prob = 1 - bafClean, type = 1)

        print(sprintf("length of ml %d and mtB %d", length(ml), length(mtB)))
        selIncluded = selIncluded & CN <= cnRange[2] & CN >= cnRange[1] & (mtSum / mnSum) < outlayerBound & abs(mnB/mnSum - 0.5) < bafThresh
        print(sprintf("Downsampled %f from %d",
                      sum(selIncluded)/length(selIncluded),
                      length(selIncluded)))


        nBMerged   = c()
        tBMerged   = c()
        nSumMerged = c()
        tSumMerged = c()
        mergeDegree = 6
        sameState <- function(states) { sum (states != states[1]) }
        for(idx in c(1:(sum(selIncluded)/mergeDegree)))
        {
            range = seq( (idx-1) * mergeDegree + 1,  idx * mergeDegree )
            if(sameState( maxRes$hidden.states[selIncluded][range]) == 0)
            #if(maxRes$hidden.states[selIncluded][idx*2-1] == maxRes$hidden.states[selIncluded][idx*2])
            {
                nBMerged   = c(nBMerged,   sum(mnB  [selIncluded][range]))
                tBMerged   = c(tBMerged,   sum(mtB  [selIncluded][range]))
                nSumMerged = c(nSumMerged, sum(mnSum[selIncluded][range]))
                tSumMerged = c(tSumMerged, sum(mtSum[selIncluded][range]))
            }
        }

        likelihoods = c()
        for(start in c(1:(length(DOA_range)-1)))
        {
            print(c(DOA_range[start], DOA_range[start +1]))
            res= segHMM2(rbind(nBMerged, nSumMerged, tBMerged, tSumMerged), 
                     length(genotypes)/2,
                     genotypes, 
                     maxiter = maxiter, 
                     DOA_range = c(DOA_range[start], DOA_range[start +1]),
                     method  = method, underate=underate, tpmType = tpmType)
            likelihoods = c(likelihoods, res$log.lik)
        }
        list(likelihoods=likelihoods, ml=ml, selIncluded = selIncluded)
}





plotPred  <- function(fluctuation, nB,nSum, tB, tSum, fullLabel, genotypes, maxRes, method = "Baum_Welch", len, chrID)
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
    rdTitle  = paste("Relative read depth",
                     "tc", roundTo(maxRes$tc, 3),
                     "DOA", roundTo(maxRes$DOA,3), sep="--")
    plot(label, tSum, main=rdTitle, xlab="loci", ylab="")
    plot(label, nSum, main=rdTitle, xlab="loci", ylab="")
    plot(label, tSum/nSum * (sum_di_n * maxRes$DOA/sum_di_t), main=rdTitle,
         xlab="loci", ylab="", ylim = c(0,5))
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

    legend("topleft", inset=.05, title="copy number", legend=c(1:10), fill=c(1:10))


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

segHMM2  <- function(allelicDepth, k, genotypes,  tc = 0.5, fix_tc = F, diag.prob = 0.9, maxiter = 10, eps = 0.0000001, DOA = 1.5, print.info = FALSE, method="Baum_Welch", DOA_range=c(1, 1.5), underate=1, tpmType =1)
{
    # Start the clock!
    ptm <- proc.time()
    effectSize = 5
    k  = length(genotypes)/2

    diag.prob = 0.5
    gamma     <- matrix((1 - diag.prob) / (k - 1), k, k)
    diag(gamma) <- diag.prob
    CN = genotypes[c(1:k)*2 - 1] + genotypes[c(1:k)*2 ]

    if(tpmType == 2)
    {
        print("Using copy number transition TPM")
        diag.prob = 0.5
        off.prob  = (1 - diag.prob)/ ( length(unique(CN)) -1)
        gamma     <- matrix(diag.prob, k, k)
        for (i in c(1:k))
        {
            for(j in 1:k)
            {
                cni  = genotypes[2*i -1] + genotypes[2*i ]
                cnj  = genotypes[2*j -1] + genotypes[2*j ]
                if(cni != cnj)
                    gamma[i,j]  <- off.prob
            }
        }

    }



    if(tpmType == 3)
    {
        print("Using symetric states TPM")
        diag.prob = 0.40

        # redundency defined as the states reverse complement. e.g. 1,0  -  0,1
        numOfNonRedundent = 0
        for (aCN in unique(CN))
        {
            num = sum(CN == aCN)
            if(num %% 2 !=0) numOfNonRedundent = (num + 1) / 2 + numOfNonRedundent
            else    numOfNonRedundent = (num) / 2 + numOfNonRedundent
        }
        off.prob  = (1 - diag.prob)/ ( numOfNonRedundent -1)
        aSymPenalty = 0.1
        for (i in c(1:k))
        {
            for(j in 1:k)
            {
                cni  = genotypes[2*i -1] + genotypes[2*i ]
                cnj  = genotypes[2*j -1] + genotypes[2*j ]
                if(cni > 3 && cni == cnj)
                {
                    highCNPenalty = 0.90
                }
                else
                {
                     highCNPenalty = 1
                }

                if(genotypes[2*i -1] ==  genotypes[2*j ] & genotypes[2*j -1] == genotypes[2*i ])
                {
                    # higher quatile, the closer to off.prob
                    #gamma[i,j]  <- (diag.prob - off.prob) / quantile + off.prob
                    #gamma[i,j]  <- gamma[i,j] * highCNPenalty
                    #gamma[i,j]  <- diag.prob * highCNPenalty * aSymPenalty
                    gamma[i,j]  <- off.prob  * 2
                    
                }
                else
                {
                    gamma[i,j]  <- off.prob  
                }
                if(i == j)
                {
                    gamma[i,j]  <- diag.prob * highCNPenalty
                }
            }
        }
    }




    initP = as.double(rep(-log(k), k))
    numobs = dim(allelicDepth)[2]
    ifSigAmp = rep(0, numobs)
    prob = rep(0.0, numobs * k)
        res <-
            .C(method,
               as.integer(numobs),
               as.integer(length(genotypes)/2),
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
               ifSigAmp = as.integer(ifSigAmp),
               as.double(underate),
               as.logical(print.info),
               #as.double(g_rcov),
               PACKAGE = "sCNAphase"
               )

        res$hidden.states <- res$hidden.states + 1
        res$gamma <- matrix(res$gamma, nr = k)

        # Stop the clock
        print(proc.time() - ptm)
        res$depths              = c()
        res$filtered.cond.probs = c()
        res
}
calcGC  <-  function(chrID, seqFile, winSize, pos)
{
    print(length(pos))
    range = cbind(pos - winSize, pos + winSize)
    fa <- open(FaFile(seqFile))
    selRange = GRanges(chrID, IRanges(range[,1], range[,2]) )
    dna <- scanFa(fa, param=selRange)
    gcContect = sapply(as.character(dna), function(seq){ (1-nchar(gsub("[GCgc]", "",seq)) / nchar(seq)) * 100})
    gcContect
}






inferCNA  <- function(anaName, nPrefix, tPrefix, chroms, doPhase = T, forceRead=F, maxCopyNum=11, mlen = 30, maxiter = 5, doFilter = F, ploidy = c(1.0, 1.4), fakeNormal = F, underate = 1, tpmType=1, runHMM = T, allelicMapability = F, reRun = T, method = "Baum_Welch")
{
    preprocess="phased"
    g_rcov = 0.5
    chrID = "W"

    #************************************************************************
    #          Possible genotypes 
    #************************************************************************
    startPos = 0
    gap      = 100000 # in bps
    label    = c()
    breakPoints =  startPos
    tB = c()
    nB = c()
    tSum = c()
    nSum = c()
    
        gcCombine = c()
        depthCombine = c()

    forceLoad = T
    if(!forceLoad)
    {
       load(paste(anaName, chrID, ".gw.dat", sep='.'))
    }else
    {
        for (AChrID in chroms)
        {
            res = loadFreq( ana = anaName,
                           doflip = doPhase,
                           AChrID, 
                           merge="perSite",
                           nPrefix = nPrefix,
                           tPrefix = tPrefix,
                           forceRead = F,
                           allelicMapability = allelicMapability
                   )
            tB       = c(tB   ,  res$tB  )
            nB       = c(nB   ,  res$nB  )
            tSum     = c(tSum ,  res$tSum)
            nSum     = c(nSum ,  res$nSum)
            notNA=!is.na(tB)
            tB=tB[notNA]
            tSum=tSum[notNA]
            nB=nB[notNA]
            nSum=nSum[notNA]

            label    = c(label,  as.numeric(names(res$tB)) + startPos + gap)
            startPos = max(label)
            breakPoints = c(breakPoints, startPos)
            if(fakeNormal)
            {
                baseDir = "/panfs/home/wenhanchen/work/DATA/humanAssembly/HG19/"
                print(sprintf("%s/chr%d.fa",baseDir, AChrID))
                xx = calcGC  (AChrID, sprintf("%s/chr%d.fa",baseDir, AChrID), winSize = 500, pos = as.numeric(names(res$tB)) )
                yy = res$tSum
                gcCombine    = c(gcCombine, xx)
                depthCombine = c(depthCombine, yy)
            }


        }
        names(tB)   = (label)
        names(nB)   = (label)  
        names(tSum) = (label)
        names(nSum) = (label)


      
        if(fakeNormal)
        {
            fit2    <- lm(depthCombine~poly(gcCombine,7))
            expectedDepth  <- predict(fit2)

            save(file="gc.dat", tB,tSum,nB,nSum, fit2, gcCombine, depthCombine)

            pdf(file=sprintf("%s.gcEffect.pdf", anaName))
            plot(gcCombine, predict(fit2), type="p", col="red", lwd=1)
            dev.off()
            print("Generating faked read depth that reflects the GC bias")

              nB   = nB/nSum * expectedDepth +1
              nB   = 0.5 * expectedDepth +1
              nSum = expectedDepth +1
              sel = nB <1 | nSum <1
              nB[sel] = 0.5
              nSum[sel]= 1

            sel   = gcCombine > 40 & gcCombine < 60
            tB    = tB   [sel]
            nB    = nB   [sel]
            tSum  = tSum [sel]
            nSum  = nSum [sel]
            label = label[sel]
            print("depth generated")

        }

        save(tB,nB,tSum,nSum, breakPoints, file=paste(anaName, chrID, ".gw.dat", sep='.'))
        ADepth <- function(tB, tSum, nB, nSum, breakPoints)
        {
            me <- list(
                tB= tB,
                tSum = tSum,
                nB = nB, 
                nSum = nSum, 
                breakPoints = breakPoints
            )
            ## Set the name for the class
            class(me) <- append(class(me),"ADepth")
            return(me)
        }

        unmergedDepth = ADepth(tB, tSum, nB, nSum, breakPoints)

    }

    #************************************************************************
    #          Possible genotypes 
    #************************************************************************
    genotypes = genGenotypes(maxCopyNum) # all the possible genotypes with copy number < maxCopyNum
    genotypes = c(0, 0, genotypes)       # adding doulbe deletion


    #************************************************************************
    #                 Merging (state shifting detection)
    #************************************************************************
    par(cex = 1.2)
    par(cex.axis = 1)
    par(cex.lab = 1)

    outfile=paste(doPhase, anaName, chrID, preprocess, ".pdf", sep=".")
    pdf(file=outfile, width=24, height=6, title=outfile)





    #************************************************************************
    #                 Merging (state shifting detection)
    #************************************************************************
    indexing<-function(idx) {as.numeric(names(idx))}
    mergedData = mergeUsingML(data.frame(cbind(nB,nSum, tB, tSum)), 
              genotypes,
              tc=0.99,
              len = mlen)  #  Merge the allelic depth and mind the allelic depth consistency.


    if(doFilter == T)
    {
        mergedData = filtering (mergedData, -80)
    }

    fluctuation = mergedData[["fluctuation"]]
    fullLabel   = mergedData[["fullLabel"]]
    mergedDepth = mergedData[["depths"]]

    tB   = ceiling(mergedDepth$tB)
    nB   = ceiling(mergedDepth$nB)
    tSum = ceiling(mergedDepth$tSum)
    nSum = ceiling(mergedDepth$nSum)

    names(tB) = names(fluctuation)

  

    #************************************************************************
    #                 Remove germ-line CNVs
    #************************************************************************
    removeCNVs = T
    if(removeCNVs)
    {
        cleanBAF = 0.01
        bafThresh = quantile(abs(nB/nSum - 0.5), prob = 1 - cleanBAF, type = 1)
        selNonCNVs   = abs(nB/nSum - 0.5) < bafThresh
        tB    = tB   [selNonCNVs]
        nB    = nB   [selNonCNVs]
        tSum  = tSum [selNonCNVs]
        nSum  = nSum [selNonCNVs]
        fluctuation = fluctuation [selNonCNVs]
        fullLabel = fullLabel[selNonCNVs,]
        label = label[selNonCNVs]
    }

 
    #************************************************************************
    #                 CN state prediction 
    #************************************************************************
    allRes = inferStates  (anaName, preprocess, chrID,  method = method, 
                           nB,nSum, tB, tSum,
                           genotypes, maxiter=maxiter, fullLabel = fullLabel, DOA_range = ploidy,
                           tpmType = tpmType,
                           underate = underate,
                           runHMM = runHMM,
                           mlen = mlen)

    extract <- function(aList, field){ aList[field] }
    logLiks = unlist( sapply(allRes, extract, field = "log.lik"))
    maxRes = allRes [[ which(max(logLiks) == logLiks)]]


    dataFile = paste("res", anaName, preprocess,"chr", chrID,"dat", sep=".")
    label = as.numeric(names(tB))
    if(reRun)
    {
        likelihoodEst = runForPloidy  (anaName, preprocess, chrID,  method = method, 
                           nB,nSum, tB, tSum,
                           maxRes = maxRes,
                           genotypes, maxiter=maxiter, fullLabel = fullLabel, DOA_range = ploidy,
                           tpmType = tpmType,
                           underate = underate,
                           runHMM = runHMM,
                           mlen = mlen, unmergedDepth = unmergedDepth)
        likelihoods  = likelihoodEst$likelihoods
        ml          = likelihoodEst$ml
        selIncluded = likelihoodEst$selIncluded

        # find the max
        maxIdx = which(likelihoods == max(likelihoods))[1]
        theDOA     = allRes[[maxIdx]]$DOA
        genotypes = maxRes$genotypes
        res = segHMM2(rbind(nB[selIncluded], nSum[selIncluded], tB[selIncluded], tSum[selIncluded]), 
                        length(genotypes)/2,
                        genotypes, 
                        maxiter = 1, 
                        DOA_range = c(theDOA, theDOA + 0.05),
                        method  = method, underate=underate, tpmType = tpmType)
        similiarity <- function(aRes, tStates ) {sum(aRes$hidden.states[selIncluded] == tStates)}
        trueMaxIdx = which.max(unlist(lapply(allRes, similiarity, tStates = res$hidden.states)))
        maxRes = allRes[[trueMaxIdx]]
        save(nB,nSum, tB, tSum, label, fullLabel, allRes,
             likelihoods, ml, selIncluded, file=dataFile)
    } else
    {
        save(nB,nSum, tB, tSum, label, fullLabel, allRes, file=dataFile)
    }




    cat(paste ("Degree of Amplification : ", maxRes$DOA, sep=":"))
    cat(paste ("Predict tumor cellularity : ", maxRes$tc, maxRes$log.lik, sep=":"))
    plotPred (fluctuation, nB,nSum, tB, tSum, fullLabel, genotypes, maxRes, len=mlen, chrID=chrID)

    depths = rbind(nB,nSum,tB,tSum)
    T = dim(depths)[2]

    prob_RD=rep(0.1, T)
    prob_AD=rep(0.1, T)

    geno = maxRes$genotypes
    cnStates = maxRes$hidden.states
    cns = geno[maxRes$hidden.states*2-1] +  geno[maxRes$hidden.states*2]
    tc = maxRes$tc
    DOA = maxRes$DOA

    res <- .C("emissionDist_Debug",
              depths = as.integer(as.vector(depths)),
               T = as.integer(T),
                geno = as.integer(geno),
                k = as.integer( length(geno)/2),
                cnStates = as.integer(cnStates),
                tc = as.double(tc),
                DOA = as.double(DOA), 
                prob_RD = as.double(prob_RD),
                prob_AD = as.double(prob_AD),
                underate = underate,
                PACKAGE="sCNAphase"   
                )

    col = rep(1, length(cns))
    col [cns>8] = 2
    plot(res$prob_AD, pch = 4, col = col)
    plot(res$prob_RD, pch = 4, col = col)
    save(file="prob.dat", res)
    dev.off()
}


