
displayParas  <- function()
{

}

# mlen =30 (default). Merge every 30 loci.
# phasedADs.          A data frame of allelic depth in the order of (nB,nSum, tB, tSum)
# genotypes.          Pairs of maternal and paternal alleles.
# tc=1,               Purity of sample. By default pure tumor, tc =1
# thresh  = -40       Threshold of the likelihood test for allelic frequency consistency
#
mergeADs  <- function(mlen = 30, phasedADs, genotypes, tc=1, thresh  = -40)
{
    mergedData=mergeUsingML(phasedADs, 
              genotypes,
              tc=tc,
              len = mlen)
    mergedData = filtering (mergedData, thresh = thresh )
    fluctuation = mergedData[["fluctuation"]]
    fullLabel   = mergedData[["fullLabel"]]
    mergedDepth = mergedData[["depths"]]
    tB   = mergedDepth$tB   
    nB   = mergedDepth$nB   
    tSum = mergedDepth$tSum 
    nSum = mergedDepth$nSum 
    names(tB) = names(fluctuation)
    mADs = list(tB,nB,tSum, nSum, fullLabel)
    names(mADs) = c("tB", "nB", "tSum", "nSum", "fullLabel")
    mADs

}
genGenotypes <- function(copyNum = 9)
{
    crash  <- function( tab )
    {
        tab=t(tab)
        dim(tab) = c(1, dim(tab)[1] * dim(tab)[2])
        tab
    }
    genotypes=c()
    for (cn in 1:copyNum)
    {
        genotypes=c(genotypes, crash(cbind(cn:0, 0:cn)))
    }
#    if(extra > 0)
#    {
#        genotypes=c(genotypes, crash(cbind(6:0, 0:6)))
#    }
    genotypes = c(genotypes)  

#    if(copyNum > 10)
#    {
#    sigAmp = c(1,19,19,1,1,29, 29, 1,1,39,39,1,1,49,49,1,1,59,59,1,1,69,69,1,1,79,79,1,1,89,89,1,1,99,99,1,1,109,109,1, 2, 108, 108, 2, 1,119, 119, 1, 2,118,118,2) # copy number of 10,20,...,120
#    genotypes = c(genotypes, sigAmp)
#    }


    genotypes
}


readChroms  <- function (dataDir, nPrefix, tPrefix, anaName="ana", chrIDs)
{
    data(chromLen) # chromLen from R package
    tB    = c()
    nB    = c()
    tSum  = c()
    nSum  = c()
    label = c()
    startPos = 0
    for (AChrID in chrIDs)
    {
        res = loadFreq(dataDir,
                        ana = anaName,
                        doflip = T,
                        AChrID,  merge="perSite",
                        nPrefix = nPrefix,
                        tPrefix = tPrefix
        )
        tB       = c(tB   ,  res$tB  )
        nB       = c(nB   ,  res$nB  )
        tSum     = c(tSum ,  res$tSum)
        nSum     = c(nSum ,  res$nSum)
        label    = c(label,  as.numeric(names(res$tB)) + startPos )
        startPos = startPos + chromLen[AChrID, 2]
    }
    names(tB)   = (label)
    names(nB)   = (label)  
    names(tSum) = (label)
    names(nSum) = (label)
    notNA=!is.na(tB)            
    tB=tB[notNA]
    tSum=tSum[notNA]
    nB=nB[notNA]
    nSum=nSum[notNA]
    ADs = cbind(nB,nSum, tB, tSum)
    #colnames(ADs) = c("nB", "nSum", "tB", "tSum")
    #as.data.frame(ADs)
    data.frame(ADs)

}



inferStatesWrap  <- function(anaName = "ana", chrIDs, mergeDepth, genotypes)
{
    g_rcov = 1
    nB   = mergeDepth[["nB"  ]]
    nSum = mergeDepth[["nSum"]]
    tB   = mergeDepth[["tB"  ]]
    tSum = mergeDepth[["tSum"]]
    chrID = paste(chrIDs[1], chrIDs[length(chrIDs)], sep="-")
    maxRes = inferStates  (anaName, preprocess = "mergeSite", chrID, runHMM = T, method = "Baum_Welch", 
           nB,nSum, tB, tSum,
           genotypes, maxiter=30, g_rcov=g_rcov)

    maxRes
}



