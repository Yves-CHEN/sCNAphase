#**************************************************************************************************************************
#  For each of the runs, sCNAphase requires a unique anaName, so that the output filenames with this anaName, 
#  do not overwrite each other. The default output file is a "res.anaName.phased.chr.W.dat" file. This contains
#  the inputs, the setting and the outputs from the run. The following output functions based on the "*.dat" files
#  can generate more readable output, including the 
#    *  "*.csv" segmentation files (function genSegFile),
#    *  ".vcf"  files (function genVCFFile), 
#    *  "*.pdf" the d.SKY plots (function produceDSKY), 
#    *  "*.pdf" the copy number plots (function produceCytoPlot).
#  
#  genSegFile  <- function(anaList = c("5","20","40","60", "80", "95", "p"), outdir = "./")
#  
#  genVCFFile  <- function( anaName, vcfFile, chrID,  outDir = "./")
#  
#  produceCytoPlot  <- function(anaList = c("5","20","40","60", "80", "95", "p"), col = 2, outDir = "./")
#  
#  produceDSKY      <- function(anaList = c("5","20","40","60", "80", "95", "p"), outDir = "./")
#  
# ******** genseg is employed by genSegFile, produceDSKY, produceCytoPlot to generates copy number segments.
#  genseg      <- function (copyNum, mCN, label, fullLabel, dir, chrID, distThresh = 1000000 )
#  

# Plot Pred generates a series of plots for debuging purposes only.
# plotPred    <- function(fluctuation, nB,nSum, tB, tSum, fullLabel, genotypes, maxRes, method = "Baum_Welch", len, chrID)
# 




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

    genoLabel = c("A", "AB", "ABB", "ABBB", "ABBBB" )
    #genoLabel = c("M", "MP", "MPP", "MPPP", "MPPPP" )
    for(ii in 1:length(genoLabel))
    {
        text(-40 , ii , labels = genoLabel[ii], font=5, offset=10.5 , adj=1)
    }
    legend("topleft", inset=.05, title="Methods", c("Climat","mine"), fill=c(2,3))
}


# 1M
genseg <-function (copyNum, mCN, label, fullLabel, dir, mergingSel, chrID, distThresh = 1000000 )
{

    # This removed the PHFs with excessive merging errors.
    copyNum    = copyNum[ mergingSel]
    mCN        = mCN[ mergingSel]
    label      = label[ mergingSel]
    fullLabel  = fullLabel[ mergingSel, ]

    preState = -1
    seg = c()
    theCN = c()
    names(mCN) = names(copyNum)

    successiveDiff  <-  function(vecData) 
    { 
        theDiff = c()
        if(length(vecData) < 2) {
            print("Error in sucessiveDiff: vector len <2 ")
        } else {
           theDiff= vecData[2:length(vecData)] - vecData [1:(length(vecData) -1)]
        }
        theDiff

    }
    distBreak = successiveDiff(label) > distThresh
    cnBreak   = successiveDiff(copyNum) != 0
    mcnBreak  = successiveDiff(mCN) != 0
    startAt = fullLabel[c(T, cnBreak | distBreak | mcnBreak), 1]
    endAt   = fullLabel[c(cnBreak | distBreak | mcnBreak, T), 2]

    if(length(fullLabel[fullLabel[,2]<0, 2]) > 0)
    {
        print(label[fullLabel[,2]<0])
        print ("negtive full label")
    }
    if(length(label[label<0]) > 0)
    {
        print ("negtive label")
    }
    tab = as.data.frame(cbind(rep(chrID, length(copyNum[startAt])), startAt, endAt, as.integer(copyNum[c(T, cnBreak | distBreak | mcnBreak)]), as.integer(mCN[c(T, cnBreak | distBreak | mcnBreak)])))

    colnames(tab) = c("chr", "start","end","CN", "mCN")
    tab$chr   = as.numeric(as.character(tab$chr))
    tab$start = as.numeric(as.character(tab$start))
    tab$end   = as.numeric(as.character(tab$end))
    tab$CN    = as.numeric(as.character(tab$CN))
    tab$mCN   = as.numeric(as.character(tab$mCN))
    tab
}




# If the probability of merging error, ml, is calculated, then remove the sites with high merging error likelihood.
# Otherwise, it does not perform screening.
genSegFile  <- function(anaList = c("5","20","40","60", "80", "95", "p"), mlRemoveLevel = 0.04, outdir = "./", ifload = T)
{

    mSmooth <- function(pos, win=5)
    {
        len=as.integer(win/2)
        startAt = max(1, pos-len)
        endAt  = min(length(CN), pos + len)
        mm=median(CN[startAt:endAt])
        mm
    }

    rescaleTC  <- function(tc, rcov, copyNum)
    {
        factor = sum(copyNum) / (2 * length(copyNum))
        theRatio = tc / (1-tc) * factor
        d_tc = theRatio / (1+theRatio)
        d_tc
    }



    extract  <- function()
    {   
        genotypes = maxRes$genotypes
        states    = maxRes$hidden.states
        DOA       = maxRes$DOA
        tc        = maxRes$tc + 0.01 
        copyNum = genotypes[states * 2] + genotypes[states * 2 -1]
        rcov=c()
        for (xx in 1:12)
        {
            sel   = copyNum == xx
            col   = c(col, rep(2, length(states[sel])))
            ratio = 2/(2 + (xx -2) * tc)
            rcov  = c(rcov, tSum[sel] * ratio/nSum[sel])
        }
        d_tc = rescaleTC(tc, mean(rcov), copyNum)
        c(tc, d_tc, maxRes$DOA*2)
    }


    if(! file.exists(outdir))
        dir.create(outdir, recursive = T)

    for( eachAna in anaList )
    {
        datFile = sprintf("res.%s.phased.chr.W.dat", eachAna)
        if(! file.exists(datFile))
            stop(sprintf("Excepting dat file : %s", datFile))
        if(ifload) { load(datFile)  }
        breakPoints = unmergedDepth$breakPoints
        names(breakPoints) = chroms
        clean = 0.01
        pos=as.numeric(names(tB))
        sel=fullLabel[, 1] %in% pos 
        fullLabel = fullLabel[sel,]
        pos=as.numeric(names(tB))
        gwTab = as.data.frame(matrix(NA, ncol=5, nrow = 0))
        colnames(gwTab) = c("chr", "start","end","CN", "mCN")

        if ( exists('ml') )
        {
            mlThresh = quantile(ml, type=1, prob = mlRemoveLevel)
        }

        pdf(file = paste(eachAna, "pdf", sep='.'), width = 18)


        print(extract())



        for (idx in 1:length(chroms))
        {
            chrID = as.character(chroms[idx])
            sel = breakPoints[idx] < pos & pos <breakPoints[idx +1]


            if ( exists('ml') )
            {
                plot(  tB[sel]/ tSum[sel], col = (ml[sel] < mlThresh) +1, main = idx)
            }

             
            FL=as.numeric(fullLabel[sel,])
            dim(FL)  = dim(fullLabel[sel,])
            #gap   = 100000
            FL[,1] = as.numeric(fullLabel[,1][sel]) - breakPoints[chrID] -gap
            FL[,2] = as.numeric(fullLabel[,2][sel]) - breakPoints[chrID] -gap
            label=pos[sel] - breakPoints[chrID]
            geno = maxRes$genotypes
            bCN = geno[maxRes$hidden.states[sel] *2]
            aCN =  geno[maxRes$hidden.states[sel] *2 -1]
            mCN = bCN
            msel = bCN > aCN
            mCN[msel] = aCN[msel]
            CN = bCN + aCN
            CN = sapply(c(1:length(CN)), mSmooth)
            if(length(CN) > 2)
            {
                if(exists('ml'))
                {
                    tab = genseg(CN, mCN, label, FL, ml[sel] > mlThresh, chrID=chrID, dir=outdir)
                } else 
                {
                    tab = genseg(CN, mCN, label, FL, CN > -1, chrID=chrID, dir=outdir)
                }
                #print(head(tab))
            } else
            {
                print(sprintf("chr%d is gone", chrID ))
            }
            gwTab = rbind(gwTab, tab)
        }
        dev.off()
        
        gwFile = sprintf("%s/cna.%s.tumor.W.csv", outdir, eachAna) 
        print(sprintf("Generating file %s", gwFile))
        write.table(gwTab, file=gwFile, sep="\t", quote=F, row.names=F)
    }
}







produceCytoPlot  <- function(anaList = c("5","20","40","60", "80", "95", "p"), outDir = "./", col =2, removeLevel = 0.005)
{
    exist = require(viewGenome)
    if(!exist)
    {
        stop("The package viewGenome is require. Find it from github.")
    }

    mSmooth <- function(pos, win=7)
    {
        len=as.integer(win/2)
        startAt = max(1, pos-len)
        endAt  = min(length(CN), pos + len)
        mm=median(CN[startAt:endAt])
        mm
    }
    if(! file.exists(outDir))
    {
        dir.create(outDir, recursive = T)
    }
        for( eachAna in anaList )
    {
        datFile = sprintf("res.%s.phased.chr.W.dat", eachAna)

 
        if(! file.exists(datFile))
            stop(sprintf("Excepting dat file : %s", datFile))
        load(datFile)
        
        mlThresh = quantile(ml, type=1, prob = removeLevel)

        breakPoints = unmergedDepth$breakPoints
        names(breakPoints) = chroms
        clean = 0.01
        pos=as.numeric(names(tB))
        sel=fullLabel[, 1] %in% pos 
        fullLabel = fullLabel[sel,]
        pos=as.numeric(names(tB))
        gwTab = as.data.frame(matrix(NA, ncol=4, nrow = 0))
        colnames(gwTab) = c("chr", "start","end","CN")

        print(sprintf("produceing : %s/%s.per.chrom.pdf", outDir, eachAna))
        pdf(sprintf("%s/%s.per.chrom.pdf", outDir, eachAna), width = 16)

        for (idx in 1:length(chroms))
        {
            chrID = as.character(chroms[idx])
            sel = breakPoints[idx] < pos & pos <breakPoints[idx +1]
            FL=as.numeric(fullLabel[sel,])
            dim(FL)  = dim(fullLabel[sel,])
            #gap   = 100000
            FL[,1] = as.numeric(fullLabel[,1][sel]) - breakPoints[chrID] -gap
            FL[,2] = as.numeric(fullLabel[,2][sel]) - breakPoints[chrID] -gap
            label=pos[sel] - breakPoints[chrID]
            geno = maxRes$genotypes
            bCN = geno[maxRes$hidden.states[sel] *2]
            aCN =  geno[maxRes$hidden.states[sel] *2 -1]
            mCN = bCN
            msel = bCN > aCN
            mCN[msel] = aCN[msel]
            CN = bCN + aCN
            CN = sapply(c(1:length(CN)), mSmooth)
            if(length(CN) > 2)
            {
                tab = genseg(CN, mCN, label, FL, ml[sel] > mlThresh, chrID = chrID,  dir=outDir)
                myplot(x = tab[,2],
                       y = tab[,4], chr=chrID, col = col,
                       graybg = T, xlab = "Position", ylab = "Copy number",
                       ylim = c(0,max(tab[,4])))
            } else
            {
                print(sprintf("chr%d is gone", chrID ))
            }
        }
        dev.off()
    }
}


# If the probability of merging error, ml, is calculated, then remove the sites with high merging error likelihood.
# Otherwise, it does not perform screening.
produceDSKY  <- function(anaList = c("5","20","40","60", "80", "95", "p"), outDir = "./", removeLevel = 0.005, ifload = T)
{
    exist = require(viewGenome)
    if(!exist)
        stop("The package viewGenome is require. Find it from github.")
    mSmooth <- function(pos, win=7)
    {
        len=as.integer(win/2)
        startAt = max(1, pos-len)
        endAt  = min(length(CN), pos + len)
        mm=median(CN[startAt:endAt])
        mm
    }
    if(! file.exists(outDir))
        dir.create(outDir, recursive = T)

    for( eachAna in anaList )
    {
        datFile = sprintf("res.%s.phased.chr.W.dat", eachAna)
        if(! file.exists(datFile))
            stop(sprintf("Excepting dat file : %s", datFile))
        if(ifload) { load(datFile) }
        
        if(exists("ml"))
        {
            mlThresh = quantile(ml, type=1, prob = removeLevel)
        }


        breakPoints = unmergedDepth$breakPoints
        names(breakPoints) = chroms
        clean = 0.01
        pos=as.numeric(names(tB))
        sel=fullLabel[, 1] %in% pos 
        fullLabel = fullLabel[sel,]
        pos=as.numeric(names(tB))
        gwTab = as.data.frame(matrix(NA, ncol=4, nrow = 0))
        colnames(gwTab) = c("chr", "start","end","CN")

        print(sprintf("produceing : %s/%s.d.SKY.pdf", outDir, eachAna))
        pdf(sprintf("%s/%s.d.SKY.pdf", outDir, eachAna), width = 16)
        par(mfrow=c(3,8),  omi=c(0.01,0.01,0.01,0.01), plt=c(0,0.99,0,0.97))

        #for (chrID in chroms)

        for (idx in 1:length(chroms))
        {
            chrID = as.character(chroms[idx])
            sel = breakPoints[idx] < pos & pos <breakPoints[idx+1]
            FL=as.numeric(fullLabel[sel,])
            dim(FL)  = dim(fullLabel[sel,])
            gap   = 100000
            FL[,1] = as.numeric(fullLabel[,1][sel]) - breakPoints[chrID] -gap
            FL[,2] = as.numeric(fullLabel[,2][sel]) - breakPoints[chrID] -gap
            label=pos[sel] - breakPoints[chrID]
            geno = maxRes$genotypes
            bCN = geno[maxRes$hidden.states[sel] *2]
            aCN =  geno[maxRes$hidden.states[sel] *2 -1]
            mCN = bCN
            msel = bCN > aCN
            mCN[msel] = aCN[msel]
            CN = bCN + aCN
            CN = sapply(c(1:length(CN)), mSmooth)
            if(length(CN) > 2)
            {

                if(exists("ml"))
                {
                    tab = genseg(CN, mCN, label, FL, ml[sel] > mlThresh,  chrID = as.numeric(chrID), dir=outDir)
                } else 
                {
                    tab = genseg(CN, mCN, label, FL, CN > -1, chrID = as.numeric(chrID), dir=outDir)
                }
                plotPloidy(as.numeric(chrID), tab[, 2:5 ])
            } else
            {
                print(sprintf("chr%d is gone", chrID ))
            }
        }
        dev.off()
    }
}


genVCFFile  <- function(anaName, vcfFile, chrID, outDir = "./")
{
    ifexist = require(VariantAnnotation)
    if(!ifexist)
        stop("Package VariantAnnotation is require. Please install it from bioconduct.")
    getInfo  <-  function(chrID)
    {
        datFile = sprintf("res.%s.phased.chr.W.dat", anaName)
        if(! file.exists(datFile))
            stop(sprintf("Excepting dat file : %s", datFile))
        load(datFile)  # This includes mlen, gap, maxRes, unmergedDepth
        breakPoints = unmergedDepth$breakPoints
        sel = breakPoints[chrID] < as.numeric(names(unmergedDepth$tB)) & 
                     as.numeric(names(unmergedDepth$tB)) < breakPoints[chrID+1]
        pAD = mapply(function(tB, tSum){list(c(tB, tSum-tB))}, 
                     unmergedDepth$tB[sel], unmergedDepth$tSum[sel])
        allPos   = as.numeric(unlist(lapply(strsplit((rownames(info(vcf))), split="[:_]"), 
                                       function(dat){c(dat[2])})))
        allAD  = sapply(c(1:length(allPos)), function(dat){list(c(-1,-1))})
        allpCN    = sapply(c(1:length(allPos)), function(dat){list(c(-1,-1))})
        allCN  = rep(-1, length(allPos))
        pos = as.numeric(names(tB))
        selAPos = (allPos %in% (as.numeric(names(unmergedDepth$tB)) [sel] - breakPoints[chrID] - gap)) & (  info(vcf)$INDEL == F )
        if(length(pAD) !=  sum(selAPos))
            stop(sprintf("%d != %d, could it be due to the wrong vcf file, or dat file", length(pAD), sum(selAPos)))
        allAD [selAPos] =  pAD
        sel   = breakPoints[chrID] < pos & pos <breakPoints[chrID+1]
        label = pos[sel] - breakPoints[chrID]  - gap
        geno  = maxRes$genotypes
        bCN = (geno[maxRes$hidden.states[sel] *2])    
        aCN = (geno[maxRes$hidden.states[sel] *2 -1])
        CN  = (bCN + aCN) 
        indices  = which(allPos[selAPos] %in% label)
        correspondingCN = CN[ label %in% allPos[selAPos] ]
        pCN = mapply(function(aCN, bCN){list(c(aCN, bCN))}, 
                 bCN[ label %in% allPos ], aCN[ label %in% allPos])
        indices = indices[1:(length(indices)-1)]
        for (inc in 1:mlen)
        {
            allCN[selAPos][indices + inc -1]  = correspondingCN[1:length(indices)]
            allpCN[selAPos][indices + inc -1] = pCN[1:length(indices)] 
        }
        copyNum = cbind(allAD, allCN, allpCN)
        colnames(copyNum) = c("pAD", "CN", "pCN")
        as.data.frame(copyNum)
    }
    vcf<-readVcf(vcfFile,"hg19")
    #----
    # The annotation in header to be appended to a vcf.
    moreAnnotation = rbind( c("2", "Integer", "phased allelic depth"), 
                         c("1", "Integer", "Copy number"), 
                         c("2", "Integer", "allelic copy number refering to phased allelic depth"))
    rownames (moreAnnotation) = c("pAD", "CN", "pCN")
    colnames (moreAnnotation) = colnames(info(header(vcf)))
    #----
    # The data to be appended to a vcf.
    extraData = getInfo(chrID)
    #----
    # This enables the rbind of the DataFrame object.
    appendDataFrame <- function ( base, appending, binding) 
    {
        DataFrame(binding(as.data.frame(base), appending))
    }
    #----
    # The annotation in the header of a vcf file and the data of vcf.
    # 1) These info(vcf) is setable. So that I can nodify the vcf object and create a modified vcf.
    info(header(vcf)) <- appendDataFrame(info(header(vcf)), moreAnnotation, binding = rbind)
    info(vcf)[["pAD"]]   <-  IntegerList (extraData$pAD)
    info(vcf)[["CN"]]    <-  as.numeric(extraData$CN)
    info(vcf)[["pCN"]]   <-  IntegerList(extraData$pCN)
    outDir = sprintf("%s/converted/", outDir)
    if(!file.exists(outDir))
    {
        dir.create(outDir, recursive = T)
    }
    print(sprintf("Generating converted vcf file to dir: %s", outDir))
    writeVcf(vcf, filename=sprintf("%s/%s", outDir, vcfFile))
}


