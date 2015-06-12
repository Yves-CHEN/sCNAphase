/************************************************************************
 *
 *  Author : Wenhan
 *  Date   : Sep 20 2014
 *  Decription : This is for getting amount of reads in a specific region 
 *    in a bam file. Dependency is on htslib which contains the function 
 *    for parsing and reading bam files.
 *
 * ************************************************************************/

#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <string>
#include "hts/sam.h"  
#include "R.h"
  
class myexception
{

public:
    std::string err;
    myexception(std::string message)
    {
        err = message;
    }
};

int qcFilter(bam1_t *b, int thresh)
{
    printf("Q : %d \t ", b->core.qual);
    if(b->core.qual < thresh)
        return 1;
    return 0;
}

class BamReader
{
    private:
        bam_hdr_t *header;
        samFile   *samStream;
        BamReader();            // forbidden
        bam1_t* oneReadBuf;
        hts_idx_t *idx; // load index
        inline int qcFilter(bam1_t *b, int thresh)
        {
            if(b->core.qual < thresh)
                return 1;
            return 0;
        }
    public:
        BamReader(const char* fileName)
        {
            this->idx= NULL;
            // create sam file stream
            this->samStream = sam_open(fileName, "r");
            if (this->samStream  == NULL) {
                Rprintf("[error] failed to open \"%s\" for reading\n", fileName);
                throw myexception("fail to read.");
            }
            // Loading the indexing for bam
            this->idx = bam_index_load(fileName);
            if (this->idx == NULL) { // index is unavailable
                Rprintf("[samview] random alignment retrieval only works for indexed BAM   or CRAM files.\n");
                throw myexception("fail to read.");
            }
            // Reading header of bam
            this->header = sam_hdr_read(samStream);
            if ( this->header == NULL) {
                Rprintf("[samview] fail to read the header from \"%s\".\n", fileName);
                throw myexception("fail to read.");
            }
            // Initiate a buf for loadin each read.
            oneReadBuf = bam_init1();

        }
        ~BamReader()
        {
            if(idx!=NULL)
            {
                bam_destroy1(oneReadBuf);
                hts_idx_destroy(idx);
            }
        }
        int read(char* aRegion, int qcThresh) /// in the format of chr8:100-1000
        {
            int result = -1;
            int count  = 0;
            hts_itr_t *iter = sam_itr_querys(idx, header, aRegion);
            if (iter == NULL) { // reference name is not found
                Rprintf("[samview] region \"%s\" specifies an unknown reference name.   Continue anyway.\n",
                        aRegion);
                throw myexception("Region specifier in the wrong format");
            }

            while ((result = sam_itr_next(samStream, iter, oneReadBuf)) >= 0) 
            {
                if(!qcFilter(oneReadBuf, qcThresh) )
                    count++;
            }
            hts_itr_destroy(iter);
            if (result < -1) {
                Rprintf("[samview] retrieval of region \"%s\" failed due to truncated   file or corrupt BAM index file\n", aRegion);
                throw myexception("Incompleted bamfile");
            }
            return count;
        }

};

extern "C"
{

    void countReads(char** bamF, char** chrID, int* start, int* windowSize, int* num, int* qcThresh, int* countVec, bool* readFailure)
    {
        try
        {
            BamReader reader(*bamF);
            for(int i=0; i < *num; i ++)
            {
                char theRegion[10000];
                sprintf(theRegion, "%s:%d-%d", *chrID, start[i], start[i] + windowSize[i]);
                countVec[i] = reader.read(theRegion, *qcThresh);
            }
       }
       catch(myexception ex)
       {

           *readFailure = true;
           return;
       }
    }
}
  
 
