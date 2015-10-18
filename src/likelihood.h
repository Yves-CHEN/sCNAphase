
//#include "hmm/viterbi_maxll.h"



extern "C"
{

    void calcLikelihood ( int* depths, int* seq_len, int* genotypes, int* _k,
                            double* _tc, double* _DOA, double* TPM, double* pi, double* underate, double* logLike, int* states) 
    {

        int    k = *_k;
        int    T = *seq_len;
        double tc = *_tc;
        double DOA = *_DOA;

        vector< vector<double> > tpm(k, vector<double>(k)),
                yll(k, vector<double>(T)), alpha(k, vector<double>(T)),
                beta(k, vector<double>(T));



        for(int i = 0; i < k; i++)
            for(int j = 0, offs = 0; j < k; j++, offs += k)
                tpm[i][j] = log(TPM[offs + i]);

        double del = 0.001;

        emissionDist(depths, T, genotypes, k, tc, DOA, yll, del, *underate);
        bool faster     = true;
        bool print_info = false;
        bool doFilter   = false;
        double* filter  = NULL;  // filter is an useless parameter. Needs to be deleted!
        double log_lik = forward_backward(yll, tpm, pi, doFilter, filter, alpha, beta, faster,print_info);
        *logLike = log_lik;


 
//Rprintf("Predicted tumor cellularity cc = %f. when ll = %e \n", tc, log_lik);
          viterbi(yll, tpm, pi, states);

    }
}


 

   

