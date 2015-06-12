// Implementation of the forwards-backwards algorithm for
// parameter estimation of Gaussian density HMM.

#include <R.h>
#include <Rmath.h>
#include <vector>
#include <map>
#include <algorithm>
#include "hmm/viterbi_maxll.h"
#include "hmm/FB_maxll.h"
#include "bamUtils/bamReader.h"
#include "optimization/mll.h"

using namespace std;

/* --------------------------------------------
 *           Global Varables
 --------------------------------------------*/
const bool printWarning = false;

//extern  "C"   void countReads(char* bamF, char* chrID, int* start, int* windowSize, int* num, int* qcThresh, int* countVec);
    //extern "C" void countReads(char** bamF, char** chrID, int* start, int* windowSize, int* num, int* qcThresh, int* countVec);
 
extern void maxll (int* tab, int* num_n, int* genotypes, int* num_k,double* cc, int* state, double* max_ll);
extern void maxl_improv (int* tab, int* num_n, int* genotypes, int* num_k, double* cc, int* state, double* max_ll);



extern "C"
{
//        genotype(k * 2, 0); //  <maternal alleles,  paternal alleles> for kth state(genotype).
//        c                   // c for tumor cellularity.
//
    void viterbi_train(int *seq_len, int *_k, int* depths,
                             int* genotypes, double* tc, double *TPM, double *pi, int *maxiter,
                             double *eps, double *_log_lik, double *filter,
                             int *hidden_states,bool* fix_tc, bool *print_info)

    {
        // Initialization
        
        bool    verbal = true;
        int k = *_k, T = *seq_len;
        vector< vector<double> > tpm(k, vector<double>(k)),
            yll(k, vector<double>(T));

        /// preserver last estimation including the tpm, hidden_states and tc.
        vector< vector<double> >* p_tpm, * p_yll;
        vector<double> old_pi(k, 1 /double(k));
        double p_tc=*tc;

        for(int i = 0; i < k; i++)
            for(int j = 0, offs = 0; j < k; j++, offs += k)
                tpm[i][j] = log(TPM[offs + i]);
        *_log_lik = R_NegInf;
        int iter = 0;

        Rprintf("\tc = %f\tseq_len = %d\tnumOfStates = %d\n", *tc, *seq_len, *_k);

        // A tuple of four values for each of the loci is in the order of:
        //   
        //           |   AD  |   RD
        //   normal  |  mi_n |  di_n
        //   tumor   |  mi_t |  di_t
        double cellularity_iter = *tc;
        p_tpm = new vector<vector<double> > (tpm);
        p_yll = new vector<vector<double> > (yll);

        while(iter < *maxiter)
        {
            double avg = 0;
            // preserving
            delete(p_yll);
            p_yll = new vector<vector<double> >(yll);
            for(int m = 0; m < T; m++)
            {
                for(int i = 0; i < k; i++)
                {
                    int g_i[2] = {0, 0};
                    g_i[0] = genotypes[i * 2];
                    g_i[1] = genotypes[i * 2 +1];
                    yll[i][m] = logProb(depths[m*4], depths[m*4 +1], depths[m*4 +2], depths[m*4 +3], 
                            cellularity_iter, g_i);
                }
            }

            // using viterbi to find best possible hidden states for para estimation
            double log_lik =  viterbi(yll, tpm, pi, hidden_states);
            Rprintf("** log_lik = %f\n", log_lik);
        
            *_log_lik = log_lik;
        //    if (log_lik < *_log_lik + *eps)
        //        break; // time to compute the filtered probs and exit
        //    else
        //        *_log_lik = log_lik;

            /// The amount of observations in state i.
            /// State i-> j transition. 
            vector<int>           stateCounts(k,0);
            vector< vector<int> > stateTransMatrix(k, vector<int>(k,1)); // assume minimal number of accurrance as 1.

            for(int m = 0; m < T; m++)
            {
                stateCounts[hidden_states[m]] ++;
                if(m < T -1)
                {
                    stateTransMatrix[hidden_states[m]][hidden_states[m +1]] ++;
                }
            }


            /// Re-estimate the parameters
            /// 1. pi
            for(int i = 0; i < k; i++)
            {
                /// preserving
                old_pi[i] = pi[i];
                if(stateCounts[i] != 0)
                    pi[i] = log(stateCounts[i]/double(T + k*k));
            
            }
            /// 2. tpm (a in Rabiner)
            for(int i = 0; i < k; i++)
            {
                for(int j = 0; j < k; j++)
                {
                    /// preserving
                    delete(p_tpm);
                    p_tpm = new vector<vector<double> >(tpm);
                    if(stateCounts[i] != 0) 
                        tpm[i][j] = log( stateTransMatrix[i][j] / double(stateCounts[i]) );
                }
            }

            /// 3. Get tumor cellularity using MLE
            //
            double cc_init  = 0.5;
            double cc_inter = -1;

            int* geno = (int*)(malloc(sizeof(int)*T*2)); //  {m,f} x T         
            for (unsigned each = 0; each < T; each ++)
            {
                if (hidden_states[each] < 0 ||  hidden_states[each] >= k) 
                    Rprintf("wrong state = %d", hidden_states[each]);
                geno[each * 2 ]   = genotypes[(hidden_states[each]) * 2];
                geno[each * 2 +1] = genotypes[(hidden_states[each]) * 2 + 1];
            }

            maxlik(depths, geno, &T, &cc_init, &cc_inter, &verbal);

            free(geno);
            Rprintf("Predicted tumor cellularity cc = %f. \n", cc_inter);
            if(*fix_tc == false)
            {
                p_tc = cellularity_iter;
                cellularity_iter  = cc_inter;
                Rprintf("Changing tc to %e \n", cellularity_iter);
            }
            else
            {

                Rprintf("Keeping tc as %e\n", cellularity_iter);
            }
            Rprintf("--------------------------------------------- \n");
            iter++;
        }
        *maxiter = iter + 1;
        for(int i = 0; i < k; i++)
            for(int j = 0, offs = 0; j < k; j++, offs += k)
                TPM[offs + i] = exp((*p_tpm)[i][j]);

        // Decode the hidden states using Viterbi alg.
        for(int i =0; i < k; i++)
            pi[i] = old_pi[i];
        *tc = p_tc;

            Rprintf("Predicted tumor cellularity cc = %f. when ll = %e \n", p_tc, *_log_lik);
        viterbi(*p_yll, *p_tpm, pi, hidden_states);
        delete(p_tpm);
        delete(p_yll);


         Rprintf("------------- TPM ---------------- \n");
         for(int i = 0; i < k; i++)
         {
            for(int j = 0, offs = 0; j < k; j++, offs += k)
            {
                Rprintf("%f \t", TPM[offs + i]);
            }
                Rprintf("\n");
         }
         Rprintf("------------- finished ---------------- \n");
    }
}

extern "C"
{
//        genotype(k * 2, 0); //  <maternal alleles,  paternal alleles> for kth state(genotype).
//        c                   // c for tumor cellularity.
//
    void
    Baum_Welch(int *seq_len, int *_k, int* depths,
                             int* genotypes, double* tc, double *TPM, double *pi, int *maxiter,
                             double *eps, double *_log_lik, double *filter,
                             int *hidden_states, double* prob, bool* fix_tc,
                             double* DOA, double* DOA_range, bool *print_info, double* globe_rcov)

    {
        printf("[notice me] g_rcov  = %f\n", *globe_rcov);
        // Initialization
        bool    verbal = true;
        int k = *_k, T = *seq_len;

        vector< vector<double> > tpm(k, vector<double>(k)),
            yll(k, vector<double>(T)), alpha(k, vector<double>(T)),
            beta(k, vector<double>(T));
        vector<vector<long double> > gamma(k, vector<long double>(T)),
            ksi(k, vector<long double>(k));
        vector<double> old_pi(k), old_means(k);

        for(int i = 0; i < k; i++)
            for(int j = 0, offs = 0; j < k; j++, offs += k)
                tpm[i][j] = log(TPM[offs + i]);

        *_log_lik = R_NegInf;
        int    iter = 0;
        double tc_range[2] = {0.00001, 0.99999};
        double del_range[2] = {0.00001, 0.99999};
        double cellularity_iter = (tc_range[0] + tc_range[1]) / 2;

//cellularity_iter = 0.90;

        *DOA = (DOA_range[0] + DOA_range[1]) / 2;
        double deletion = del_range[0] ;

        double g_rcov[2] = {*globe_rcov - 0.1, *globe_rcov};

        while(iter < *maxiter)
        {
            
             /// 1. Get tumor cellularity using MLE

            double cc_iter = -1;
            //double rr_range[2] = {0.5, 2.5};
            /// converting tpm to 2d array
            Rprintf("[info] Optimizing ... \n");
            double **i_tpm = AllocateDynamicArray<double>(k,k);
            for (int i=0; i < k; i++)
                for (int j=0; j < k; j++)
                    i_tpm[i][j] = tpm[i][j];

            if(iter  == 0)
            {
                // global search
                getTC_Ratio_2(filter, depths, genotypes, i_tpm, pi, &T, &k, 
                    tc_range,   &cellularity_iter,
                    DOA_range,  DOA,
                    del_range,  &deletion,
                    g_rcov, 1);
            }
            else
            {
                //local search
                getTC_Ratio_2(filter, depths, genotypes, i_tpm, pi, &T, &k, 
                    tc_range,   &cellularity_iter,
                    DOA_range,  DOA,
                    del_range,  &deletion,
                    g_rcov, 0);

            }
    ///       getDel(filter, depths, genotypes, i_tpm, pi, &T, &k, 
    ///               cellularity_iter,
    ///               *DOA,
    ///               del_range,  &deletion,
    ///               g_rcov);



            FreeDynamicArray<double>(i_tpm);

            printf("degree of deletion is : %f \n", deletion);
            Rprintf("Predicted tumor cellularity ratio = %f. \n", *DOA);
            Rprintf("Predicted tumor cellularity cc = %f. \n", cellularity_iter);


            emissionDist(depths, T, genotypes, k, cellularity_iter, *DOA, yll, deletion);
            bool faster = false;
            Rprintf("******* Forward-backward alpha-beta.\n");
            double log_lik = forward_backward(yll, tpm, pi, false, filter,
                                              alpha, beta, faster, *print_info);
            Rprintf("@ log_lik = %f when tc = %f, ratio = %f, del = %f\n", 
                    log_lik, cellularity_iter, *DOA, deletion);
            // two successive run changes < 1
            if(fabs(log_lik - *_log_lik) < 1)
            {
                *_log_lik = log_lik;
                Rprintf("@ two successive run changes < 1 \n");
                break;
            }
            *_log_lik = log_lik;

           // Calculate gamma and ksi
            vector<long double> sum_gamma(k, 0);
            vector< vector<long double> > sum_ksi(k, vector<long double>(k, 0));
            for(int m = 0; m < T; m++)
            {
                //double avfb = 0, avksi = 0;
                long double t_avfb = 0, t_avksi = 0;

#pragma omp parallel for reduction(+:t_avksi,t_avfb)
                for(int i = 0; i < k; i++)
                {
                    t_avfb += pArithmetic::myExp(gamma[i][m] = alpha[i][m] + beta[i][m], "t_avfb");
                    if (m < T - 1)
                        for(int j = 0; j < k; j++)
                            t_avksi += pArithmetic::myExp(ksi[i][j] = alpha[i][m] + tpm[i][j] + yll[j][m + 1]
                                        + beta[j][m + 1], "t_avksi");
                }
                if(boost::math::isinf(1/t_avksi) && m < T - 1)
                {
                    logging("[info] t_avksi == 0 when updating TPM.","error.txt");
                    exit(0);
                }
#pragma omp parallel for
                for(int i = 0; i < k; i++)
                {
                    long double t_gamma_i_m = pArithmetic::myExp(gamma[i][m] , "gamma[i][m]") / t_avfb;
                    //gamma[i][m] = pArithmetic::myExp(gamma[i][m] , "gamma[i][m]") / t_avfb;
                    if (m < T - 1)
                    {
#pragma omp atomic
                        sum_gamma[i] +=  t_gamma_i_m;

                        for(int j = 0; j < k; j++)
                        {
                            long double t_sum_ksi_i_j = pArithmetic::myExp(ksi[i][j] , "sum_ksi[i][j]") / t_avksi;
#pragma omp atomic
                            sum_ksi[i][j] += t_sum_ksi_i_j;
                            //sum_ksi[i][j] += pArithmetic::myExp(ksi[i][j] , "sum_ksi[i][j]") / t_avksi;
                        }
                    }
                }
            }


            /// 2. Updating pi
            /// stop updating when pi[i] is tiny ...
            Rprintf("[INFO] Updating Pi ...\n" );
            bool stopUpdate = false;
            for(int i = 0; i < k; i++)
            {
                old_pi[i] = pi[i];
                // pi[i] = gamma[i][0];
                pi[i] = pArithmetic::myLog(gamma[i][0], "pi");
                if(pi[i] < -1000) stopUpdate =true;
            }
            if(stopUpdate)
                for(int i = 0; i < k; i++)
                     pi[i] = old_pi[i];
            /// 3. Updating tpm 
            Rprintf("[INFO] Updating TPM ...\n" );
            for(int i = 0; i < k; i++)
                for(int j = 0; j < k; j++)
                    tpm[i][j] = pArithmetic::myLog(sum_ksi[i][j] / sum_gamma[i], "tpm[i][j]");
            Rprintf("--------------------------------------------- \n");
            iter++;
        } // EM iteration

        for(int t = 0; t < T; t++)
        {
            for(int i = 0; i < k; i++)
            {
                    prob[t * k + i ] = gamma[i][t];
            }
        }
        *tc = cellularity_iter;
        // Decode the hidden states using Viterbi alg.
        Rprintf("Predicted tumor cellularity cc = %f. when ll = %e \n", *tc, *_log_lik);
        viterbi(yll, tpm, pi, hidden_states);

        Rprintf("------------- finished ---------------- \n");
        
    }

}



