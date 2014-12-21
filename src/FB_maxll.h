/************************************************************************
 *  Author: Wenhan CHEN
 *  Date  : 22 Sep 2014
 *  Last_Modified : 04 Dec 2014 14:53:49
 *  Description: This includes 
 *    1) Implementation of forward-backward algorithm.
 *    2) NLopt functions for calculation of tc against likelihood of the
 *    sequential obsersions from forward-backward algorithm.
 *
 *    NLopt package is MIT implementation of a set of optimization algorithms.
 *    Can be found at:
 *      http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#BOBYQA
 * 
************************************************************************/

#include <R.h>
#include <Rmath.h>
#include <vector>
#include "utils.h"
#include "optimizer.h"

using namespace std;

extern double logProb(unsigned int mi_n , unsigned int  di_n , unsigned int  mi_t , unsigned int  di_t,
        double c,vector<int> genotype);
extern double forward_backward (vector< vector<double> > &yll, 
                              vector< vector<double> > &tpm,
                              double *pi, bool do_filter, double *filter,
                              vector< vector<double> > &alpha,
                              vector< vector<double> > &beta,
                              bool faster,
                              bool print_info);


class  TFunc_test : public ObjectiveFunc
{

public:
    /// vector<double> initParameters; // inherited from ObjectiveFunc
    int*      genotypes;
    int*      depth;
    double*   pi;
    double*   filter;
    int       T;
    int       k;
    bool      doFilter;

    vector< vector<double> >   tpm  ;
    vector< vector<double> >   yll  ;
    vector< vector<double> >   alpha;
    vector< vector<double> >   beta ;
    inline TFunc_test (int* depth, int* genotypes, double* pi, double** tpm, double* filter,
            int* numOfObser, int* numOfStates, vector<double> initPara)
    :depth(depth), genotypes(genotypes), pi(pi), T(*numOfObser), k(*numOfStates), filter(filter), ObjectiveFunc(initPara)
    {
        doFilter = false;
        this->tpm.resize   (k, vector<double>(k, 0));
        for(int i = 0; i < k; i++)
        {
            for(int j = 0; j < k; j++)
            {
                this->tpm[i][j] = tpm[i][j] ;
            }
        }

        yll.resize   (k, vector<double>(T, 0));
        alpha.resize (k, vector<double>(T, 0));
        beta.resize  (k, vector<double>(T, 0));
    };

    inline double objective(const double* paras)
    {
        double tc  = paras[0];
        /// Initialization of the observations
        for(int m = 0; m < T; m++)
        {
            int   m_i   = depth[m * 4 + 0];
            int   d_i   = depth[m * 4 + 1];
            int   m_i_t = depth[m * 4 + 2];
            int   d_i_t = depth[m * 4 + 3];
            for(int i = 0; i < k; i++)
            {
//int cnGain = genotypes[i];
                vector<int> g_i(2,0);
                g_i[0] = genotypes[i * 2];
                g_i[1] = genotypes[i * 2 +1];

                yll[i][m] = logProb(m_i, d_i, m_i_t, d_i_t, tc, g_i);
            }
        }
        bool faster = false;
        bool print_info = false;
        double log_lik = forward_backward(yll, tpm, pi, doFilter, filter, alpha, beta, faster,print_info);
        return  log_lik;
    };

    inline double applyConstraint(const double* paras)    
    {
        return 0;
    }
    vector<double> getSolution()
    {
        return ObjectiveFunc::initParameters;
    };
};



class MaxLikelihood_test : public Optimizer
{

public :
    /// Followings are parameters inherited from super.
    //ObjectiveFunc*  targetF;
    inline MaxLikelihood_test(ObjectiveFunc* targetF)
        :Optimizer(targetF)
    {
  //targetF = new TFunc(depth, genotypes, pi, tpm, filter, numOfStates, numOfStates);
    };

    inline void doOptimize()
    {
        double lb[1] = {0.00001};             // lower bounds.
        double ub[1] = {0.99999};              // upper bounds.
        double paraMiniChange[1] = {1e-8};      // optimization parameter minimal changes each step.
        double maxTimes  = 500;
        opt = nlopt_create(NLOPT_GN_DIRECT_L, 1);   // Choose algorithm
        nlopt_set_lower_bounds(opt, lb);            // set stop criteria
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_xtol_rel(opt, 1e-12);             // this would not apply to Globle search.
        nlopt_set_xtol_abs(opt, paraMiniChange);    // as above
        nlopt_set_maxeval(opt, maxTimes);           // evaluate optimazation for < 100;
        Optimizer::doOptimize();

    }
};

void getTC(double* filter, int* depth, int* genotypes, double** tpm, double* pi,
        int* numOfObser, int* numOfStates, double* cc, double* res_tc)

{
    vector<double> paras ;
    paras.push_back(*cc);
    TFunc_test ff (depth, genotypes, pi, tpm, filter, numOfStates, numOfStates, paras);
    MaxLikelihood_test mm (&ff);
    mm.doOptimize();
    paras = ff.initParameters;
    *res_tc = paras[0];
}





/* --------------------------------------------
 *           Main Logics
 *-------------------------------------------*/

double forward_backward (vector< vector<double> > &yll, 
                              vector< vector<double> > &tpm,
                              double *pi, bool do_filter, double *filter,
                              vector< vector<double> > &alpha,
                              vector< vector<double> > &beta,
                              bool faster,
                              bool print_info)
{
    // Initialization
    int T = yll[0].size(), k = yll.size();


    double sum_log_lik = 0;
    long double ll = 0;
    for(int i = 0; i < k; i++)
    {
        ll += pArithmetic::myExp(alpha[i][0] = yll[i][0] + pi[i], string("ll"));
        beta[i][T - 1] = 0;
    }
    sum_log_lik =  pArithmetic::myLog(ll, "sum_log_lik");
    for(int i = 0; i < k; i++)
    {
        alpha[i][0] -= sum_log_lik;
    }
    // Calculate the observed log-likelihood using scaling variables and
    // a dynamic programming table.
    //
    //


    for(int m = 1, t = T - 2, offs = k; m < T; m++, t--)
    {
        long double avf = 0, avb = 0;
        for(int i = 0; i < k; i++)
        {
            long double t_alpha_i_m =0,   t_beta_i_t  =0;
            for(int j = 0; j < k; j++)
            {
                
                t_alpha_i_m += pArithmetic::myExp(yll[i][m]     + alpha[j][m - 1] + tpm[j][i], "t_alpha_i_m");
                if(!faster)
                    t_beta_i_t  += pArithmetic::myExp(yll[j][t + 1] +  beta[j][t + 1] + tpm[i][j], "t_beta_i_t");
                if (isnan(t_alpha_i_m) || isnan(t_beta_i_t))
                {
                    Rprintf("t_alpha_i_m == nan");
                    exit(0);
                }
            }
//                if(t < 2)
//                   {
//                       Rprintf("beta[%d][%d] = %Le\n", t, i, t_beta_i_t);
//                   }

            avf += t_alpha_i_m;
            if(!faster)
                avb += t_beta_i_t;
            alpha[i][m] = pArithmetic::myLog(t_alpha_i_m, "alpha[i][m]");
            if(!faster)
                beta[i][t]  = pArithmetic::myLog(t_beta_i_t, "beta[i][t]");
        }

        avf = pArithmetic::myLog(avf, "avf"), avb = pArithmetic::myLog(avb, "avb"), ll = 0;
        for(int i = 0; i < k; i++)
        {
/// if alpha is a tiny value. e.g. -1.1e100, log(exp(alpha)) will not give a -inf then -1.1e100.
            alpha[i][m] = alpha[i][m] - avf;
            if(!faster)
                beta[i][t]  = beta[i][t] - avb;

        }

        double b_sum = sum_log_lik;
        sum_log_lik += double(avf);
        if(sum_log_lik - b_sum > 0)
        {
            Rprintf("@ overflowing b_sum = %e\t sumloglik = %e\n", double(b_sum), double(sum_log_lik));
        }

    }
    return sum_log_lik;

}



