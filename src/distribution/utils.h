/************************************************************************
 *  Author: Wenhan CHEN
 *  Date  : 22 Sep 2014
 *  Last_Modified : 19 Jun 2015 00:14:11
 *  Description: This includes 
 *    1) Boost implementation of Binomial_distribution.
 *    2) My implementation of log and exp at the long double precision.
 *    3) Dynamic 2D array creation and destruction.
 *
 * 
************************************************************************/
#ifndef PARITH
#define PARITH 


#include <omp.h>
#include <string>
#include <fstream>
#include <R.h>
#include "boost/math/distributions/binomial.hpp"

#include <sys/times.h>

#define CSTACK_DEFNS 7  /// These are for modify the default stack usages when using OpenMP.
#include "Rinterface.h" ///    by setting the default to unlimited.  R_CStackLimit=(uintptr_t)-1;



using namespace std;
using boost::math::binomial;

/* --------------------------------------------
 *           Class utilities
 --------------------------------------------*/


extern const bool printWarning;


template<typename T>
int logging(T con, const char* filename)
{
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    outfile << con;
    return 0;
}



/// Handle precision of log and exp operations when the value is tiny.
class pArithmetic
{
    public :
    static const long double logInf_ld ;
    static const long double Inf_ld    ;

    static const double logInf ;
    static const double Inf    ;

    ///  This gives a trancation warning message if the value is too small for expl(logInf).
    ///  The return value is double, because the smallest value in exp space is e-2500, which 
    ///  can be represented in double in log space.
    inline static double myLog(long double val, string who)
    {
        if(val < Inf_ld)
        {
            if(printWarning)
                Rprintf("[warning] trancating value when log(%Le) performing %s.\n", val, who.c_str());
            return logInf_ld;
        }
        return logl(val);
    }
    inline static long double myExp(double val, string who)
    {
        if(val < logInf_ld)
        {
            if(printWarning)
                Rprintf("[warning] exp(%e) generates a 0 performing %s.\n", val, who.c_str());
            return Inf_ld;
        }
        return expl(val);
    }
};
const long double pArithmetic::logInf_ld = logl(1.000e-4500L) ;
const long double pArithmetic::Inf_ld    = expl(logl(1.000e-4500L)) ;
const double pArithmetic::logInf = log(1e-310) ;
const double pArithmetic::Inf    = exp(log(1e-310)) ;


/* --------------------------------------------
 *           Globel functions
 --------------------------------------------*/
// double log_beta(long double a, long double b)
// {
//     return (a+b-2) * log(0.5) - log(boost::math::ibeta_derivative(a,b,0.5));
// }

double log_beta(long double a, long double b)
{
    double p0 = a / (a+b);
    return (a -1)* log(p0) + (b-1) * log(1-p0) - log(boost::math::ibeta_derivative(a,b,p0));
}


double log_beta_binom(double k, double n, double a, double b)
{
    double res = log_beta(a + k ,  n-k + b) - log_beta(a, b) - log_beta(k+1, n-k+1) - log (n+1);
    return res;
}



double logBinomial_boost(unsigned int mi_t , unsigned int di_t, double prob)
{
    binomial model_1(di_t, prob);
    double emp=pdf(model_1, mi_t);
    if(emp < pArithmetic::Inf)
        return pArithmetic::logInf;
    else
        return log(emp);
}
double logBinomial(unsigned int mi_t , unsigned int di_t, double prob)
{   
    return  (mi_t * log(prob) + (di_t - mi_t) * log (1-prob) );
}

/// c for tumor cellularity.
/// genotype is a vector of 2 values for the number of maternal alleles and paternal alleles.
double logProb(unsigned int mi_n , unsigned int  di_n , unsigned int  mi_t , unsigned int  di_t,
                double c, int genotype[])
{
    unsigned int x = genotype[0];
    unsigned int y = genotype[1];
    double mix_m = (1-c) * mi_n + c * mi_n * x;  // mixed maternal alleles
    double sum   = (1-c) * di_n + c * mi_n * x + c * (di_n - mi_n) * y;   // all alleles
    if(sum == 0 )
    {
        // not a fitted a binomial, this happens when genotype is {}, and tumor
        // cellularity is close to 100%. The  number of reads (di_t) are
        // misaligned reads in this case. Assume the misalign rate is 1/1000.
        return di_n * log(1/1000);
    }

    //if(sum == 0) {printf("[error] sum == 0 in logProb \n"); exit(-1);}
    double prob  = mix_m/sum;
    if(prob == 0) {Rprintf("[warning] prob == 0 in logProb. Assigned to 1e-10 \n"); prob = 1e-10;}
    if(prob == 1) {Rprintf("[warning] prob == 0 in logProb. Assigned to 1e-10 \n"); prob = 1- 1e-10;}
    /// using beta-binomial dist
    //return log_beta_binom (mi_t ,  di_t, mix_m, sum - mix_m);


    /// using binomial dist
    return logBinomial_boost (mi_t ,  di_t, prob);
    //return logBinomial(mi_t ,  di_t, prob);
}
double logProb_depth(unsigned int di_n , double  sum_di_n , unsigned int  di_t , double  sum_di_t,
                double c, double ratio, double cnGain)
{
    double mix_m = (1-c) * di_n + c * di_n *( 0.5 * cnGain  + 1);  // mixed maternal alleles
    double sum   =  sum_di_n * (ratio * c + 1 - c) ;   // all alleles
    if(sum == 0)  {Rprintf("[error] sum == 0 in logProb \n"); exit(-1);}
    double prob  = mix_m/sum;
    if(prob == 0) {Rprintf("[warning] prob == 0 in logProb. Assigned to 1e-10 \n"); prob = 1e-10;}
    /// using beta-binomial dist
    //return log_beta_binom (di_t ,  sum_di_t, mix_m, sum - mix_m);
    //
    /// using binomial dist
    return logBinomial_boost (di_t ,  sum_di_t, prob);
    //return logBinomial(di_t ,  sum_di_t, prob);
    //return log_beta_binom(di_t ,  sum_di_t, mix_m, sum - mix_m);
}

double logProb_depth_del(unsigned int di_n , double  sum_di_n , unsigned int  di_t , double  sum_di_t,
                double c, double ratio, double cnGain, double deletion)
{
    double mix_m = (1-c) * di_n + c * di_n *( 0.5 * cnGain  + 1);  // mixed maternal alleles

    //printf("here deletion : %f \n", deletion);

    double sum   = (1-deletion) * sum_di_n * (ratio * c + 1 - c) ;   // all alleles
    if(sum == 0)  {Rprintf("[error] sum == 0 in logProb \n"); exit(-1);}
    double prob  = mix_m/sum;

    if(prob == 0) {Rprintf("[warning] prob == 0 in logProb. Assigned to 1e-10 \n"); prob = 1e-10;}
    /// using beta-binomial dist
    //return log_beta_binom (di_t ,  sum_di_t, mix_m, sum - mix_m);
    //
    /// using binomial dist
    if(prob >= 1)
    {
        return  R_NegInf;
    }
    return logBinomial_boost (di_t ,  sum_di_t, prob);

    //return logBinomial(di_t ,  sum_di_t, prob);
    //return log_beta_binom(di_t ,  sum_di_t, mix_m, sum - mix_m);
}


//       g_rcov  = \frac{tumor}{( 1 - del) *  (normal * (tc*DOA + 1-tc))} when 0<= del  <=1
double get_rcov(int* depths, int T, double tc, double DOA,  double del)
{
    double sum_di_normal = 0;
    double sum_di_tumor  = 0;
    for(int m = 0; m < T; m++)
    {
        sum_di_normal += depths[m*4 +1];
        sum_di_tumor  += depths[m*4 +3];
    }
    double g_rcov = sum_di_tumor / ( (1 - del) * sum_di_normal *  (tc * DOA + 1 - tc));
    return g_rcov;
 
}

double get_del(int* depths, int T, double tc, double DOA,  double rcov)
{
    double sum_di_normal = 0;
    double sum_di_tumor  = 0;
    for(int m = 0; m < T; m++)
    {
        sum_di_normal += depths[m*4 +1];
        sum_di_tumor  += depths[m*4 +3];
    }
    double remain = sum_di_tumor / ( rcov * sum_di_normal *  (tc * DOA + 1 - tc));
    return 1 - remain;
 
}




void emissionDist(int* depths, int T, int* genotypes, int k, double tc, double DOA, vector< vector<double> >&   yll, double deletion)
{

    double sum_di_normal = 0;
    double sum_di_tumor  = 0;

    R_CStackLimit=(uintptr_t)-1;

    for(int m = 0; m < T; m++)
    {
        sum_di_normal += depths[m*4 +1];
        sum_di_tumor  += depths[m*4 +3];
    }

#pragma omp parallel
    {

#pragma omp for
        for(int m = 0; m < T; m++)
        {
            int largerIdx = 0;

            for(int i = 0; i < k; i++)
            {
                int g_i[2] = {0, 0};
                
                g_i[0] = genotypes[i * 2];
                g_i[1] = genotypes[i * 2 +1];
                yll[i][m] = 0;
                yll[i][m] = logProb(depths[m*4], depths[m*4 +1], depths[m*4 +2], depths[m*4 +3],
                                tc, g_i);
#pragma omp critical 
                 {
                     if(yll[i][m] > yll[largerIdx][m])
                     {
                         largerIdx = i;
                     }
                 }

                int cnGain = g_i[0] + g_i[1] -2;
#pragma omp atomic
                yll[i][m] += logProb_depth_del( depths[m*4 +1], sum_di_normal, depths[m*4 +3], sum_di_tumor,
                                tc, DOA, cnGain, deletion);

            }
            if(yll[largerIdx][m] < pArithmetic::logInf_ld/2)
            {
                yll[largerIdx][m] = pArithmetic::logInf_ld / 3;
                logging(m, "tiny.Prob.unsolved.loci.txt");
                logging("\n", "tiny.Prob.unsolved.loci.txt");
            }
        }
    }
}


void emissionDist_back(int* depths, int T, int* genotypes, int k, double tc, double DOA, vector< vector<double> >&   yll)
{
    double sum_di_normal = 0;
    double sum_di_tumor  = 0;

   // for(int m = 0; m < T; m++)
   // {
   //     sum_di_normal += depths[m*4 +1];
   //     sum_di_tumor  += depths[m*4 +3];
   // }
    for(int m = 0; m < T; m++)
    {
        for(int i = 0; i < k; i++)
        {
            vector<int> g_i(2,0);
            g_i[0] = genotypes[i * 2];
            g_i[1] = genotypes[i * 2 +1];
            
            yll[i][m] = 0;
            //depths[m*4] =1;
            //depths[m*4 + 1] =2;
        //    yll[i][m] = logProb(depths[m*4], depths[m*4 +1], depths[m*4 +2], depths[m*4 +3], 
        //            tc, g_i);

           // int cnGain = g_i[0] + g_i[1] -2;
           // yll[i][m] += logProb_depth( depths[m*4 +1], sum_di_normal, depths[m*4 +3], sum_di_tumor, 
           //         tc, DOA, cnGain);

            int cnGain = genotypes[i];


       //     if(i < int(k/2 ))
       //     {

       //         yll[i][m] += -1000000;
       //     }
       //     else
       //     {

               yll[i][m] += logProb_depth( depths[m*4 ], depths[m*4 +1], depths[m*4 +2], depths[m*4 +3], 
                   tc, DOA, cnGain);
 
         ///     yll[i][m] += logProb_depth( depths[m*4 +1], sum_di_normal, depths[m*4 +3], sum_di_tumor, 
         ///         tc, DOA, cnGain);

            //}
        

        }
        //sum_di_n += depths[m*4 +1];
        //sum_di_t += depths[m*4 +3];
    }



}



// tab : nB, nSum, tB, tSum
extern "C"
{
void maxll_improv (int* tab, int* num_n, int* genotypes, int* num_k, double* cc, int* state, double* max_ll)
{
    int const minusInf = -10000000;
    double mll = minusInf;
    double p_mll = minusInf;
    double q_mll = minusInf;
    double tc = *cc;
    std::vector<double> best( *num_n, mll);
    double bestPossible =0;
    *state  = -1;
    for(int k=0; k < *num_k; k ++)
    {
        double sum = 0;
        double p_sum = 0;
        double q_sum = 0;
        double x = genotypes[k*2];
        double y = genotypes[k*2 + 1];

        for(int n =0; n < *num_n; n ++)
        {
           int nB = tab[n*4], nSum = tab[n*4 + 1], tB = tab[n*4 + 2], tSum = tab[n*4 + 3];
           double prob = 0.5;
           if( !(x == 0 && y == 0 && (1-tc) < 0.01) )
           {
               prob = (nB * x * tc + nB * (1-tc))
                        / ((nSum -nB) * y * tc + nB * x * tc + nSum * (1-tc));
           }

           if(prob ==0)
           {
                prob = 0.00001;
           }
           else if(prob == 1)
           {
               prob = 0.9999;
           }
           //sum = sum + tB * log(prob) + (tSum - tB) * log(1-prob);
           double res=logBinomial_boost (tB ,  tSum, prob);
           sum = sum + res;
           if(res > best[n]) best[n] = res;
           if(k == *num_k -1)
           {
               bestPossible += best[n];
           }

           //sum = sum + log(pbinom(each$tB, each$tSum, prob))

           if(n < ((*num_n)/2) )
           {
               p_sum += res;
           }
           else
           {
               q_sum += res;
           }
        }

        if(p_sum > p_mll)
        {
            p_mll = p_sum;
        }
        if(q_sum > q_mll)
        {
            q_mll = q_sum;
        }
        if(sum > mll)
        {
            mll = sum;
            *state = k;
        }

    }

    *max_ll = 2 * mll - bestPossible - (p_mll + q_mll);
    //*max_ll = mll;
}
}

extern "C"
{
void maxll (int* tab, int* num_n, int* genotypes, int* num_k, double* cc, int* state, double* max_ll)
{
    double mll = -10000000;
    double tc = *cc;
    std::vector<double> best( *num_n, mll);
    double bestPossible =0;
    *state  = -1;
    for(int k=0; k < *num_k; k ++)
    {
        double sum = 0;
        double x = genotypes[k*2];
        double y = genotypes[k*2 + 1];
        for(int n =0; n < *num_n; n ++)
        {
           int nB = tab[n*4], nSum = tab[n*4 + 1], tB = tab[n*4 + 2], tSum = tab[n*4 + 3];
           if( ((nSum -nB) * y + nB * x) == 0) printf("[error] %f\n", ((nSum -nB) * y + nB * x) );
           double prob= (nB * x * tc + nB * (1-tc))
                        / ((nSum -nB) * y * tc + nB * x * tc + nSum * (1-tc));

           if(prob ==0)
           {
                prob = 0.00001;
           }
           else if(prob == 1)
           {
               prob = 0.9999;
           }
           //sum = sum + tB * log(prob) + (tSum - tB) * log(1-prob);
           double res=logBinomial_boost (tB ,  tSum, prob);

           sum = sum + res;
           if(res > best[n]) best[n] = res;
           if(k == *num_k -1)
           {
               bestPossible += best[n];
           }
//printf("Pdf from binomial(%d,%d,%f) = %f\n", tB, tSum, prob, (res) );
           //sum = sum + log(pbinom(each$tB, each$tSum, prob))
        }
        
        if(sum > mll)
        {
            mll = sum;
            *state = k;
        }

    }

    *max_ll = mll - bestPossible;
    //*max_ll = mll;
}
}




// double logProb2(unsigned int di_n , unsigned int  sum_di_n , unsigned int  di_t , unsigned int  sum_di_t,
//                 double c, vector<int> genotype)
// {
//     unsigned int x = genotype[0];
//     unsigned int y = genotype[1];
//     double mix_m = (1-c) * di_n + c * di_n * x;  // mixed maternal alleles
//     double sum   = (1-c) * di_n + c * mi_n * x + c * (di_n - mi_n) * y;   // all alleles
//     if(sum == 0) {Rprintf("[error] sum == 0 in logProb \n"); exit(-1);}
//     double prob  = mix_m/sum;
//     if(prob == 0) {Rprintf("[warning] prob == 0 in logProb. Assigned to 1e-10 \n"); prob = 1e-10;}
//     /// using beta-binomial dist
//     return log_beta_binom (mi_t ,  di_t, mix_m, sum - mix_m);
// 
// 
//     /// using binomial dist
//     //return logBinomial_boost (mi_t ,  di_t, prob);
//     //return logBinomial(mi_t ,  di_t, prob);
// }



template <typename T> 
T **AllocateDynamicArray( int nRows, int nCols)
{
      T **dynamicArray;

      dynamicArray = new T*[nRows];
      for( int i = 0 ; i < nRows ; i++ )
      dynamicArray[i] = new T [nCols];

      return dynamicArray;
}
template <typename T>
void FreeDynamicArray(T** dArray)
{
      delete [] *dArray;
      delete [] dArray;
}


#endif
