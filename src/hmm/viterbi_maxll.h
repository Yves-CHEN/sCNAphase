/************************************************************************
 *  Author: Wenhan CHEN
 *  Date  : 26 Aug 2014
 *  Last_Modified : 12 Sep 2015 18:29:02
 *  Description: This is for estimation of tumor cellularity given the 
 *    genotypes at every heterozygous loci using maximum likelihood 
 *    funciton. NLOPT package is used for the purpose of optimation.
 *
 * 
************************************************************************/


# include "nlopt.h"
# include <stdio.h>
# include <R.h>
# include <Rmath.h>
# include <vector>
//# include "global_vars.h"

using namespace std;


/* -------------------------------------------------------
*
*    DEFINE TYPEs.
*
 ------------------------------------------------------- */
typedef struct {
        double a, b;
} my_constraint_data;


typedef struct {
        unsigned m, p;
} Genotype;



typedef struct {
//        | Maternal      |
// ------------------------------------
//        | Alleic Depth  | Read Depth
// normal | $m_i^n$       | $d_i^n$ 
// tumor  | $m_i^t$       | $d_i^t$
    
        unsigned m_n, d_n,
                 m_t, d_t;
} DepthTable;

typedef struct {
    Genotype*   genoSet;
    DepthTable* dTableSet;
    int         T;
} MChain ;


//  typedef struct {
//      int         count  = 0;
//      bool        verbal = false;
//  } IOControl ;
//
//
bool ML_VERBAL = false;

/* -------------------------------------------------------
*
*    DECLARE FUNCTIONs.
*
 ------------------------------------------------------- */


// a) objective function definition
//double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
// b) definition constraints
double myconstraint(unsigned n, const double *x, double *grad, void *data);

// Main logic of optimization
double maxll(MChain* obersvations, double cc, int T);


// Wrapper  for maxll
void maxlik(int* depth, int* genotypes, int* numOfObser, double* cc, double* res);


/* -------------------------------------------------------
*
*    DECLARE FUNCTIONs.
*
 ------------------------------------------------------- */
double viterbi_probOnly(vector< vector<double> > &yll, vector< vector<double> > &tpm,
             double *pi)
{

    int T = yll[0].size(), k = yll.size();
    vector< vector<double> > log_lik_table(k, vector<double>(T));
    vector< vector<unsigned int> > ptrs(k, vector<unsigned int>(T));
    // Initialization
    vector<double> v(k);

    for(int i = 0; i < k; i++)
        log_lik_table[i][0] = yll[i][0] + pi[i];
    for(int m = 1; m < T; m++)
        for(int i = 0; i < k; i++)
        {

            int best_ind = 0;
            double s = R_NegInf;
            
            for(int j = 0; j < k; j++)
            {
                
                v[j] = log_lik_table[j][m - 1] + tpm[j][i];
                if (v[j] > s)
                    s = v[j], best_ind = j;

            }
            log_lik_table[i][m] = yll[i][m] + v[ptrs[i][m] = best_ind];

        }

    //	Track back the Viterbi path using the pointers stored in the
    //	previous step

    double s = log_lik_table[0][T - 1];

    for(int i = 1; i < k; i++)
        if (log_lik_table[i][T - 1] > s)
            s = log_lik_table[i][T - 1];

    return s;

}


double viterbi(vector< vector<double> > &yll, vector< vector<double> > &tpm,
             double *pi, int *states)
{

    int T = yll[0].size(), k = yll.size();
    vector< vector<double> > log_lik_table(k, vector<double>(T));
    vector< vector<unsigned int> > ptrs(k, vector<unsigned int>(T));
    // Initialization
    vector<double> v(k);

    for(int i = 0; i < k; i++)
        log_lik_table[i][0] = yll[i][0] + pi[i];
    for(int m = 1; m < T; m++)
        for(int i = 0; i < k; i++)
        {

            int best_ind = 0;
            double s = R_NegInf;
            
            for(int j = 0; j < k; j++)
            {
                
                v[j] = log_lik_table[j][m - 1] + tpm[j][i];
                if (v[j] > s)
                    s = v[j], best_ind = j;

            }
            log_lik_table[i][m] = yll[i][m] + v[ptrs[i][m] = best_ind];

        }

    //	Track back the Viterbi path using the pointers stored in the
    //	previous step

    double s = log_lik_table[0][T - 1];

    states[T - 1] = 0;
    for(int i = 1; i < k; i++)
        if (log_lik_table[i][T - 1] > s)
            s = log_lik_table[i][T - 1], states[T - 1] = i;
    // 	Rprintf("%d\n", hidden_states[T - 1]);
    for(int m = T - 2; m >= 0; m--)
        states[m] = ptrs[states[m + 1]][m + 1];

    return s;
    //return log_lik_table[states[T - 1]][T - 1];

}


   
double myfunc(unsigned n, const double *cc, double *grad, void *obser)
{
    MChain* obs = (MChain*) (obser);
    Genotype*    genoSet   = obs -> genoSet;
    DepthTable*  dTableSet = obs -> dTableSet;
    int    const T = obs -> T;
    double const C = *cc;
    double ll   = 0;
    double gradient = 0;

    if(ML_VERBAL)
        Rprintf("\nCellularity : %f . \n", *cc);
    for(unsigned i = 0; i < T; i ++)
    {
        int     x = genoSet[i].m;
        int     y = genoSet[i].p;
        int   m_i = dTableSet[i].m_n;
        int   d_i = dTableSet[i].d_n;
        int   m_i_t = dTableSet[i].m_t;
        int   d_i_t = dTableSet[i].d_t;


        //Rprintf("d_i = %d, m_i = %d, x=%d, y = %d\n", d_i,m_i,x,y);
        double Di = (1 - C) * d_i + C * m_i * x + C * (d_i - m_i) *y;
        double N_i       = ((1-C) * m_i + C * m_i * x) / Di;
        double N_i_prime = m_i * (x-1) / Di + 
            (m_i - C * m_i + C * m_i * x) * (d_i - m_i*x -d_i*y + m_i*y) / (Di*Di);
        ll        += m_i_t * log(N_i)        + (d_i_t - m_i_t) * log(1-N_i);
        gradient  += m_i_t / N_i * N_i_prime - (d_i_t - m_i_t) / (1-N_i) * N_i_prime;
        //Rprintf("Di = %f, Ni= %f, Ni_p=%f, ll = %f, grad = %f\n", Di,N_i,N_i_prime, ll,gradient);

    }

    if (grad) {
        grad[0] = gradient;
    }
    if(ML_VERBAL == true)
        Rprintf("ll = %10.8g, gradiet = %10.8g\n", ll, gradient);
    return ll;
}


// constraint for C:  C>0 && C <1
// double conC(unsigned n, const double *cc, double *grad, void *data)
// {
//     if (grad) {
//     grad[0] = 1;
//     }
// 
//     cc[0] -1
//     
//     my_constraint_data *d = (my_constraint_data *) data;
//     double a = d->a, b = d->b;
//     if (grad) {
//         grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
//         grad[1] = -1.0;
//     }
//     return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
// }
// 
// constaint for N_i:  N_i > 0 && N_i < 1
//  double conNi(unsigned n, const double *x, double *grad, void *data)
//  {
//          
//      my_constraint_data *d = (my_constraint_data *) data;
//      double a = d->a, b = d->b;
//      if (grad) {
//      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
//      grad[1] = -1.0;
//      }
//      return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
//  }

//
//  depth:
// A tuple of four values for each of the loci is in the order of:
//   
//           |   AD  |   RD
//   normal  |  mi_n |  di_n
//   tumor   |  mi_t |  di_t
//
//   genotypes:
//     A tuple of <maternal AD, paternal AD>

void maxlik(int* depth, int* genotypes, int* numOfObser, double* cc, double* res, bool* verbal)
{
    int T = *numOfObser;
    Genotype*   genoSet   = (Genotype*)   (malloc(sizeof(Genotype) * T));
    DepthTable* dTableSet = (DepthTable*) (malloc(sizeof(DepthTable) * T));

    MChain observ = {genoSet, dTableSet, T};

    //ML_VERBAL = true;
    ML_VERBAL = *verbal;


    //
    // Initialization of the observations
    for(unsigned i = 0; i < T; i ++)
    {
        dTableSet[i].m_n = depth[i*4 + 0];
        dTableSet[i].d_n = depth[i*4 + 1];
        dTableSet[i].m_t = depth[i*4 + 2];
        dTableSet[i].d_t = depth[i*4 + 3];
        genoSet[i].m     = genotypes[i*2 + 0];
        genoSet[i].p     = genotypes[i*2 + 1];

    //      Rprintf(" dTableSet[i].m_n   = %d \n",  depth[i*4 + 0]);
    //      Rprintf(" dTableSet[i].d_n   = %d \n",  depth[i*4 + 1]);
    //      Rprintf(" dTableSet[i].m_t   = %d \n",  depth[i*4 + 2]);
    //      Rprintf(" dTableSet[i].d_t   = %d \n",  depth[i*4 + 3]);
    //      Rprintf(" genoSet[i].m       = %d \n",  genotypes[i*2 + 0]);
    //      Rprintf(" genoSet[i].p       = %d \n",  genotypes[i*2 + 1]);
    }

    *res = maxll(&observ, *cc, T);
    free(genoSet);
    free(dTableSet);
}

double maxll(MChain* observations, double ccInit, int T)
{


    double lb[1] = {0.0000001}; // lower bounds.
    double ub[1] = {0.999999};  // upper bounds.
    double paraMiniChange[1] = 
                   {1e-8};   // optimization parameter minimal changes each step.
    double cc = ccInit;         // initial guess.
    double minf;                // the minimum objective value, upon return.
    int    error = 1;           // no errors by default.

    nlopt_opt opt;
    //opt = nlopt_create(NLOPT_LD_MMA, 1); /* algorithm and dimensionality */
    //opt = nlopt_create(NLOPT_LD_LBFGS, 1); /* algorithm and dimensionality */
    opt = nlopt_create(NLOPT_LN_COBYLA, 1); /* algorithm and dimensionality */


    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_max_objective(opt, myfunc, observations);
    //   my_constraint_data data[2] = { {2,0}, {-1,1} };
    //   nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
    //   nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
    nlopt_set_xtol_rel(opt, 1e-12);
    nlopt_set_xtol_abs(opt, paraMiniChange);

    int errCode = nlopt_optimize(opt, &cc, &minf);
    if ( errCode < 0) {
        Rprintf("nlopt failed!  errorCode == %d\n", errCode);
        Rprintf(" NLOPT_ROUNDOFF_LIMITED = -4\n");
    }
    else {
        if(ML_VERBAL)
            Rprintf("found maximum at f(%g) = %0.10g\n", cc, minf);
    }


    nlopt_destroy(opt);
    return cc;
}



