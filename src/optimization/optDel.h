#ifndef OPTDEL
#define OPTDEL


#include <R.h>
#include <Rmath.h>
#include "utils.h"
#include "optimizer.h"
using namespace std;

class  TFuncDel : public ObjectiveFunc
{

public:
    /// vector<double> initParameters; // inherited from ObjectiveFunc
    int*      genotypes;
    int*      depths;
    double*   pi;
    double*   filter;
    int       T;
    int       k;
    bool      doFilter;
    double*   g_rcov_range;
    double    tc;
    double    DOA;

    vector< vector<double> >   tpm  ;
    vector< vector<double> >   yll  ;
    vector< vector<double> >   alpha;
    vector< vector<double> >   beta ;


    inline TFuncDel (int* depths, int* genotypes, double* pi, double** tpm, double* filter,
            int* numOfObser, int* numOfStates, vector<double> initPara, double* g_rcov_range, double tc, double DOA)
    :depths(depths), genotypes(genotypes), pi(pi), T(*numOfObser), k(*numOfStates), filter(filter), ObjectiveFunc(initPara), g_rcov_range(g_rcov_range), tc(tc), DOA(DOA)
    {

        printf("g_rcov_range : %f \n", g_rcov_range[0]);
        doFilter = false;
        this->tpm.resize   (k, vector<double>(k, 0));
        for(int i = 0; i < k; i++)
            for(int j = 0; j < k; j++)
                this->tpm[i][j] = tpm[i][j] ;
        yll.resize   (k, vector<double>(T, 0));
        alpha.resize (k, vector<double>(T, 0));
        beta.resize  (k, vector<double>(T, 0));
    };

      
    inline double objective(const double* paras)
    {
        double del = paras[0];

        emissionDist(depths, T, genotypes, k, tc, DOA, yll, del);
        double est_rcov = get_rcov(depths, T, tc, DOA, del);

        bool faster = false;
        bool print_info = false;
        double log_lik = R_NegInf;
        if(est_rcov > g_rcov_range[0] && est_rcov < g_rcov_range[1])
            log_lik = forward_backward(yll, tpm, pi, doFilter, filter, alpha, beta, faster,print_info);

        Rprintf("\t [info] tc = %f\t DOA = %f when ll = %f, del = %f, est_rcov = %f \n",
                tc, DOA, log_lik, del, est_rcov);
        return  log_lik;
    };
    
    inline double applyConstraint(const double* paras)
    {
        double del = paras[0];
        double est_rcov = get_rcov(depths, T, tc, DOA, del);

        printf("invoked ... deletion = %f \t est_rcov = %f  \n", del, est_rcov);

        return g_rcov_range[0] - est_rcov;
        //return  (deletion - 0.999) && (0.0001 - deletion);
    };


    vector<double> getSolution()
    {
        return ObjectiveFunc::initParameters;
    };
};

class MaxLikelihoodDel : public Optimizer
{

public :
    /// Followings are parameters inherited from super.
    //ObjectiveFunc*  targetF;
    double lb[1];
    double ub[1];
    inline MaxLikelihoodDel(ObjectiveFunc* targetF, double lb[1], double ub[1])
        :Optimizer(targetF)
    {
        this->lb[0] = lb[0];
        this->ub[0] = ub[0];
    };

    inline void doOptimize()
    {
        double paraMiniChange[1] = {1e-8};      // optimization parameter minimal changes each step.
        double maxTimes  = 1000;
        //opt = nlopt_create(NLOPT_GN_DIRECT_L, 1);   // Choose algorithm
        //opt = nlopt_create(NLOPT_LN_BOBYQA, 1);   // 
        opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1);   // 
        //opt = nlopt_create(NLOPT_LN_COBYLA, 3); 
        nlopt_set_lower_bounds(opt, lb);            // set stop criteria
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_xtol_rel(opt, 1e-12);             // this would not apply to Globle search.
        nlopt_set_xtol_abs(opt, paraMiniChange);    // as above
        nlopt_set_maxeval(opt, maxTimes);           // evaluate optimazation for < 100;
        Optimizer::doOptimize();

    }
};


void getDel(double* filter, int* depth, int* genotypes, double** tpm, double* pi,
        int* numOfObser, int* numOfStates,
        double tc, double DOA,
        double* del_range, double* res_del,
        double* g_rcov_range)
{
    vector<double> paras ;
    /// set the mean value as the initial
    double lb[1] ;             // lower bounds.
    double ub[1] ;              // upper bounds.

    lb[0] = del_range[0];
    ub[0] = del_range[1];

    double e_lb = get_del(depth, *numOfObser, tc, DOA, g_rcov_range[0]);
    double e_ub = get_del(depth, *numOfObser, tc, DOA, g_rcov_range[1]);

    if(e_lb > del_range[0] && e_lb < del_range[1]) lb[0] = e_lb;
    if(e_ub < del_range[1] && e_ub > del_range[0]) ub[0] = e_ub;


    double initDel = (lb[0] + ub[0])/2;
    paras.push_back(initDel);
    TFuncDel ff (depth, genotypes, pi, tpm, filter, numOfObser, numOfStates, paras, g_rcov_range, tc, DOA);
    MaxLikelihoodDel mm (&ff, lb, ub);
    mm.doOptimize();
    paras = ff.initParameters;
    *res_del = paras[0];

   // double est_rcov = get_rcov(depth, *numOfObser, tc, DOA, *res_del);
   // Rprintf("\t @ [info] finish opt  tc = %f\t DOA = %f del = %f, est_rcov = %f\n",
   //             tc, DOA, *res_del, est_rcov);
}


#endif
