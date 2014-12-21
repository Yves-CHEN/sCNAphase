#ifndef Optimizer_BASE
#define Optimizer_BASE
#include <vector>
#include "nlopt.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;



/************************************************************************
Templates :

class SMaker : public Optimizer
{
    public:
    inline SMaker (ObjectiveFunc* objF):Optimizer(objF){}; // Don't modifiy this!
    void doOptimize()
    {
        /// set the optimazation parameters of nlopt from here
        /// nlopt_opt opt;


        
        Optimizer::doOptimize();
    }
};
==================================================================


class Tfunction : public ObjectiveFunc
{
    public:
    /// Define the parameters used in obj funciton.
    //  e.g.
    //    double factor;
    inline Tfunction(vector<double> initGuess, double factor)
        :ObjectiveFunc(initGuess), factor(factor)
    {
            // set the parameters.
    };
    inline double objective(const double* paras)
    {
        /// write the object funtions with the parameter to be estimated, and
        ///    all other known parameters set in member vars.
        /// e.g.
        //   double x = paras[0];
        //   return factor * pow(x, 2);
    };
};

 *
 *
 * ************************************************************************/

class Constraint
{
    public:
        virtual double applyConstraint(const double* paras) =  0;

       
};

/// 0. Inherit ObjectiveFunc. e.g.  class targetFunc : public ObjectiveFunc.
/// 1. Set the parameters for objective function as the class member.
/// 2. Set all the parameters in the constructor.
/// 3. Implement objective(vector<double> initParameters).
///
class ObjectiveFunc
{
    private:
    ObjectiveFunc();
    public:

    vector<double> initParameters;
    vector<void*>  constraints;

    inline ObjectiveFunc(vector<double> initParameters):initParameters(initParameters) {};
    virtual  double objective (const double* paras) =0 ;
    double applyConstraint(const double* paras) {return 0;} ;

};




extern "C"
{
    double targetFun_stub(unsigned n, const double* initPara, double *grad, void *op )
    {
        ObjectiveFunc* theOp  = (ObjectiveFunc*) (op);
        printf("[notice] pos initPara: %f \n", (theOp->initParameters[0]));
        printf("[notice] pos initPara: %f \n", (theOp->initParameters[1]));


        printf("[notice] initPara: %f \n", initPara[0]);
        printf("[notice] initPara: %f \n", initPara[1]);
        return theOp->objective(initPara);
    }
    double targetFun(unsigned n, const double* initPara, double *grad, void *op )
    {
        ObjectiveFunc* theOp  = (ObjectiveFunc*) (op);
        return theOp->objective(initPara);
    }
    double constrain(unsigned n, const double* initPara, double *grad, void *aCon )
    {
        Constraint* con  = (Constraint*) (aCon);
        return con->applyConstraint(initPara);
    }
}


/// A set of optimization parameters are specific to a objective function which
///   is the business logic of the optimization problem. This means Optimizer
///   has a otrong relationship with ObjectiveFunc.
/// Howto:
///  1. inherit Optimizer
///  2. Set the optimization parameters in doOptimize(), then invoke Optimizer::doOptimize().

class Optimizer
{

private:
    Optimizer();

public :
 
    void*      targetF;  /// Pointer to ObjectiveFunc object (abstract class). 
    inline     Optimizer(ObjectiveFunc* targetF):targetF(targetF){};

    nlopt_opt  opt;
    void doOptimize()
    {
        int errCode = 0;
        double  objVal = 0;
        ObjectiveFunc* theFunc = (ObjectiveFunc*) targetF;
        /// dirty way of converting vector to double* array.

        printf("[notice] pos initPara: %f \n", (theFunc->initParameters[0]));
        double* paras = &(theFunc->initParameters[0]);

        printf("[notice] initPara: %f \n", paras[0]);
        for(int i =0; i < theFunc->constraints.size(); i ++)
        {
            printf("adding %d\n", i );
            nlopt_add_inequality_constraint(opt, constrain, theFunc->constraints[i], 1e-8);
        }

        vector<double> paras_test = theFunc->initParameters;

        nlopt_set_max_objective(opt, targetFun_stub, theFunc);


        errCode = nlopt_optimize(opt, paras , &objVal);
        //errCode = nlopt_optimize(opt, &(paras_test[0]) , &objVal);
        if ( errCode < 0) {
            printf("[error] nlopt failed!  errorCode == %d\n", errCode);
            exit(0);
        }

        nlopt_destroy(opt);
    }

    friend double targetFun_stub(unsigned n, const double* initPara, double *grad, void *op );
    virtual ~Optimizer()
    {

    }
};

#endif

