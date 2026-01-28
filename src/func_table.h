/*
 * The table of implimented splitting functions
 *
 */

extern int totinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                   double *train_to_est_ratio);

extern void tot(int n, double *y[], double *x, int nclass,
                  int edge, double *improve, double *split, int *csplit,
                  double myrisk, double *wt, double *treatment, double propensity,int minsize);
extern void totss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                    double *wt, double *treatment, double max_y, double propensity);
extern double totpred(double *y, double wt, double treatment, double *yhat, double p);

extern int totDinit(int n, double *y[], int maxcat, char **error,
             int *size, int who, double *wt, double *treatment,
             int bucketnum, int bucketMax, double *train_to_est_ratio);
extern void totDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt,
                   double *treatment, double max_y, double propensity);
extern void totD(int n, double *y[], double *x, int nclass, int edge, double *improve,
                      double *split, int *csplit, double myrisk, double *wt,
                      double *treatment, double propensity, int minsize, int bucketnum,
                      int bucketMax);
extern double totDpred(double *y, double wt, double treatment, double *yhat, double p);

extern int CTinit(int n, double *y[], int maxcat, char **error,
           int *size, int who, double *wt, double *treatment, int bucketnum,
           int bucketMax, double *train_to_est_ratio);
extern void CTss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                 double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void CT(int n, double *y[], double *x, int nclass,int edge, double *improve, double *split,
               int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
               double train_to_est_ratio);
extern double CTpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int CTDinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                   double *train_to_est_ratio);

extern void CTDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt,
                  double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void CTD(int n, double *y[], double *x, int nclass,
        int edge, double *improve, double *split, int *csplit,
        double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
        int bucketMax, double train_to_est_ratio);
extern double CTDpred(double *y, double wt, double treatment, double *yhat,
                      double propensity);

extern int fitinit(int n, double *y[], int maxcat, char **error,
            int *size, int who, double *wt, double *treatment, int bucketnum,
            int bucketMax, double *train_to_est_ratio);

extern void fitss(int n, double *y[], double *value, double *tr_mean, double *con_mean,
                  double *risk, double *wt, double *treatment, double max_y,
                  double alpha, double train_to_est_ratio);

extern void fit(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize,
                double alpha, double train_to_est_ratio);
extern double fitpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int fitDinit(int n, double *y[], int maxcat, char **error,
                   int *size, int who, double *wt, double *treatment, int bucketnum,
                   int bucketMax, double *train_to_est_ratio);

extern void fitDss(int n, double *y[], double *value, double *tr_mean, double *con_mean,
                  double *risk, double *wt, double *treatment, double max_y,
                  double alpha, double train_to_est_ratio);

extern void fitD(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize,
                int bucketnum, int bucketMax, double alpha, double train_to_est_ratio);
extern double fitDpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int tstatsinit(int n, double *y[], int maxcat, char **error,
                      int *size, int who, double *wt, double *treatment,
                      int bucketnum, int bucketMax, double *train_to_est_ratio);

extern void tstatsss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                     double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void tstats(int n, double *y[], double *x, int nclass,
                   int edge, double *improve, double *split, int *csplit,
                   double myrisk, double *wt, double *treatment, int minsize,
                   double alpha, double train_to_est_ratio);
extern double tstatspred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int tstatsDinit(int n, double *y[], int maxcat, char **error,
                int *size, int who, double *wt, double *treatment,
                int bucketnum, int bucektMax, double *train_to_est_ratio);
extern void tstatsDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                      double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void tstatsD(int n, double *y[], double *x, int nclass,
            int edge, double *improve, double *split, int *csplit,
            double myrisk, double *wt, double *treatment, int minsize, double alpha,
            int bucketnum, int bucketMax, double train_to_est_ratio);
extern double tstatsDpred(double *y, double wt, double treatment, double *yhat, double propensity);


extern int usersplit_init(int n, double *y[], int maxcat, char **error,
                  int *size, int who, double *wt, double *treatment, int bucketnum,
                  int bucketMax, double *train_to_est_ratio);
extern void usersplit_eval(int n, double *y[], double *value, double *tr_mean, double *con_mean,
                 double *risk, double *wt, double *treatment, double max_y, double alpha,
                 double train_to_est_ratio);
extern void usersplit(int n, double *y[], double *x, int nclass,int edge, double *improve,
               double *split, int *csplit, double myrisk, double *wt, double *treatment,
               int minsize, double alpha, double train_to_est_ratio);
extern double usersplit_pred(double *y, double wt, double treatment, double *yhat, double propensity);


extern int userDinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                   double *train_to_est_ratio);

extern void userDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt,
                  double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void userD(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
                int bucketMax, double train_to_est_ratio);
extern double userDpred(double *y, double wt, double treatment, double *yhat,
                      double propensity);

///// adding policy
extern int policyinit(int n, double *y[], int maxcat, char **error,
                    int *size, int who, double *wt, double *treatment, int bucketnum,
                    int bucketMax, double *train_to_est_ratio);
extern void policyss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                   double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void policy(int n, double *y[], double *x, int nclass,int edge, double *improve, double *split,
                 int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
                 double train_to_est_ratio);
extern double policypred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int policyDinit(int n, double *y[], int maxcat, char **error, int *size,
                     int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                     double *train_to_est_ratio);

extern void policyDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt,
                    double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void policyD(int n, double *y[], double *x, int nclass,
                  int edge, double *improve, double *split, int *csplit,
                  double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
                  int bucketMax, double train_to_est_ratio);
extern double policyDpred(double *y, double wt, double treatment, double *yhat,
                        double propensity);


extern int medinit(int n, double *y[], int maxcat, char **error,
                 int *size, int who, double *wt, double *treatment, int bucketnum,
                 int bucketMax, double *train_to_est_ratio);
extern void medss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void med(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split,
              int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
              double train_to_est_ratio);
extern double medpred(double *y, double wt, double treatment, double *yhat, double propensity);


extern int medDinit(int n, double *y[], int maxcat, char **error,
                 int *size, int who, double *wt, double *treatment, int bucketnum,
                 int bucketMax, double *train_to_est_ratio);
extern void medDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void medD(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split,
              int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
              double train_to_est_ratio);
extern double medDpred(double *y, double wt, double treatment, double *yhat, double propensity);


/* ------------------------------------------------------ --------------------------- */
extern double tot_xpred(double *y, double wt, double treatment, double *yhat, double propensity);
extern double matching_xpred(double *y1, double *y2, double wt1, double wt2, double treatment1,
                                  double treatment2, double pred1, double pred2);
extern double fitH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                         double trs, double cons, double alpha, double xtrain_to_est_ratio);
extern double fitA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean);
extern double CTH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double trs, double cons, double alpha, double xtrain_to_est_ratio,
                        double propensity);
extern double CTA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                 double tree_tr_mean, double tree_con_mean, double alpha);
extern double userH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double trs, double cons, double alpha, double xtrain_to_est_ratio,
                        double propensity);
extern double userA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double tree_tr_mean, double tree_con_mean, double alpha);

///// adding policy
extern double policyH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                          double trs, double cons, double alpha, double xtrain_to_est_ratio,
                          double propensity);
extern double policyA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                          double tree_tr_mean, double tree_con_mean, double alpha, double gamma);


//extern double totxeval(int *unique_leaf, int **val_leaf_mat, int cp_id, int t, int *sorts, double *wt,
//                      double *treatment,  double *y[], double propensity, int k, int nobs,
//                       double val_sum_wt, int val_count);



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
extern int anovainit(int n, double *y[], int maxcat, char **error,
  	     double *parm, int *size, int who, double *wt, double *treatment);


extern void anovass(int n, double *y[], double *value, double *risk,
  	    double *wt, double *treatment, double max_y);


extern void anova(int n, double *y[], double *x, int nclass,
  	  int edge, double *improve, double *split, int *csplit,
		  double myrisk, double *wt, double *treatment, int minsize);

extern double anovapred(double *y, double wt, double treatment, double *yhat, double p);


static struct {
    int (*init_split) (int n, double *y[], int maxcat, char **error,
         double *parm, int *size, int who, double *wt, double *treatment);
    void (*choose_split) (int n, double *y[], double *x, int nclass,
          int edge, double *improve, double *split, int *csplit,
          double myrisk, double *wt, double *treatment, int minsize);
    void (*eval) (int n, double *y[], double *value, double *risk,
          double *wt, double *treatment, double max_y);
    double (*error) (char);
} func_table[] = {
    {anovainit, anova, anovass, (double (*)(char))anovapred}, // test for TOT:
};
#define NUM_METHODS 1        /* size of the above structure */

static struct {
    int (*init_split) (int n, double *y[], int maxcat, char **error,int *size, int who, double *wt, double *treatment,int bucketnum, int bucketMax, double *train_to_est_ratio);
    void (*choose_split) (int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...);
    void (*eval) (int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...);
    double (*error) (char);
} split_func_table[] = { // important changes in medinit
    {totinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))tot, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))totss, (double (*)(char))totpred},
    {CTinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))CT, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))CTss, (double (*)(char))CTpred},
    {fitinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))fit, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))fitss, (double (*)(char))fitpred},
    {tstatsinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))tstats, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))tstatsss, (double (*)(char))tstatspred},
    {totDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))totD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))totDss, (double (*)(char))totDpred},
    {CTDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))CTD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))CTDss, (double (*)(char))CTDpred},
    {fitDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))fitD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))fitDss, (double (*)(char))fitDpred},
    {tstatsDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))tstatsD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))tstatsDss, (double (*)(char))tstatsDpred},
    {usersplit_init, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))usersplit, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))usersplit_eval, (double (*)(char))usersplit_pred},
    {userDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))userD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))userDss, (double (*)(char))userDpred},
    {policyinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))policy, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))policyss, (double (*)(char))policypred},
    {policyDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))policyD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))policyDss, (double (*)(char))policyDpred},
    {medinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))med, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))medss, (double (*)(char))medpred},
    {medDinit, (void (*)(int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...))medD, (void (*)(int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...))medDss, (double (*)(char))medDpred},
};

static struct {
        double (* xeval)(double *y, ...);
} cv_func_table[] = {
    {(double (*)(double *y, ...))tot_xpred},
    {(double (*)(double *y, ...))matching_xpred},
    {(double (*)(double *y, ...))fitH_xpred},
    {(double (*)(double *y, ...))fitA_xpred},
    {(double (*)(double *y, ...))CTH_xpred},
    {(double (*)(double *y, ...))CTA_xpred},
    {(double (*)(double *y, ...))userH_xpred},
    {(double (*)(double *y, ...))userA_xpred},
    {(double (*)(double *y, ...))policyH_xpred},
    {(double (*)(double *y, ...))policyA_xpred},
};

#define NUM_SPLIT_RULE 14 // total number of split rules
#define NUM_CROSSMETH 11
