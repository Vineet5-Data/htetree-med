/*
 * split.Rule = user (Median Absolute Deviation)
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"
#include "Medianfunction.h"

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;

int
userinit(int n, double *y[], int maxcat, char **error,
		int *size, int who, double *wt, double *treatment, 
		int bucketnum, int bucketMax, double *train_to_est_ratio)
{
	if (who == 1 && maxcat > 0) {
		graycode_init0(maxcat);
		countn = (int *) ALLOC(2 * maxcat, sizeof(int));
		tsplit = countn + maxcat;
		treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
		wts = treatment_effect + maxcat;
		trs = wts + maxcat;
		sums = trs + maxcat;
		wtsums = sums + maxcat;
		trsums = wtsums + maxcat;
		wtsqrsums = trsums + maxcat;
		trsqrsums = wtsqrsums + maxcat;
	}
	*size = 1;
	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
	return 0;
}


void
userss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, 
       double *risk, double *wt, double *treatment, double max_y,
       double alpha, double train_to_est_ratio)
{
	int i,r=0, k=0;
	double temp0 = 0., temp1 = 0., twt = 0., temp = 0.; /* sum of the weights */ 
	double ttreat = 0.;
	double effect;
	
	double tr_var, con_var;
	double con_sqr_sum = 0., tr_sqr_sum = 0.;
	
	double medianestimator; /* Hodges-Lehmann-Estimator */
  int tr_n = 0; 
  for(i = 0; i < n; i++) { /* # of treatment set  */
    tr_n+= treatment[i];
  }
  double y_treat[tr_n]; /* outcome variable, which received treatment */
  double y_con[(n-tr_n)]; /* outcome variable, which does not received treatment */

for (i = 0; i < n; i++) {
  if (treatment[i] == 0) {
  y_con[r]= *y[i] * wt[i]; /* fill entries of y_con*/
r++;
} else {
  y_treat[k]= *y[i]* wt[i]; /* fill entries of y_treat*/
k++;
}
}

	for (i = 0; i < n; i++) {
		temp1 += *y[i] * wt[i] * treatment[i]; 
		temp0 += *y[i] * wt[i] * (1 - treatment[i]);
		twt += wt[i];
		ttreat += wt[i] * treatment[i];
		
		tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
	}
//	medianestimator = testmean( y_treat,  ttreat) - testmean( y_con, (twt-ttreat));
//		medianestimator = quick_select5( y_treat,  ttreat) - quick_select5( y_con, (twt-ttreat));
  medianestimator = hodges_lehmann_estimator(y_treat,y_con,ttreat,(twt-ttreat)); 
// double ss=0;
//  for (i = 0; i < n; i++) {
  //  temp = fabs(( (*y[i] - medianestimator) * treatment[i]) + ( (*y[i] - medianestimator) * (1-treatment[i]) ));
//   ss += temp * wt[i];
// }
  
effect = temp1 / ttreat - temp0 / (twt - ttreat);
  
tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));
  
  *tr_mean = temp1 / ttreat;
  *con_mean = temp0 / (twt - ttreat);
  *value = effect;
// *risk = 4 * twt * max_y * max_y- alpha *ss;
 *risk = 4 * twt * max_y * max_y- alpha * fabs(effect * twt - medianestimator); // have to take a closer look 
}

void user(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, 
		int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
		double train_to_est_ratio)
{
	int i, j,m=0,p=0, r=0, k=0;
	double temp;
	double left_sum, right_sum;
	double left_tr_sum, right_tr_sum;
	double left_tr, right_tr;
	double left_wt, right_wt;
	int left_n, right_n;
	double best;
	int direction = LEFT;
	int where = 0;
	double node_effect, left_effect, right_effect;
	double left_temp, right_temp;
	int min_node_size = minsize;
	
	int tr_n = 0;
	for(i = 0; i < n; i++) {
	  tr_n+= treatment[i];
	}

	double tr_var, con_var;
	double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
	double left_tr_var, left_con_var, right_tr_var, right_con_var;
	
	double medianeffect_estimator, right_medianeffect_estimator, left_medianeffect_estimator;
	double y_right_all_tr[tr_n], y_right_all_con[n-tr_n], y_right_tr[tr_n], y_right_con[n-tr_n];
	double y_left_tr[tr_n];
	double y_left_con[n-tr_n];
 
	right_wt = 0.;
	right_tr = 0.;
	right_sum = 0.;
	right_tr_sum = 0.;

	right_sqr_sum = 0.;
	right_tr_sqr_sum = 0.;
	
	right_n = n;
	
	for (i = 0; i < n; i++) {
		right_wt += wt[i];
		right_tr += wt[i] * treatment[i];
		right_sum += *y[i] * wt[i];
		right_tr_sum += *y[i] * wt[i] * treatment[i];
		right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
		right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		
		if (treatment[i] == 0) {  /* Vielleicht ineffizient, muss ggf überarbeitet werden */
    y_right_all_con[r]= *y[i] * wt[i];
		  r++;
		} else {
		  y_right_all_tr[k]= *y[i] * wt[i];
		  k++;
		}
	}

for (i = 0; i < tr_n; i++) {
    y_right_tr[p]= y_right_all_tr[tr_n-1-i];
    p++;
    }
for (i = 0; i < n-tr_n; i++) {
  y_right_con[m]= y_right_all_con[(n-tr_n)-1-i];
  m++;
}

	temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
	tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
	con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
	  - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
	  / ((right_wt - right_tr) * (right_wt - right_tr));
	//medianeffect = findMedian(y_right_tr, right_tr) - findMedian(y_right_con, (right_wt - right_tr));
	
//	medianeffect_estimator = testmean( y_right_all_tr,  tr_n) - testmean( y_right_all_con, (n - tr_n));
	
//	medianeffect_estimator = quick_select5( y_right_all_tr,  tr_n) - quick_select5( y_right_all_con, (n - tr_n));
	medianeffect_estimator = hodges_lehmann_estimator(y_right_all_tr,y_right_all_con, tr_n,(n - tr_n));
	
	node_effect = alpha * fabs(temp * right_wt - medianeffect_estimator ); //??

	if (nclass == 0) {
		/* continuous predictor */
		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		best = 0;
		r=0;
		k=0;
		
		for (i = 0; right_n > edge; i++) {
			left_wt += wt[i];
			right_wt -= wt[i];
			
			left_tr += wt[i] * treatment[i];
			right_tr -= wt[i] * treatment[i];
			
			left_n++;
			right_n--;
			
			temp = *y[i] * wt[i] * treatment[i];
			
			left_tr_sum += temp;
			right_tr_sum -= temp;
			
			left_sum += *y[i] * wt[i];
			right_sum -= *y[i] * wt[i];
			
		    y_left_con[r]= y_right_all_con[i];
			  y_left_tr[k]= y_right_all_tr[i];
			 
			  r++;
			  k++;

			////////////////////
			if (x[i + 1] != x[i] && left_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size) {

				left_temp = left_tr_sum / left_tr - 
					(left_sum - left_tr_sum) / (left_wt - left_tr);
			
		
	//	left_medianeffect_estimator = quick_select5( y_left_tr,  left_tr) - quick_select5( y_left_con, (left_wt - left_tr));
	//	left_medianeffect_estimator = testmean( y_left_tr,  left_tr) - testmean( y_left_con, (left_wt - left_tr));
		
		left_medianeffect_estimator= hodges_lehmann_estimator(y_left_tr,y_left_con,left_tr,(left_wt - left_tr));
		left_effect = alpha * fabs(left_temp * left_wt -  left_medianeffect_estimator );

				right_temp = right_tr_sum / right_tr -
					(right_sum - right_tr_sum) / (right_wt - right_tr);
		
	//	right_medianeffect_estimator = testmean( y_right_tr,  right_tr) - testmean( y_right_con, (right_wt - right_tr));
//	right_medianeffect_estimator = quick_select5( y_right_tr,  right_tr) - quick_select5( y_right_con, (right_wt - right_tr));
		
				right_medianeffect_estimator = hodges_lehmann_estimator(y_right_tr,y_right_con,right_tr,(right_wt - right_tr));
					right_effect = alpha *fabs(right_temp * right_wt - right_medianeffect_estimator );

				temp = left_effect + right_effect - node_effect;
				if (temp > best) {
					best = temp;
					where = i;               
					if (left_temp < right_temp)
						direction = LEFT;
					else
						direction = RIGHT;
				}             

			}
		}

		*improve = best;
		if (best > 0) {         /* found something */
			csplit[0] = direction;
			*split = (x[where] + x[where + 1]) / 2; 
		}
	}

	/*
	 * Categorical predictor
	 */
	else {
	  for (i = 0; i < nclass; i++) {
			countn[i] = 0;
			wts[i] = 0;
			trs[i] = 0;
			sums[i] = 0;
			wtsums[i] = 0;
			trsums[i] = 0;
			wtsqrsums[i] = 0;
			trsqrsums[i] = 0;
		}

		/* rank the classes by treatment effect */
		for (i = 0; i < n; i++) {
			j = (int) x[i] - 1;
			countn[j]++;
			wts[j] += wt[i];
			trs[j] += wt[i] * treatment[i];
			sums[j] += *y[i];
			wtsums[j] += *y[i] * wt[i];
			trsums[j] += *y[i] * wt[i] * treatment[i];
			wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
			trsqrsums[j] +=  (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		}

		for (i = 0; i < nclass; i++) {
			if (countn[i] > 0) {
				tsplit[i] = RIGHT;
				treatment_effect[i] = trsums[j] / trs[j] - (wtsums[j] - trsums[j]) / (wts[j] - trs[j]);
			} else
				tsplit[i] = 0;
		}
		graycode_init2(nclass, countn, treatment_effect);

		/*
		 * Now find the split that we want
		 */

		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		left_sqr_sum = 0.;
		left_tr_sqr_sum = 0.;
		r=0;
		k=0;

		best = 0;
		where = 0;
		while ((j = graycode()) < nclass) {
			tsplit[j] = LEFT;
			left_n += countn[j];
			right_n -= countn[j];

			left_wt += wts[j];
			right_wt -= wts[j];

			left_tr += trs[j];
			right_tr -= trs[j];

			left_sum += wtsums[j];
			right_sum -= wtsums[j];

			left_tr_sum += trsums[j];
			right_tr_sum -= trsums[j];
			
			y_left_con[r]= y_right_all_con[i];
			y_left_tr[k]= y_right_all_tr[i];
			
			r++;
			k++;

			left_sqr_sum += wtsqrsums[j];
			right_sqr_sum -= wtsqrsums[j];

			left_tr_sqr_sum += trsqrsums[j];
			right_tr_sqr_sum -= trsqrsums[j];

			if (left_n >= edge && right_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size) {

				left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) 
					/ (left_wt - left_tr);
			  
	//		  left_medianeffect_estimator = testmean( y_left_tr,  left_tr) - testmean( y_left_con, (left_wt - left_tr));
//			 		  left_medianeffect_estimator = quick_select5( y_left_tr,  left_tr) - quick_select5( y_left_con, (left_wt - left_tr));

			  left_medianeffect_estimator= hodges_lehmann_estimator(y_left_tr,y_left_con,left_tr,(left_wt - left_tr));
			  
				left_tr_var = left_tr_sqr_sum / left_tr 
					- left_tr_sum  * left_tr_sum / (left_tr * left_tr);
				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
					/ ((left_wt - left_tr) * (left_wt - left_tr));       
				/*left_effect = alpha * left_temp * left_temp * left_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt * 
					(left_tr_var / left_tr + left_con_var / (left_wt - left_tr));*/
			left_effect = alpha * fabs(left_temp * left_wt -  left_medianeffect_estimator );
				
				right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) 
					/ (right_wt - right_tr);
		
	//	right_medianeffect_estimator = testmean( y_right_tr,  right_tr) - testmean( y_right_con, (right_wt - right_tr));
		
		//		right_medianeffect_estimator = quick_select5( y_right_tr,  right_tr) - quick_select5( y_right_con, (right_wt - right_tr));
				
				right_medianeffect_estimator = hodges_lehmann_estimator(y_right_tr,y_right_con,right_tr,(right_wt - right_tr));
				right_effect = alpha *fabs(right_temp * right_wt - right_medianeffect_estimator );
				
				right_tr_var = right_tr_sqr_sum / right_tr 
					- right_tr_sum * right_tr_sum / (right_tr * right_tr);
				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
					/ ((right_wt - right_tr) * (right_wt - right_tr));
				/*right_effect = alpha * right_temp * right_temp * right_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt *
					(right_tr_var / right_tr + right_con_var / (right_wt - right_tr)); */
				
				temp = left_effect + right_effect - node_effect;


				if (temp > best) {
					best = temp;

					if (left_temp > right_temp)
						for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
					else
						for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
				}
			}
		}
		*improve = best;
	}
}


double
userpred(double *y, double wt, double treatment, double *yhat, double propensity)
{
	double ystar;
	double temp;

	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
	temp = ystar - *yhat;
	return temp * temp * wt;
}

