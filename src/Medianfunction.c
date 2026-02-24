#include <math.h>
#include <stdio.h>
#include "Medianfunction.h"
#include <stdlib.h> // for qsort
#include <R_ext/Memory.h>
#define SWAP(x, y) {double temp = x; x = y; y = temp; }


/*int findmiddle(double a[],int n){
  
double m; 
// check for even case
if (n % 2 == 0){
  m = a[n/2];
} else{
  m = a[n/2+1];
}
return 0;

}
void swap(double *p,double *q) {
  double t;
  
  t=*p;
  *p=*q;
  *q=t;
}

*/


/*
void sort(double a[], int n) {
  int i,j;
  
  for(i = 0;i < n-1;i++) {
    for(j = 0;j < n-i-1;j++) {
      if(a[j] > a[j+1])
        swap(&a[j],&a[j+1]);
    }
  }
}
 
 */


/*void printArray(double a[], int size)
 {
 int i;
 for (i = 0; i < size; i++)
 printf("%f ", a[i]);
 printf("\n");
 }*/

/*
double findMedian(double a[], int n)
{ double m;
  // First we sort the array
  sort(a, n);
  // check for even case
  if (n % 2 == 0){
    m = (a[n/2] + a[n/2 - 1])/ 2.0;
  } else{
    m = a[n/2];
  }
  //printf("Median = %f ", m);
  return 0;
  
}

*/



// compare-function for qsort (changed to double instead of int)
static int cmpfunc(const void *a, const void *b) {
  const double da = *(const double*)a;
  const double db = *(const double*)b;
  return (da > db) - (da < db);
}

 
double findMedian(double a[], int n)
{ 
  // to store median
  double m; 
  // First we sort the array
  qsort(a, n, sizeof(double), cmpfunc);
  
  // check for even case
  if (n % 2 == 0){
    m = (a[n/2] + a[n/2 - 1])/ 2.0;
  } else{
    m = a[n/2];
  }
  //printf("Median = %f ", m);
 return m;
  
}
double quick_select5( double arr[],  int n)
{
   int low, high ;
   int median;
   int middle, ll, hh;
  
  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) /* One element only */
return arr[median] ;
    
    if (high == low + 1) {  /* Two elements only */
if (arr[low] > arr[high])
  SWAP(arr[low], arr[high]) ;
return arr[median] ;
    }
    
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     SWAP(arr[middle], arr[low]) ;
    
    /* Swap low item (now in position middle) into position (low+1) */
    SWAP(arr[middle], arr[low+1]) ;
    
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh]  > arr[low]) ;
      
      if (hh < ll)
        break;
      
      SWAP(arr[ll], arr[hh]) ;
    }
    
    /* Swap middle item (in position low) back into correct position */
    SWAP(arr[low], arr[hh]) ;
    
    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}

double hodges_lehmann_estimator(double *sample1, double *sample2,
                                int sample_1, int sample_2)
{
    int i, j;
    int samplesize = sample_1 * sample_2;        

    double *differences = (double*) R_alloc((size_t)samplesize, sizeof(double));

    int idx = 0;
    for (i = 0; i < sample_1; i++) {
        for (j = 0; j < sample_2; j++) {
            differences[idx++] = sample1[i] - sample2[j];
        }
    }

    return findMedian(differences, samplesize);
}

/*
double med_array_function(double *y[],int n, double *wt, double *treatment,double HL)
  {
  int i; 
  for (i=0; i < n; i++) {
    med_array_right_all[i] = (*y[i]) * (*y[i]) * wt[i] * treatment[i] + (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]) + (*y[i] * wt[i] * treatment[i] + *y[i] * wt[i] * (1- treatment[i])) * medianeffect_estimator + medianeffect_estimator * medianeffect_estimator; 
  }
  }

*/


/*
  
  // tested on https://www.tutorialspoint.com/compile_c_online.php
 
 */ 

/*
#include <stdio.h>
#include <stdlib.h>
 */
 
 /*
  // alternative way of comp func
  
  int cmpfunc (const void * a, const void * b)
  {
  if (*(double*)a > *(double*)b)
  return 1;
  else if (*(double*)a < *(double*)b)
  return -1;
  else
  return 0;  
  }
  */
 
 /*
 
 // used way of comp func
 int cmpfunc (const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
 }


double main () {
  int n = 6;
  double values[] = { 88.2, 56.3, 100.2, 2.0, 25.3, 30.2};
  
  double m;
  
  printf("Before sorting the list is: \n");
  for( n = 0 ; n < 6; n++ ) {
    printf("%f ", values[n]);
  }
  
  //int len = sizeof(values) /sizeof(values[0]);
  printf("\nLength: %d \n", n);
  
  
  qsort(values, 6, sizeof(double), cmpfunc);
  
  printf("\nAfter sorting the list is: \n");
  for( n = 0 ; n < 6; n++ ) {   
    printf("%f ", values[n]);
  }
  
  if (n % 2 == 0){
    m = (values[n/2] + values[n/2 - 1])/ 2.0;
  } else{
    m = values[n/2];
  }
  printf("\nMedian = %f \n", m);
  
  return(0);
}
 
 */



