#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#include <netcdf.h>
#include <omp.h>
#include "mex.h"


/* **********************************************************************************/
/* ****************************** ESTIMATE_DISTANCE *********************************/
/* **********************************************************************************/
/*
 * estimate_distance: estimate the distances between ref patch and moving patches (prev and after)
 * patches are of fixed size (2*SIM_PATCH_S_HAFTSIZE+1) x (2*SIM_PATCH_S_HAFTSIZE+1)
 * reference patch are centered at [center_ref_idy, center_ref_idz]
 * moving patches are centered at [center_moving_idy, center_moving_idz]
 * dist_all contain 2 elements: distances to moving patches in the prev and after plane
 * x_ref: reference plane
 * x_prev: previous plane
 * x_after: plane after
 * ref_ids_y(z): indices of points in reference patch
 * moving_ids_y(z): indices of points in moving patch
 */
void fusion(double *x_NLM, double *weight_NLM, double *x_ref, double *x_moving, double *x_fusion,
              int *gridpatches_y, int *gridpatches_z, int *acc_ids, int est_idy, int est_idz,
              int simpatch_haftsize, int neighbor_haftsize, int accpatch_fullsize, int Ny, int Nz, double tau)
{
    int dim_simpatch = 2*simpatch_haftsize + 1;
    int dim_neighbor = 2*simpatch_haftsize + 1;
    double norm_fact = 1.0/((double) (dim_simpatch*dim_simpatch));
    double dist[dim_neighbor*dim_neighbor];
    int ri = dim_neighbor*dim_neighbor/2;
    for (int i = 0; i < dim_neighbor; i++)
    {
        for (int j = 0; j < dim_neighbor; j++)
        {
            double d = 0.0;
            double du;
            int ni = i * dim_neighbor + j;

            int si, ref_idy, ref_idz, moving_idy, moving_idz;
            for (int m = 0; m < dim_simpatch; m++)
            {
                for (int n = 0; n < dim_simpatch; n++)
                {
                    si = i * dim_simpatch + j;
                    ref_idy = gridpatches_y[(est_idy * Nz + est_idz) * dim_neighbor * dim_simpatch + ri * dim_simpatch + si];
                    moving_idy = gridpatches_y[(est_idy * Nz + est_idz) * dim_neighbor * dim_simpatch + ni * dim_simpatch + si];
                    ref_idz = gridpatches_z[(est_idy * Nz + est_idz) * dim_neighbor * dim_simpatch + ri * dim_simpatch + si];
                    moving_idz = gridpatches_z[(est_idy * Nz + est_idz) * dim_neighbor * dim_simpatch + ni * dim_simpatch + si];

                    //compute distance btw reference patch and fusion patch
                    du = x_ref[ref_idy * Nz + ref_idz] - x_moving[moving_idy * Nz + moving_idz];
                    d = d + norm_fact*du*du;
                    //mexPrintf("\n d =  %f \n", d);
                    //mexEvalString("drawnow;");
                }
            }
            if(i==1)
                mexPrintf("\n d = %i \n", (est_idy * Nz + est_idz) * dim_neighbor * dim_simpatch + ri * dim_simpatch + si);

            /*
            double w = exp(-d/(2.0*tau*tau));
            for(int m = 0; m < accpatch_fullsize; m++)
            {
                int si =  acc_ids[m];
                int ref_idy = gridpatches_y[est_idy][est_idz][ri][si];
                int moving_idy = gridpatches_y[est_idy][est_idz][ni][si];
                int ref_idz = gridpatches_z[est_idy][est_idz][ri][si];
                int moving_idz = gridpatches_z[est_idy][est_idz][ni][si];

                x_NLM[ref_idy * Nz + ref_idz] = x_NLM[ref_idy * Nz + ref_idz] + w*x_fusion[moving_idy * Nz + moving_idz];
                weight_NLM[ref_idy * Nz + ref_idz] = weight_NLM[ref_idy * Nz + ref_idz] + w;
            }*/

        }

    }
}


/* **********************************************************************************/
/* ************************** MAIN FUNCTION MEXFUNCTION *****************************/
/* **********************************************************************************/
/*
 * all arguments should be as default
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    #define x_rec_OUT plhs[0]
    #define W_OUT plhs[1]

    #define x_ref_IN prhs[0]
    #define x_moving_IN prhs[1]
    #define x_fusion_IN prhs[2]
    #define gridpatches_y_IN prhs[3]
    #define gridpatches_z_IN prhs[4]
    #define simpatch_haftsize_IN prhs[5]
    #define neighbor_haftsize_IN prhs[6]
    #define acc_ids_IN prhs[7]
    #define tau_IN prhs[8]

    /* Get constants */
    int simpatch_haftsize = (int) mxGetScalar(simpatch_haftsize_IN);
    int neighbor_haftsize = (int) mxGetScalar(neighbor_haftsize_IN);
    double tau = (double) mxGetScalar(tau_IN);

    int simpatch_fullsize = 2 * simpatch_haftsize + 1;
    int neighbor_fullsize = 2 * neighbor_haftsize + 1;

    int Ny = (int) mxGetM(x_ref_IN);
    int Nz = (int) mxGetN(x_ref_IN);
    int accpatch_fullsize = mxGetNumberOfElements(acc_ids_IN);
    long int temp = mxGetNumberOfElements(acc_ids_IN);
    int dim_acc = accpatch_fullsize * accpatch_fullsize;
    int dim_sim = simpatch_fullsize * simpatch_fullsize;
    int dim_neighbor = neighbor_fullsize * neighbor_fullsize;

    /* Get the pointer to the data. Note that mxGetPr only get double type */
    double *x_ref = mxGetPr(x_ref_IN);
    double *x_moving = mxGetPr(x_moving_IN);
    double *x_fusion = mxGetPr(x_fusion_IN);
    int *acc_ids = (int *)mxGetPr(acc_ids_IN);
    int *gridpatches_y = (int *) mxGetPr(gridpatches_y_IN);
    int *gridpatches_z = (int *)mxGetPr(gridpatches_z_IN);

    /* Create the output matrix */
    x_rec_OUT = mxCreateDoubleMatrix(Ny, Nz, mxREAL);
    W_OUT = mxCreateDoubleMatrix(Ny, Nz, mxREAL);
    double *x_rec = mxGetPr(x_rec_OUT);
    double *W = mxGetPr(W_OUT);

    /*initialize super-resolved result and weight matrix with zeros */
    for (int k = 0; k < Ny * Nz; k++)
    {
        x_rec[k] = 0.0;
        W[k] = 0.0;
    }
    mexPrintf("\n gridpatches_y = %i \n", gridpatches_z[0]);
    mexPrintf("\n gridpatches_y = %i \n", gridpatches_z[1]);
    mexPrintf("\n gridpatches_y = %i \n", gridpatches_z[2]);


    //for (int i = 0; i < Ny + 1; i++)
        //for (int j = 0; j < Nz + 1; j++)
            //fusion(x_rec, W, x_ref, x_moving, x_fusion, gridpatches_y, gridpatches_z, acc_ids, i, j, simpatch_haftsize, neighbor_haftsize, accpatch_fullsize, Ny, Nz, tau);
}































/* **********************************************************************************/
/* ****************************** ACCUMULATION **************************************/
/* **********************************************************************************/

/*
 * accumulate_estimation: accumulate contribution in a small patch around center of moving patch [center_moving_idy, center_moving_idz]
 * accumulative patch is of size (2*ACCUM_PATCH_S_HAFTSIZE + 1) x (2*ACCUM_PATCH_S_HAFTSIZE + 1)
 * x_NLM: accumulative field
 * weight_NLM: weight matrix
 * x_prev, x_after: plane contain information to fuse
 * dist_all contains 2 distances of reference patch and moving patches from prev and after plane
 */

/*void accumulate_estimation(double *x_NLM, double *weight_NLM, double *x_fusion, double *dist
                           int *gridpatches_y, int *gridpatches_z, int *acc_ids, int accpatch_fullsize,
                           int simpatch_haftsize, int neighbor_haftsize, int est_idy, int est_idz, int Nz, double tau)
{
    double w = exp(-dist/(2.0*tau*tau));
    int dim_neighbor = 2*simpatch_haftsize + 1;
    int ri = dim_neighbor*dim_neighbor/2;
    for (int i = 0; i < dim_neighbor; i++)
    {
        for (int j = 0; j < dim_neighbor; j++)
        {
            int ni = i * dim_neighbor + j;
            double w = exp(-dist[i * dim_neighbor + j]/(2.0*tau*tau));

            for(int m = 0; m < accpatch_fullsize; m++)
            {
                int si =  acc_ids[m];
                int ref_idy = gridpatches_y[est_idy][est_idz][ri][si];
                int moving_idy = gridpatches_y[est_idy][est_idz][ni][si];
                int ref_idz = gridpatches_z[est_idy][est_idz][ri][si];
                int moving_idz = gridpatches_z[est_idy][est_idz][ni][si];

                x_NLM[ref_idy * Nz + ref_idz] = x_NLM[ref_idy * Nz + ref_idz] + w*x_fusion[moving_idy*NUM_COLS + moving_idz];
                weight_NLM[ref_idy*NUM_COLS + ref_idz] = weight_NLM[ref_idy*NUM_COLS + ref_idz] + w;
            }
    }
}*/
