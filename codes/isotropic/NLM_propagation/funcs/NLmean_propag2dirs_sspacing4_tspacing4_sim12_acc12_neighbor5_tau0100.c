/*
 * compile: gcc -O3 -std=c99 -o [filename_out] -fopenmp [filename].c -lm -I/usr/include/netcdf-3/ -L/usr/lib64/ -lnetcdf -lnetcdf_c++
 * in the terminal: export OMP_NUM_THREADS=3
*/

#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <omp.h>

/* This is the name of the data file we will read. */
#define FILENAME_RD "/data/PhDworks/isotropic/NLM/Udiff_spacespacing4.nc"
#define FILENAME_WR "/data/PhDworks/isotropic/NLM/NLmean_propag2dirs_sspacing4_tspacing4_sim12_acc12_neighbor5_tau0100.nc"

/* all constants */
#define N_HR 96

#define SCALE_FACTOR_SPACE 4
#define SCALE_FACTOR_TIME 4

#define SIM_HAFTSIZE 12
#define ACC_HAFTSIZE 12
#define NEIGHBOR_HAFTSIZE 5

#define SIM_FULLSIZE (2 * SIM_HAFTSIZE + 1)
#define ACC_FULLSIZE (2 * ACC_HAFTSIZE + 1)
#define NEIGHBOR_FULLSIZE (2 * NEIGHBOR_HAFTSIZE + 1)

#define TAU 0.1

#define NUM_VARS 1
#define NUM_SCALES 2


#define NUM_3DSNAPS 37  /* #3D snapshots */
#define NUM_BLOCKS N_HR/SCALE_FACTOR_TIME - 1 /* #(1:SCALE_FACTOR_TIME:N_HR) - 1*/
#define NUM_2DSNAPS (SCALE_FACTOR_TIME * NUM_BLOCKS + 1) /* #2D snapshots in each 3D block */
#define NDIMS 4

/* Handle errors by printing an error message and exiting with a non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

/* **********************************************************************************/
/* ****************************** USEFUL FUNCTIONS **********************************/
/* **********************************************************************************/

/*
 * get_onesnap: take part of a big array(arr1) and put to small one (arr2): arr2 = arr1[id_start:id_end]
 */
void get_onesnap(double *arr1,double *arr2, int id_start, int id_end)
{
    for (int i = id_start; i < id_end + 1; i++)
        arr2[i - id_start] = arr1[i];
}


/*
 * put_onesnap: assign small array (arr2) into biger one (arr1): arr1[id_start:id_end] = arr2
 */
void put_onesnap(double *arr1,double *arr2, int id_start, int id_end)
{
    for (int i = id_start; i < id_end + 1; i++)
        arr1[i] = arr2[i - id_start];
}


/*
 * norm_by_weight: normalize x[dim] by weight W[dim]
 */

void norm_by_weight(int dim, double *x, double *W)
{
    for (int k = 0; k < dim; k++)
        x[k] = x[k]/W[k];
}

void add_mat(int dim, double *sum, double *x1, double *x2)
{
    for (int k = 0; k < dim; k++)
        sum[k] = x1[k] + x2[k];
}

void initialize(int dim, double *x, double val)
{
    for (int k = 0; k < dim; k++)
        x[k] = val;
}

/* **********************************************************************************/
/* ****************************** NETCDF UTILS **************************************/
/* **********************************************************************************/

/*
 * creat_netcdf: create the netcdf file [filename] contain [num_vars] variables
 * variable names are [varname]
*/

void create_netcdf(char *filename, int num_vars, char *varname[num_vars])
{
    int ncid_wr, retval_wr;
    int vel_varid_wr;
    int Nt, Nx, Ny, Nz;
    int dimids[NDIMS];

    /* Create the file. */
    if ((retval_wr = nc_create(filename, NC_CLOBBER, &ncid_wr)))
       ERR(retval_wr);

    /* Define the dimensions. The record dimension is defined to have
     * unlimited length - it can grow as needed.*/
    if ((retval_wr = nc_def_dim(ncid_wr, "Ny", N_HR, &Ny)))
        ERR(retval_wr);
    if ((retval_wr = nc_def_dim(ncid_wr, "Nz", N_HR, &Nz)))
        ERR(retval_wr);
    if ((retval_wr = nc_def_dim(ncid_wr, "Nt", NC_UNLIMITED, &Nt)))
        ERR(retval_wr);

    /* Define the netCDF variables for the data. */
    dimids[0] = Nt;
    dimids[1] = Nx;
    dimids[2] = Ny;
    dimids[3] = Nz;

    for (int i = 0; i<num_vars; i++)
    {
        if ((retval_wr = nc_def_var(ncid_wr, varname[i], NC_FLOAT, NDIMS, dimids, &vel_varid_wr)))
            ERR(retval_wr);
    }

    /* End define mode (SHOULD NOT FORGET THIS!). */
    if ((retval_wr = nc_enddef(ncid_wr)))
        ERR(retval_wr);

    /* Close the file. */
    if ((retval_wr = nc_close(ncid_wr)))
        ERR(retval_wr);
    printf("\n *** SUCCESS creating file: %s!\n", filename);
}


/*
 * write_netcdf:
 * write into [filename], variable [varname] [snap_end - snap_start + 1 ] snapshots [snaps] started at [snap_start]
*/

void write_netcdf(char *filename, char *varname, size_t *start, size_t *count, double *snaps)
{
    int ncid_wr, retval_wr;
    int vel_varid_wr;

    /* Open the file. NC_WRITE tells netCDF we want read-only access to the file.*/
    if ((retval_wr = nc_open(filename, NC_WRITE, &ncid_wr)))
        ERR(retval_wr);

    /* Get variable*/
    if ((retval_wr = nc_inq_varid(ncid_wr, varname, &vel_varid_wr)))
        ERR(retval_wr);;

    /* Put variable*/
    if ((retval_wr = nc_put_vara_double(ncid_wr, vel_varid_wr, start, count, &snaps[0])))
        ERR(retval_wr);

    /* Close the file. */
    if ((retval_wr = nc_close(ncid_wr)))
        ERR(retval_wr);

    printf("\n *** SUCCESS writing variables \"%s\" to \"%s\"!\n", varname, filename);
}


/*
 * read_netcdf: read from [filename], variable [varname] [snap_end - snap_start + 1 ] snapshots [snaps]
 * started at [snap_start]
*/

void read_netcdf(char *filename, char *varname, size_t *start, size_t *count, double *snaps)
{
    int ncid_rd, retval_rd;
    int vel_varid_rd;

    /*  ******** PREPARE TO READ ************* */
    /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
    if ((retval_rd = nc_open(filename, NC_NOWRITE, &ncid_rd)))
        ERR(retval_rd);

    /* Get the varids of the velocity in netCDF */
    if ((retval_rd = nc_inq_varid(ncid_rd, varname, &vel_varid_rd)))
        ERR(retval_rd);

    if ((retval_rd = nc_get_vara_double(ncid_rd, vel_varid_rd, start, count, &snaps[0])))
        ERR(retval_rd);

    /* Close the file, freeing all resources. */
    if ((retval_rd = nc_close(ncid_rd)))
        ERR(retval_rd);

    printf("\n *** SUCCESS reading variables \"%s\" from \"%s\" \n", varname, filename);
}


/* **********************************************************************************/
/* ****************************** ESTIMATE_DISTANCE *********************************/
/* **********************************************************************************/
/*
 * estimate_distance: estimate the distances between ref patch and moving patches (prev and after)
 * patches are of fixed size (2*SIM_HAFTSIZE+1) x (2*SIM_HAFTSIZE+1)
 * reference patch are centered at [center_ref_idy, center_ref_idz]
 * moving patches are centered at [center_moving_idy, center_moving_idz]
 * dist_all contain 2 elements: distances to moving patches in the prev and after plane
 * x_ref: reference plane
 * x_prev: previous plane
 * x_after: plane after
 * ref_ids_y(z): indices of points in reference patch
 * moving_ids_y(z): indices of points in moving patch
 */
void generate_grids(int *gridpatches_y, int *gridpatches_z, int * acc_ids)
{
    int neighbor_id, sim_id;
    int gridyoffset_neighbor[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE], gridzoffset_neighbor[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE];

    for (int m = 0; m < NEIGHBOR_FULLSIZE; m++)
    {
        for (int n = 0; n < NEIGHBOR_FULLSIZE; n++)
        {
            gridyoffset_neighbor[m * NEIGHBOR_FULLSIZE + n] = m - NEIGHBOR_HAFTSIZE;
            gridzoffset_neighbor[m * NEIGHBOR_FULLSIZE + n] = n - NEIGHBOR_HAFTSIZE;
        }
    }

    int gridyoffset_sim[SIM_FULLSIZE * SIM_FULLSIZE], gridzoffset_sim[SIM_FULLSIZE * SIM_FULLSIZE];
    for (int p = 0; p < SIM_FULLSIZE; p++)
    {
        for (int q = 0; q < SIM_FULLSIZE; q++)
        {
            gridyoffset_sim[p * SIM_FULLSIZE + q] = p - SIM_HAFTSIZE;
            gridzoffset_sim[p * SIM_FULLSIZE + q] = q - SIM_HAFTSIZE;
        }
    }

    int grid_sim[SIM_FULLSIZE][SIM_FULLSIZE];
    for (int p = 0; p < SIM_FULLSIZE; p++)
        for (int q = 0; q < SIM_FULLSIZE; q++)
            grid_sim[p][q] = p * SIM_FULLSIZE + q;
    for (int p = 0; p < ACC_FULLSIZE; p++)
        for (int q = 0; q < ACC_FULLSIZE; q++)
            acc_ids[p * ACC_FULLSIZE + q] = grid_sim[SIM_HAFTSIZE - ACC_HAFTSIZE + p][SIM_HAFTSIZE - ACC_HAFTSIZE + q];

    int valy, valz;
    long int grid_id;
    for (int i = 0; i < N_HR; i++)
    {
        for (int j = 0; j < N_HR; j++)
        {
            for (int neighbor_id = 0; neighbor_id < NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE; neighbor_id++)
            {
                for (int sim_id = 0; sim_id < SIM_FULLSIZE * SIM_FULLSIZE; sim_id++)
                {
                    grid_id = i * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                            + j * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                            + neighbor_id * SIM_FULLSIZE * SIM_FULLSIZE + sim_id;

                    valy = i + gridyoffset_neighbor[neighbor_id] + gridyoffset_sim[sim_id];
                    valz = j + gridzoffset_neighbor[neighbor_id] + gridzoffset_sim[sim_id];

                    if (valy < 0)
                        gridpatches_y[grid_id] = (N_HR - 1) + valy;
                    else if (valy > (N_HR - 1))
                        gridpatches_y[grid_id] = valy - (N_HR - 1);
                    else
                        gridpatches_y[grid_id] = valy;

                    if (valz < 0)
                        gridpatches_z[grid_id] = (N_HR - 1) + valz;
                    else if (valz > (N_HR - 1))
                        gridpatches_z[grid_id] = valz - (N_HR - 1);
                    else
                        gridpatches_z[grid_id] = valz;
                }
            }
        }
    }
    //printf("\n  gridpatches_z: %i \n", gridpatches_y[0]);
}


/* **********************************************************************************/
/* ****************************** NLMEAN *********************************/
/* **********************************************************************************/
/*
 * estimate_distance: estimate the distances between ref patch and moving patches (prev and after)
 * patches are of fixed size (2*SIM_HAFTSIZE+1) x (2*SIM_HAFTSIZE+1)
 * reference patch are centered at [center_ref_idy, center_ref_idz]
 * moving patches are centered at [center_moving_idy, center_moving_idz]
 * dist_all contain 2 elements: distances to moving patches in the prev and after plane
 * x_ref: reference plane
 * x_prev: previous plane
 * x_after: plane after
 * ref_ids_y(z): indices of points in reference patch
 * moving_ids_y(z): indices of points in moving patch
 */
/*void fusion(double *x_NLM, double *weight_NLM, double *x_ref, double *x_moving, double *x_fusion,
              int gridpatches_y[N_HR][N_HR][NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE][SIM_FULLSIZE * SIM_FULLSIZE],
              int gridpatches_z[N_HR][N_HR][NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE][SIM_FULLSIZE * SIM_FULLSIZE],
              int acc_ids[NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE], int est_idy, int est_idz)*/
void NLmean(double *x_NLM, double *weight_NLM, double *x_ref, double *x_moving, double *x_fusion, int *gridy, int *gridz, int *accids)
{
    double norm_fact = 1.0/((double) (SIM_FULLSIZE * SIM_FULLSIZE));
    int ri = NEIGHBOR_HAFTSIZE * NEIGHBOR_FULLSIZE + NEIGHBOR_HAFTSIZE;

    int est_idy;
    #pragma omp parallel for private (est_idy)
    for (est_idy = 0; est_idy < N_HR; est_idy++)
        for (int est_idz = 0; est_idz < N_HR; est_idz++)
            for (int ni = 0; ni < NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE; ni++)
            {
                int ref_idy, ref_idz, moving_idy, moving_idz;
                double du;
                double d = 0.0;
                long int grid_rid, grid_nid;

                for (int si = 0; si < SIM_FULLSIZE * SIM_FULLSIZE; si++)
                {
                    grid_rid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ri * SIM_FULLSIZE * SIM_FULLSIZE + si ;
                    grid_nid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ni * SIM_FULLSIZE * SIM_FULLSIZE + si;

                    ref_idy = gridy[grid_rid];
                    moving_idy = gridy[grid_nid];
                    ref_idz = gridz[grid_rid];
                    moving_idz = gridz[grid_nid];

                    //compute distance btw reference patch and fusion patch
                    du = x_ref[ref_idy * N_HR + ref_idz] - x_moving[moving_idy * N_HR + moving_idz];
                    d = d + norm_fact*du*du;
                }

                double w = exp(-d/(2.0*TAU*TAU));
                for(int k = 0; k < ACC_FULLSIZE * ACC_FULLSIZE; k++)
                {
                    int ai =  accids[k];
                    grid_rid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ri * SIM_FULLSIZE * SIM_FULLSIZE + ai ;
                    grid_nid = est_idy * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE
                             + est_idz * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE + ni * SIM_FULLSIZE * SIM_FULLSIZE + ai;

                    ref_idy = gridy[grid_rid];
                    moving_idy = gridy[grid_nid];
                    ref_idz = gridz[grid_rid];
                    moving_idz = gridz[grid_nid];

                    x_NLM[ref_idy * N_HR + ref_idz] = x_NLM[ref_idy * N_HR + ref_idz] + w*x_fusion[moving_idy * N_HR + moving_idz];
                    weight_NLM[ref_idy * N_HR + ref_idz] = weight_NLM[ref_idy * N_HR + ref_idz] + w;
                }
                //printf("\n w=%f\n ",w);
            }
}




void propag_forward(double *Xrec, double *Xlf, int *gridy, int *gridz, int *accids, int t_first, int t_bound1, int t_offset)
{
    for (int t_est = t_first + 1; t_est <= t_bound1; t_est++)
    {
        int t_prev = t_est - 1;
        double xref_lf[N_HR * N_HR], xref_hf[N_HR * N_HR], xmov_lf[N_HR * N_HR], xmov_hf[N_HR * N_HR], w[N_HR * N_HR];
        get_onesnap(Xlf, xref_lf, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);
        get_onesnap(Xlf, xmov_lf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);

        //Initialize with zeros
        initialize(N_HR * N_HR, xref_hf, 0.0);
        initialize(N_HR * N_HR, w, 0.0);

        // Propagation from previous planes
        NLmean(xref_hf, w, xref_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        // Normalize and put back
        norm_by_weight(N_HR*N_HR, xref_hf, w);
        put_onesnap(Xrec, xref_hf, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);
    }
}

void propag_backward(double *Xrec, double *Xlf, int *gridy, int *gridz, int *accids, int t_last, int t_bound2, int t_offset)
{
    for (int t_est = t_last - 1; t_est >= t_bound2; --t_est)
    {
        int t_prev = t_est + 1;
        double xref_lf[N_HR * N_HR], xref_hf[N_HR * N_HR], xmov_lf[N_HR * N_HR], xmov_hf[N_HR * N_HR], w[N_HR * N_HR];

        get_onesnap(Xlf, xref_lf, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);
        get_onesnap(Xlf, xmov_lf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);

        //Initialize with zeros
        initialize(N_HR * N_HR, xref_hf, 0.0);
        initialize(N_HR * N_HR, w, 0.0);

        // Propagation from previous planes
        NLmean(xref_hf, w, xref_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        // Normalize and put back
        norm_by_weight(N_HR*N_HR, xref_hf, w);
        put_onesnap(Xrec, xref_hf, t_offset + t_est * N_HR * N_HR, t_offset + (t_est + 1) * N_HR * N_HR - 1);
    }
}


void propag_2planes(double *Xrec, double *Xlf, int *gridy, int *gridz, int *accids, int t_mid, int t_offset)
{
    double xref_lf[N_HR * N_HR], xref_hf[N_HR * N_HR], xmov_lf[N_HR * N_HR], xmov_hf[N_HR * N_HR], w[N_HR * N_HR];
    int t_prev = t_mid - 1;
    int t_after = t_mid + 1;
    //Initialize with zeros
    initialize(N_HR * N_HR, xref_hf, 0.0);
    initialize(N_HR * N_HR, w, 0.0);

    get_onesnap(Xlf, xref_lf, t_offset + t_mid * N_HR * N_HR, t_offset + (t_mid + 1) * N_HR * N_HR - 1);

    get_onesnap(Xlf, xmov_lf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
    get_onesnap(Xrec, xmov_hf, t_offset + t_prev * N_HR * N_HR, t_offset + (t_prev + 1) * N_HR * N_HR - 1);
    NLmean(xref_hf, w, xref_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

    get_onesnap(Xlf, xmov_lf, t_offset + t_after * N_HR * N_HR, t_offset + (t_after + 1) * N_HR * N_HR - 1);
    get_onesnap(Xrec, xmov_hf, t_offset + t_after * N_HR * N_HR, t_offset + (t_after + 1) * N_HR * N_HR - 1);
    NLmean(xref_hf, w, xref_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

    // Normalize and put back
    norm_by_weight(N_HR*N_HR, xref_hf, w);
    put_onesnap(Xrec, xref_hf, t_offset + t_mid * N_HR * N_HR, t_offset + (t_mid + 1) * N_HR * N_HR - 1);
}



void propag_towardcenter(double *Xrec, double *Xlf, int *gridy, int *gridz, int *accids, int t_first, int t_offset)
{
    double xref1_lf[N_HR * N_HR], xref2_lf[N_HR * N_HR], xmov_lf[N_HR * N_HR], xmov_hf[N_HR * N_HR];
    double xref1_hf[N_HR * N_HR], w1[N_HR * N_HR], xref2_hf[N_HR * N_HR], w2[N_HR * N_HR];

    int tc = (int)SCALE_FACTOR_TIME/2;
    if (SCALE_FACTOR_TIME % 2) { tc = (int)SCALE_FACTOR_TIME/2 + 1; }

    for (int td = 1; td < tc; td++)
    {
        int t1 = t_first + td; // bound on left side
        int t2 = t_first + SCALE_FACTOR_TIME - td; // bound on right side

        // Initialize with zeros
        initialize(N_HR * N_HR, xref1_hf, 0.0);
        initialize(N_HR * N_HR, w1, 0.0);
        initialize(N_HR * N_HR, xref2_hf, 0.0);
        initialize(N_HR * N_HR, w2, 0.0);

        get_onesnap(Xlf, xref1_lf, t_offset + t1 * N_HR * N_HR, t_offset + (t1 + 1) * N_HR * N_HR - 1);
        get_onesnap(Xlf, xref2_lf, t_offset + t2 * N_HR * N_HR, t_offset + (t2 + 1) * N_HR * N_HR - 1);

        //Propagate from left bound
        get_onesnap(Xlf, xmov_lf, t_offset + (t1 - 1) * N_HR * N_HR, t_offset + t1 * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + (t1 - 1) * N_HR * N_HR, t_offset + t1 * N_HR * N_HR - 1);
        NLmean(xref1_hf, w1, xref1_lf, xmov_lf, xmov_hf, gridy, gridz, accids);
        NLmean(xref2_hf, w2, xref2_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        //Propagate from right bound
        get_onesnap(Xlf, xmov_lf, t_offset + (t2 + 1) * N_HR * N_HR, t_offset + (t2 + 2) * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + (t2 + 1) * N_HR * N_HR, t_offset + (t2 + 2) * N_HR * N_HR - 1);
        NLmean(xref1_hf, w1, xref1_lf, xmov_lf, xmov_hf, gridy, gridz, accids);
        NLmean(xref2_hf, w2, xref2_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        // Normalize and put back
        norm_by_weight(N_HR*N_HR, xref1_hf, w1);
        put_onesnap(Xrec, xref1_hf, t_offset + t1 * N_HR * N_HR, t_offset + (t1 + 1) * N_HR * N_HR - 1);
        norm_by_weight(N_HR*N_HR, xref2_hf, w2);
        put_onesnap(Xrec, xref2_hf, t_offset + t2 * N_HR * N_HR, t_offset + (t2 + 1) * N_HR * N_HR - 1);
    }

    // Last plane in the center
    if (SCALE_FACTOR_TIME % 2 == 0)
    {
        initialize(N_HR * N_HR, xref1_hf, 0.0);
        initialize(N_HR * N_HR, w1, 0.0);

        get_onesnap(Xlf, xref1_lf, t_offset + (t_first + tc) * N_HR * N_HR, t_offset + (t_first + tc + 1) * N_HR * N_HR - 1);
        get_onesnap(Xlf, xmov_lf, t_offset + (t_first + tc - 1) * N_HR * N_HR, t_offset + (t_first + tc) * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + (t_first + tc - 1) * N_HR * N_HR, t_offset + (t_first + tc) * N_HR * N_HR - 1);
        NLmean(xref1_hf, w1, xref1_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        get_onesnap(Xlf, xmov_lf, t_offset + (t_first + tc + 1) * N_HR * N_HR, t_offset + (t_first + tc + 2) * N_HR * N_HR - 1);
        get_onesnap(Xrec, xmov_hf, t_offset + (t_first + tc + 1) * N_HR * N_HR, t_offset + (t_first + tc + 2) * N_HR * N_HR - 1);
        NLmean(xref1_hf, w1, xref1_lf, xmov_lf, xmov_hf, gridy, gridz, accids);

        norm_by_weight(N_HR*N_HR, xref1_hf, w1);
        put_onesnap(Xrec, xref1_hf, t_offset + (t_first + tc) * N_HR * N_HR, t_offset + (t_first + tc + 1) * N_HR * N_HR - 1);
    }
}


/* **********************************************************************************/
/* ********************************** MAIN FUNCTION *********************************/
/* **********************************************************************************/

int main()
{
    /* Creat the file to save results */
    char *varnames[NUM_VARS] = {"x_rec_all"};
    create_netcdf(FILENAME_WR, NUM_VARS, varnames);

    /* Allocate memory */
    double *x_fusion_lf_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));
    double *x_fusion_hf_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));
    double *x_rec_all = (double*)malloc(NUM_3DSNAPS * NUM_2DSNAPS * N_HR * N_HR * sizeof(double));

    /* read all snapshots */
    size_t start_ids[4] = {0, 0, 0, 0};
    size_t count_ids[4] = {NUM_3DSNAPS, NUM_2DSNAPS, N_HR, N_HR };
    read_netcdf(FILENAME_RD, "Uinterp_all", start_ids, count_ids, x_fusion_lf_all);
    read_netcdf(FILENAME_RD, "Udiff_all", start_ids, count_ids, x_fusion_hf_all);

    double time_all_start = omp_get_wtime();

    double *x_current_lf = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_current_hf = (double*)malloc(N_HR * N_HR * sizeof(double));
    double *x_rec = (double*)malloc(N_HR * N_HR * sizeof(double));

    long int grid_size = N_HR * N_HR * NEIGHBOR_FULLSIZE * NEIGHBOR_FULLSIZE * SIM_FULLSIZE * SIM_FULLSIZE;
    int *gridpatches_y = (int*)malloc(grid_size * sizeof(int));
    int *gridpatches_z = (int*)malloc(grid_size * sizeof(int));
    int *acc_ids = (int*)malloc(ACC_FULLSIZE * ACC_FULLSIZE * sizeof(int));
    generate_grids(gridpatches_y, gridpatches_z, acc_ids);


    for(int snap3d_id = 0; snap3d_id < NUM_3DSNAPS; snap3d_id++)
    {
        int t_offset = snap3d_id * NUM_2DSNAPS * N_HR*N_HR;

        // put first PIV
        get_onesnap(x_fusion_hf_all, x_current_hf, t_offset + 0 * N_HR * N_HR, t_offset + 1 * N_HR * N_HR - 1);
        put_onesnap(x_rec_all, x_current_hf, t_offset + 0 * N_HR * N_HR, t_offset + 1 * N_HR * N_HR - 1);

        int block_id;
        for(block_id = 0; block_id < NUM_BLOCKS; block_id++)
        {
            double time_start = omp_get_wtime();

            int t_first = SCALE_FACTOR_TIME*block_id;
            int t_last = SCALE_FACTOR_TIME*(block_id+1);

            // Put last PIV of the block
            get_onesnap(x_fusion_hf_all, x_current_hf, t_offset + t_last * N_HR * N_HR, t_offset + (t_last + 1) * N_HR * N_HR - 1);
            put_onesnap(x_rec_all, x_current_hf, t_offset + t_last * N_HR * N_HR, t_offset + (t_last + 1) * N_HR * N_HR - 1);

            propag_towardcenter(x_rec_all, x_fusion_lf_all, gridpatches_y, gridpatches_z, acc_ids, t_first, t_offset);
            printf("\n Estimated block %i (total 23) in 3D snapshot %i (total 37) in %f seconds \n", block_id, snap3d_id, (double)omp_get_wtime() - time_start);
        }
    }

    // Write to file
    write_netcdf(FILENAME_WR, "x_rec_all", start_ids, count_ids, x_rec_all);

    /* free memory */
    free(x_rec); free(x_current_lf); free(x_current_hf);
    free(x_rec_all); free(x_fusion_lf_all); free(x_fusion_hf_all);
    free(gridpatches_y); free(gridpatches_z); free(acc_ids);
    printf("\n FINISH ALL COMPUTATION IN %f SECONDS \n", (double)omp_get_wtime() - time_all_start);

    return 1;
}
