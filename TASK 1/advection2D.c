/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - inital values of u(x,y)
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files
**********************************************************************/

#include <stdio.h>
#include <math.h>
#include <omp.h>

/*********************************************************************
                      Main function
**********************************************************************/

int main()
{

    /* Grid properties */
    const int NX = 1000;    // Number of x points
    const int NY = 1000;    // Number of y points
    const float xmin = 0.0; // Minimum x value
    const float xmax = 1.0; // Maximum x value
    const float ymin = 0.0; // Minimum y value
    const float ymax = 1.0; // Maximum y value

    /* Parameters for the Gaussian initial conditions */
    const float x0 = 0.1;                  // Centre(x)
    const float y0 = 0.1;                  // Centre(y)
    const float sigmax = 0.03;             // Width(x)
    const float sigmay = 0.03;             // Width(y)
    const float sigmax2 = sigmax * sigmax; // Width(x) squared
    const float sigmay2 = sigmay * sigmay; // Width(y) squared

    /* Boundary conditions */
    const float bval_left = 0.0;  // Left boudnary value
    const float bval_right = 0.0; // Right boundary value
    const float bval_lower = 0.0; // Lower boundary
    const float bval_upper = 0.0; // Upper bounary

    /* Time stepping parameters */
    const float CFL = 0.9;   // CFL number
    const int nsteps = 1500; // Number of time steps

    /* Velocity */
    const float velx = 0.01; // Velocity in x direction
    const float vely = 0.01; // Velocity in y direction

    /* Arrays to store variables. These have NX+2 elements
       to allow boundary values to be stored at both ends */
    float x[NX + 2];            // x-axis values
    float y[NX + 2];            // y-axis values
    float u[NX + 2][NY + 2];    // Array of u values
    float dudt[NX + 2][NY + 2]; // Rate of change of u

    float x2; // x squared (used to calculate iniital conditions)
    float y2; // y squared (used to calculate iniital conditions)

    /* Calculate distance between points */
    float dx = (xmax - xmin) / ((float)NX);
    float dy = (ymax - ymin) / ((float)NY);

    /* Calculate time step using the CFL condition */
    /* The fabs function gives the absolute value in case the velocity is -ve */
    float dt = CFL / ((fabs(velx) / dx) + (fabs(vely) / dy));

    /*** Report information about the calculation ***/
    printf("Grid spacing dx     = %g\n", dx);
    printf("Grid spacing dy     = %g\n", dy);
    printf("CFL number          = %g\n", CFL);
    printf("Time step           = %g\n", dt);
    printf("No. of time steps   = %d\n", nsteps);
    printf("End time            = %g\n", dt * (float)nsteps);
    printf("Distance advected x = %g\n", velx * dt * (float)nsteps);
    printf("Distance advected y = %g\n", vely * dt * (float)nsteps);

/*** Place x points in the middle of the cell ***/
// can be parallelised as calculation is not sequential, rather populating x[i] with calculations
/* LOOP 1 */
#pragma omp parallel for default(none) shared(x, dx, NX)
    for (int i = 0; i < NX + 2; i++)
    {
        x[i] = ((float)i - 0.5) * dx;
    }

    /*** Place y points in the middle of the cell ***/
    // can be parallelised too! same as before
    /* LOOP 2 */
#pragma omp parallel for default(none) shared(y, dy, NY)
    for (int j = 0; j < NY + 2; j++)
    {
        y[j] = ((float)j - 0.5) * dy;
    }

    /*** Set up Gaussian initial conditions ***/
    // can be parallelised as each population of u only utilises variables x2 or y2 which are calculated inside of the for loop and does not reuse u inside
    /* LOOP 3 */
#pragma omp parallel for collapse(2) default(none) shared(u, x, y, x0, y0, sigmax2, sigmay2, NY, NX) private(x2, y2)
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            x2 = (x[i] - x0) * (x[i] - x0);
            y2 = (y[j] - y0) * (y[j] - y0);
            u[i][j] = exp(-1.0 * ((x2 / (2.0 * sigmax2)) + (y2 / (2.0 * sigmay2))));
        }
    }

    /*** Write array of initial u values out to file ***/
    /* cannot be parallelised file order here must be sequential */
    FILE *initialfile;
    initialfile = fopen("initial.dat", "w");
    /* LOOP 4 */
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(initialfile);

    /*** Update solution by looping over time steps ***/
    /* cannot be parallelised as contents inside must be executed in order before LOOP 8 is executed*/
    /* LOOP 5 */
    {
        for (int m = 0; m < nsteps; m++)
        {

            /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
            /* can parellalised as no dependency of variable u between iterations*/
            /* LOOP 6 */
#pragma omp parallel for default(none) shared(u, bval_left, bval_right, NX, NY)
            for (int j = 0; j < NY + 2; j++)
            {
                u[0][j] = bval_left;
                u[NX + 1][j] = bval_right;
            }

            /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
            /* can be parellalised as no dependency of u between iterations */
            /* LOOP 7 */
#pragma omp parallel for default(none) shared(u, bval_lower, bval_upper, NX, NY)
            for (int i = 0; i < NX + 2; i++)
            {
                u[i][0] = bval_lower;
                u[i][NY + 1] = bval_upper;
            }

            /*** Calculate rate of change of u using leftward difference ***/
            /* Loop over points in the domain but not boundary values */
            /* can be parallelised as reading only from */
            /* LOOP 8 */
#pragma omp parallel for collapse(2) default(none) shared(u, dx, dy, NX, NY, velx, vely, dudt)
            for (int i = 1; i < NX + 1; i++)
            {
                for (int j = 1; j < NY + 1; j++)
                {
                    dudt[i][j] = -velx * (u[i][j] - u[i - 1][j]) / dx - vely * (u[i][j] - u[i][j - 1]) / dy;
                }
            }

/*** Update u from t to t+dt ***/
/* Loop over points in the domain but not boundary values */
/* Can be parallelised as only uses the same dependency as it is updating*/
/* LOOP 9 */
#pragma omp parallel for collapse(2) default(none) shared(u, dudt, dt, NX, NY)
            for (int i = 1; i < NX + 1; i++)
            {
                for (int j = 1; j < NY + 1; j++)
                {
                    u[i][j] = u[i][j] + dudt[i][j] * dt;
                }
            }

        } // time loop
    }
    /*** Write array of final u values out to file ***/
    /* cannot be parellised as operations inside the double for loop must be sequential */
    FILE *finalfile;
    finalfile = fopen("final.dat", "w");
    /* LOOP 10 */
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(finalfile);

    return 0;
}

/* End of file ******************************************************/