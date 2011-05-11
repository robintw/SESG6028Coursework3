
/* 

Code that solves Laplace's equation in 3D for Dirichlet boundary conditions set in
grid_set_boundary.

Note this only work for NON-periodic boundary conditions

Written by I.J.Bush March 2011 - if it breaks your computer or aught else, tough.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grid.h"
#include "array_alloc.h"
#include "timer.h"
#include <mpi.h>

int grid_init( int ng[ 3 ], struct grid *g )
{

  /* Initialise the grid 
     Arguments are:
     ng        - The total number of points in the grid
     grid      - The derived type holding all the data about the grid

  */

  int i;
  int dim_size[3];
  int periods[3];
  int coords[3];
  int npx = 3;
  int npy = 2;
  int npz = 1;
  int rank;
  
  MPI_Comm cart_comm; 
  /* Calculate the dimensions of the process grid to use */
  
    /* Store the size of the full grid */
  for( i = 0; i < 3; i++ ){
    g->ng[ i ] = ng[ i ];
  }
  
  g->npx = npx;
  g->npy = npy;
  g->npz = npz;
  
  /* Create the cartesian communicator with the sizes of the processor grid
  Set all dimensions to be NON-PERIODIC as we don't need/want to loop back across the grid
  as we have to not update the dimensions anyway! */
  dim_size[0] = g->npz;
  dim_size[1] = g->npy;
  dim_size[2] = g->npx;
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Create the communicator */	
  MPI_Cart_create(MPI_COMM_WORLD, 3, dim_size, periods, 1, &cart_comm);
	
  /* Get our co-ordinates within that communicator */
  MPI_Cart_coords(cart_comm, rank, 3, coords);
  
  /* Who are my neighbours in each direction? */
  MPI_Cart_shift( cart_comm, 2, 1, &g->north, &g->south    );
  MPI_Cart_shift( cart_comm, 1, 1, &g->west, &g->east );
  MPI_Cart_shift( cart_comm, 0, 1, &g->up, &g->down );

  g->nx = ceil(ng[0] / (float) npx);
  g->ny = ceil(ng[1] / (float) npy);
  g->nz = ceil(ng[2] / (float) npz);
		
  if (coords[0] == (npz - 1))
  {
		/* We're at the far end of z */
		g->nz = ng[0] - (g->nz * (npz - 1));
  }
  if (coords[1] == (npy - 1))
  {
		/* We're at the far end of y */
		g->ny = ng[1] - (g->ny * (npy - 1));
  }
  if (coords[2] == (npx - 1))
  {
		/* We're at the far end of x */
		g->nx = ng[2] - (g->nx * (npx - 1));
  }
  
  /* Store the coords in the g structure */
  g->px = coords[0];
  g->py = coords[1];
  g->pz = coords[2];
  
  printf("Rank: %d at location (%d, %d, %d) with sizes x = %d, y = %d, z = %d. N: %d, S: %d, E: %d, W: %d, U: %d, D: %d\n", rank, coords[0], coords[1], coords[2], g->nx, g->ny, g->nz, g->north, g->south, g->east, g->west, g->up, g->down);

  /* Allocate the grid. Note we will need two versions of the data on the grid, one to hold
     the current values, and one to write the results into when we are updating the grid. We swap between
     the two versions as we go from iteration to iteration, the 'current' member of the derived type
     indicating which version is the most up to date.
      
     We're making the array 2 cells bigger in each dimension so that we have a spare cell at the end of every row etc. so that we
     can send and receive into/from that extra bit. */

 for( i = 0; i < 2; i++ ){

    g->data[ i ] = alloc_3d_double(g->nz + 2, g->ny + 2, g->nx +2 ); 
    if( g->data[ i ] == NULL )
    {
      printf("Failed to allocated memory.\n");
      return EXIT_FAILURE;
    }
  }
  
  /* Which version of the grid is the "current" version. The other we will write
     the next result into */
  g->current = 0;

  /* Zero the timer and the iteration counter*/
  g->t_iter = 0.0;
  g->n_iter = 0;

  printf("Returning\n");
  return EXIT_SUCCESS;

}

void grid_finalize( struct grid *g ) 
{

  /* Free the grid */

  free_3d_double( g->data[ 1 ], g->nx );
  free_3d_double( g->data[ 0 ], g->nx );

}

void grid_initial_guess( struct grid *g )
{

  /* Initial guess at the solution. We'll use a very simple one ... */

  int i, j, k, l;

  printf("Inside grid_initial_guess\n");
  
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < g->nx - 1; j++ ){
	for( k = 0; k < g->ny - 1; k++ ){
	  for( l = 0; l < g->nz - 1; l++ ){
	    g->data[ i ][ j ][ k ][ l ] = 0.0;
	  }
	}
      } 
    }

}

void grid_set_boundary( struct grid *g )
{
  /* Set the boundary conditions for the grid. We're using Dirichlet
     boundary conditions, so we need to set every point on the edge of the grid.

     As a simple example set all the faces to zero except the bottom xy face, which we set 
     to unity */

  int i, j, k;
  
  printf("Inside grid_set_boundary: px = %d / %d, py = %d / %d, pz = %d / %d\n", g->px, g->npx, g->py, g->npy, g->pz, g->npz);

  /* Set each face of the cuboid in turn */
  /* Also remember that we need to do do it for both versions */
  
  /* Our coords in the process grid are stored in g as px, py and pz. If any of these = 0 or = npx
  (or npy or npz) then we're on the edge of the cuboid, and where we are. We then need to choose
  the correct bit of the array to set. Remember we're doing this for both versions of the grid 
  
  We're also setting the upper and lower bounds for later loops (ub and lb variables). These should
  be 1 and (n - 2) in the normal case (because the array is 1 bigger in each dimension). We initialise
  them all first to something that will be correct for anything that isn't changed below. */
  
  g->lb_x = 1;
  g->ub_x = (g->nx - 2);
  
  g->lb_y = 1;
  g->ub_y = (g->ny - 2);
  
  g->lb_z = 1;
  g->ub_z = (g->nz - 2);
  
  if (g->px == 0)
  {
  	printf("At x = 0\n");
  	/* We're at the x = 0 face of the grid. So we need to set the entire face equal to zero */
  	for (k = 0; k < 2; k++)
  	{
  		for (i = 0; i < g->nz; i++)
  		{
  			for (j = 0; j < g->ny; j++)
  			{
  				g->data[k][i][j][1] = 0.0;
  			}
  		}
  	}
  	g->lb_x = 2;
  }
  if (g->px == (g->npx - 1))
  {
    printf("At x = npx\n");
	/* We're at the high x face of the grid. So we need to set the entire face equal to zero */
	for (k = 0; k < 2; k++)
  	{
  		for (i = 0; i < g->ny; i++)
  		{
  			for (j = 0; j < g->nz; j++)
  			{
  				g->data[k][j][i][g->nx-2] = 0.0;
  			}
  		}
  	}
  	g->ub_x = (g->nx - 3);
  }
  	
  	
  	if (g->py == 0)
  	{
  		printf("At y = 0\n");
		/* We're at the y = 0 face of the grid. So we need to set the entire face equal to zero */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < g->nx; i++)
			{
				for (j = 0; j < g->nz; j++)
				{
					g->data[k][j][1][i] = 0.0;
				}
			}
		}
		g->lb_y = 2;
	}
  	if (g->py == (g->npy - 1))
  	{
  		printf("At y = npy\n");
		/* We're at the high y face of the grid. So we need to set the entire face equal to zero */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < g->nx; i++)
			{
				for (j = 0; j < g->nz; j++)
				{
					g->data[k][j][g->ny-2][i] = 0.0;
				}
			}
		}
		g->ub_y = (g->ny - 3);
  	}
  	
  	if (g->pz == 0)
  	{
  		printf("At z = 0\n");
		/* We're at the z = 0 face of the grid. So we need to set the entire face equal to ONE (unity) */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < g->nx; i++)
			{
				for (j = 0; j < g->ny; j++)
				{
					g->data[k][1][j][i] = 1.0;
				}
			}
		}
		g->lb_z = 2;
	}
  	if (g->pz == (g->npz - 1))
  	{
  		printf("At z = npz\n");
		/* We're at the high z face of the grid. So we need to set the entire face equal to zero */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < g->nx; i++)
			{
				for (j = 0; j < g->ny; j++)
				{
					g->data[k][g->nz-2][j][i] = 0.0;
				}
			}
		}
		g->ub_z = (g->nz-2);
	}
	
	printf("Finished grid_set_boundary\n");
}

double grid_update( struct grid *g ) {

  /* Perform the grid update */

  #define ONE_SIXTH 0.166666666666666666666666666666666666666666666

  double dg, diff;
  double start, finish;

  int current, update;
  int i, j, k;
  
  int tag, rank;
  
  MPI_Request send_req, recv_req;
  MPI_Request requests[100];
  int req_num = 0;
  
  MPI_Datatype face;

  /* Work out which version of the grid hold the current values, and
     which we will write the update into */
  switch( g->current ) {
  case 0:  
    current = 0;
    update  = 1; 
    break;  
  case 1:  
    current = 1;
    update  = 0; 
    break;
  default: 
    fprintf( stderr, "Internal error: impossible value for g->current\n" );
    exit( EXIT_FAILURE );
  }
  

  /* Perform the update and check for convergence  */
  
  dg = 0.0;
  
  /* Get the neighbouring faces from the grid neighbours */
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (g->west >= 0)
  {
  	tag = (g->west + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }
  
  /* Exchange for the western face of the chunk */
  printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny);
  MPI_Isend(&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag, MPI_COMM_WORLD, &send_req);
  printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny);
  MPI_Irecv(&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag, MPI_COMM_WORLD, &recv_req);
  /* MPI_Sendrecv(&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag,
  		&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  if (g->east >= 0)
  {
  	tag = (g->east + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }		
  
  printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny);
  MPI_Isend(&g->data[current][g->nz-1][0][0], g->nx*g->ny, MPI_DOUBLE, g->east, tag, MPI_COMM_WORLD, &send_req);
  printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny);
  MPI_Irecv(&g->data[current][g->nz-1][0][0], g->nx*g->ny, MPI_DOUBLE, g->east, tag, MPI_COMM_WORLD, &recv_req);
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  printf("Waiting...\n");
  MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
  
  
  
  
  
  
  /* Exchange North/South faces */
  req_num = 0;
  MPI_Type_vector(2*g->ny, 1, g->nx, MPI_DOUBLE, &face);
  MPI_Type_commit(&face);
  
  if (g->north >= 0)
  {
  	tag = (g->north + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }
  
  /* Exchange for the northern face of the chunk */
  /* printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny); */
  MPI_Isend(&g->data[current][0][0][0], 1, face, g->north, tag, MPI_COMM_WORLD, &send_req);
  /* printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny); */
  MPI_Irecv(&g->data[current][0][0][0], 1, face, g->north, tag, MPI_COMM_WORLD, &recv_req);
  /* MPI_Sendrecv(&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag,
  		&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  if (g->south >= 0)
  {
  	tag = (g->south + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }		
  
  /* printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny); */
  MPI_Isend(&g->data[current][0][0][g->nx-1], 1, face, g->south, tag, MPI_COMM_WORLD, &send_req);
  /* printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny); */
  MPI_Irecv(&g->data[current][0][0][g->nx-1], 1, face, g->south, tag, MPI_COMM_WORLD, &recv_req);
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  printf("Waiting...\n");
  MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
  
    /* Exchange Up/Down faces */
  req_num = 0;
  MPI_Type_vector(g->nz, g->nx, g->nx * g->ny, MPI_DOUBLE, &face);
  MPI_Type_commit(&face);
  
  if (g->up >= 0)
  {
  	tag = (g->up + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }
  
  /* Exchange for the up  face of the chunk */
  /* printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny); */
  MPI_Isend(&g->data[current][0][0][0], 1, face, g->up, tag, MPI_COMM_WORLD, &send_req);
  /* printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->west, tag, g->nx*g->ny); */
  MPI_Irecv(&g->data[current][0][0][0], 1, face, g->up, tag, MPI_COMM_WORLD, &recv_req);
  /* MPI_Sendrecv(&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag,
  		&g->data[current][0][0][0], g->nx*g->ny, MPI_DOUBLE, g->west, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  if (g->down >= 0)
  {
  	tag = (g->down + 1) * (rank + 1);
  }
  else
  {
  	tag = 0;
  }		
  
  /* printf("Send: From %d to %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny); */
  MPI_Isend(&g->data[current][0][g->ny-1][0], 1, face, g->down, tag, MPI_COMM_WORLD, &send_req);
  /* printf("Recv: At %d from %d. Tag = %d. Size = %d\n", rank, g->east, tag, g->nx*g->ny); */
  MPI_Irecv(&g->data[current][0][g->ny-1][0], 1, face, g->down, tag, MPI_COMM_WORLD, &recv_req);
  
  requests[req_num] = send_req;
  requests[req_num+1] = recv_req;
  
  req_num += 2;
  
  printf("Waiting...\n");
  MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
  
  
  
  
  /* Loop through and do the calculations */
  for( i = g->lb_x; i <= g->ub_x; i++ )
  {
    for( j = g->lb_y; j <= g->ub_y; j++ )
    {
      for( k = g->lb_z; k <= g->ub_z; k++ )
      {
		g->data[ update ][ i ][ j ][ k ] = ONE_SIXTH * ( g->data[ current ][ i + 1 ][ j     ][ k     ] +
		g->data[ current ][ i - 1 ][ j     ][ k     ] +
		g->data[ current ][ i     ][ j + 1 ][ k     ] +
		g->data[ current ][ i     ][ j - 1 ][ k     ] +
		g->data[ current ][ i     ][ j     ][ k + 1 ] +
		g->data[ current ][ i     ][ j     ][ k - 1 ] );
		diff = fabs( g->data[ update ][ i ][ j ][ k ] - g->data[ current ][ i ][ j ][ k ] );
		dg = dg > diff ? dg : diff;
      }
    }
  }

  /* Update the iteration counter */
  g->n_iter++;

  /* The updated grid is now the current grid, so swap over */
  g->current = update;

  return dg;
}

double grid_checksum( struct grid g ){

  /* Simple checksum over the grid to help checking solution
     i..e simply add up all the grid points */

  double sum;

  int i, j, k;

  sum = 0.0;
  for( i = 1; i < g.nz - 1; i++ ) {
    for( j = 1; j < g.ny - 1; j++ ) {
      for( k = 1; k < g.nx - 1; k++ ) {
	sum += g.data[ g.current ][ i ][ j ][ k ];
      }
    }
  }

  return sum;

}

void grid_print_times( struct grid g ){

  /* Report the measured timing data */

  fprintf( stdout, "Timing breakdown for the grid operations:\n" );
  fprintf( stdout, "                                  Total     Per iteration\n" );
  fprintf( stdout, "Iteration time        :      %12.4f   %12.8f\n", g.t_iter, g.t_iter / g.n_iter );
  fprintf( stdout, "\n" );

}
