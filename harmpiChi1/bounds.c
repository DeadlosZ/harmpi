//Modified by Alexander Tchekhovskoy: MPI+3D
/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/


#include "decs.h"

void bound_mpi_dim(int dim, int ispack, double prim[][N2M][N3M][NPR], int pflag[][N2M][N3M]);
int pack_prim(int ispack, int dim, int isup, double prim[][N2M][N3M][NPR], double *mpi_buf );
int pack_pflag(int ispack, int dim, int isup, int pflag[][N2M][N3M], int *mpi_buf );
void bound_x1dn(double prim[][N2M][N3M][NPR] );
void bound_x2dn(double prim[][N2M][N3M][NPR] );
void bound_x2dn_polefix(double prim[][N2M][N3M][NPR] );
void bound_x3dn(double prim[][N2M][N3M][NPR] );
void bound_x1up(double prim[][N2M][N3M][NPR] );
void bound_x2up(double prim[][N2M][N3M][NPR] );
void bound_x2up_polefix(double prim[][N2M][N3M][NPR] );
void bound_x3up(double prim[][N2M][N3M][NPR] );

// YOAV: wind boundry conditions functions
void bound_wind( double prim[][N2M][N3M][NPR] );
void bound_disk_noMPI( double prim[][N2M][N3M][NPR]);
void bound_disk_odd( double prim[][N2M][N3M][NPR]);
void bound_disk_even_bot( double prim[][N2M][N3M][NPR]);
void bound_disk_even_top( double prim[][N2M][N3M][NPR]);
void beta_profile( double prim[][N2M][N3M][NPR], int ii, int jj, int kk );
void rho_profile( double prim[][N2M][N3M][NPR], int ii, int jj, int kk );
void cart_vel_to_BL( double r, double theta, double phi, double vel[NDIM] );
void coord_transform( double *pr,int ii, int jj, int kk );

/* bound array containing entire set of primitive variables */

void bound_x1dn(double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  void inflow_check(double *pr, int ii, int jj, int kk, int type );
  struct of_geom geom ;
  
  int iNg, jNg, kNg;
  
#if(N1!=1)
  if (!is_physical_bc(1, 0)) return;

  /* inner r boundary condition */
  for(j=0;j<N2;j++)
  {
    for(k=0; k<N3;k++)  //!!!ATCH: make sure don't need to expand to transverse ghost cells
    {
#if( RESCALE )
      get_geometry(0,j,k,CENT,&geom) ;
      rescale(prim[0][j][k],FORWARD, 1, 0,j,k,CENT,&geom) ;
#endif
      
      for (iNg=-N1G; iNg<0; iNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[iNg][j][k][m] = prim[N1+iNg][j][k][m];
        pflag[iNg][j][k] = pflag[N1+iNg][j][k];
#elif(DIRICHLET==1)
        PLOOP prim[iNg][j][k][m] = pbound[iNg][j][k][m];
        pflag[iNg][j][k] = pflag[0][j][k] ;
#else //outflow
        PLOOP prim[iNg][j][k][m] = prim[0][j][k][m];
        pflag[iNg][j][k] = pflag[0][j][k] ;
#endif
      }
      
#if( RESCALE )
      for (iNg = -N1G; iNg<=0; iNg++)
      {
        get_geometry(iNg,j,k,CENT,&geom) ;
        rescale(prim[iNg][j][k],REVERSE, 1, iNg,j,k,CENT,&geom) ;
      }
#endif
    }
  }

  /* make sure there is no inflow at the inner boundary */
  if(1!=INFLOW) {
    for(i=-N1G;i<=-1;i++)  for(j=-N2G;j<N2+N2G;j++) for(k=-N3G;k<N3+N3G;k++)
    {
      //!!!ATCH: eliminated one loop that seemed excessive. verify.
      inflow_check(prim[i][j][k],i,j,k,0) ; //0 stands for -x1 boundary
    }
  }
#endif
  
}

void bound_x1up(double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  void inflow_check(double *pr, int ii, int jj, int kk, int type );
  struct of_geom geom ;
  
  int iNg, jNg, kNg;
#if(N1!=1)
  if (!is_physical_bc(1, 1)) return;

  /* Outer r boundary condition */
  for(j=0;j<N2;j++)
  {
    for(k=0; k<N3;k++)  //!!!ATCH: make sure don't need to expand to transverse ghost cells
    {
#if( RESCALE )
      get_geometry(N1-1,j,k,CENT,&geom) ;
      rescale(prim[N1-1][j][k],FORWARD, 1, N1-1,j,k,CENT,&geom) ;
#endif
      
      for (iNg=0; iNg<N1G; iNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[N1+iNg][j][k][m] = prim[iNg][j][k][m];
        pflag[N1+iNg][j][k] = pflag[iNg][j][k] ;
#elif(DIRICHLET==1)
        PLOOP prim[N1+iNg][j][k][m] = pbound[N1+iNg][j][k][m];
        pflag[N1+iNg][j][k] = pflag[N1-1][j][k] ;
#else //outflow
        PLOOP prim[N1+iNg][j][k][m] = prim[N1-1][j][k][m];
        pflag[N1+iNg][j][k] = pflag[N1-1][j][k] ;
#endif
      }
#if( RESCALE )
      for (iNg= 0; iNg<N1G+1; iNg++) //!!!ATCH: added +1 to N1G to ensure that all ghost cells are looped over
      {
        get_geometry(N1-1+iNg,j,k,CENT,&geom) ;
        rescale(prim[N1-1+iNg][j][k],REVERSE, 1, N1-1+iNG,j,k,CENT,&geom) ;
      }
#endif
    }
  }
  /* make sure there is no inflow at the outer boundary */
  if(1!=INFLOW) {
    for(i=N1;i<=N1+N1G-1;i++)  for(j=-N2G;j<N2+N2G;j++) for(k=-N3G;k<N3+N3G;k++)
    {
      //!!!ATCH: eliminated one loop that seemed excessive. verify.
      inflow_check(prim[i][j][k],i,j,k,1) ; //1 stands for +x1 boundary
    }
  }

#endif
}

void bound_x2dn_polefix( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;
  
#if(POLEFIX && POLEFIX < N2/2 && BL)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 0)) return;

  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  for(i=-N1G;i<N1+N1G;i++) {
    for(k=-N3G;k<N3+N3G;k++) {
      for(j=0;j<jref;j++) {
        PLOOP {
          if(m==B1 || m==B2 || (N3>1 && m==B3))
            //don't touch magnetic fields
            continue;
          else if(m==U2) {
            //linear interpolation of transverse velocity (both poles)
            prim[i][j][k][m] = (j+0.5)/(jref+0.5) * prim[i][jref][k][m];
          }
          else {
            //everything else copy (both poles)
            prim[i][j][k][m] = prim[i][jref][k][m];
          }
        }
      }
    }
  }
#endif
}

void bound_x2up_polefix( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(POLEFIX && POLEFIX < N2/2 && BL)
  //only do anything if physically in an MPI process that touches the outer pole
  if (!is_physical_bc(2, 1)) return;

  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  for(i=-N1G;i<N1+N1G;i++) {
    for(k=-N3G;k<N3+N3G;k++) {
      for(j=0;j<jref;j++) {
        PLOOP {
          if(m==B1 || m==B2 || (N3>1 && m==B3))
            //don't touch magnetic fields
            continue;
          else if(m==U2) {
            //linear interpolation of transverse velocity (both poles)
            prim[i][N2-1-j][k][m] = (j+0.5)/(jref+0.5) * prim[i][N2-1-jref][k][m];
          }
          else {
            //everything else copy (both poles)
            prim[i][N2-1-j][k][m] = prim[i][N2-1-jref][k][m];
          }
        }
      }
    }
  }
#endif
  
}

/* polar BCs */
//inner theta boundary
void bound_x2dn( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

  //YOAV: generate wind
#if (WIND==1)
  bound_wind( prim );
#endif

#if(N2!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 0)) return;

  bound_x2dn_polefix(prim);
  
  for (i=-N1G; i<N1+N1G; i++)
  {
    for (k=-N3G; k<N3+N3G; k++)
    {
      for (jNg=-N2G; jNg<0; jNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[i][jNg][k][m] = prim[i][N2+jNg][k][m];
        pflag[i][jNg][k] = pflag[i][N2+jNg][k];
#elif(DIRICHLET==1)
        PLOOP prim[i][jNg][k][m] = pbound[i][jNg][k][m];
        pflag[i][jNg][k] = pflag[i][-jNg-1][k];
#elif (OUTFLOW==1)
        PLOOP prim[i][jNg][k][m] = prim[i][0][k][m];
        pflag[i][jNg][k] = pflag[i][0][k];
#else //symmetric/asymmetric
        PLOOP prim[i][jNg][k][m] = prim[i][-jNg-1][k][m];
        pflag[i][jNg][k] = pflag[i][-jNg-1][k];
#endif
        
      }
    }
  }
  
  /* polar BCs */
  /* make sure b and u are antisymmetric at the poles */
  if(BL){
    for(i=-N1G;i<N1+N1G;i++) {
      for(k=-N3G;k<N3+N3G;k++) {
        for(j=-N2G;j<0;j++) {
          prim[i][j][k][U2] *= -1. ;
          prim[i][j][k][B2] *= -1. ;
        }
      }
    }
  }

#endif

}

/* polar BCs */
//outer theta boundary
void bound_x2up( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(N2!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 1)) return;

  bound_x2up_polefix(prim);

  for (i=-N1G; i<N1+N1G; i++)
  {
    for (k=-N3G; k<N3+N3G; k++)
    {
      //outer theta boundary
      for (jNg=0; jNg<N2G; jNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[i][N2+jNg][k][m] = prim[i][jNg][k][m];
        pflag[i][N2+jNg][k] = pflag[i][jNg][k];
#elif(DIRICHLET==1)
        PLOOP prim[i][N2+jNg][k][m] = pbound[i][N2+jNg][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-jNg-1][k];
#elif (OUTFLOW==1)
        PLOOP prim[i][N2+jNg][k][m] = prim[i][N2-1][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-1][k];
#else //symmetric/asymmetric
        PLOOP prim[i][N2+jNg][k][m] = prim[i][N2-jNg-1][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-jNg-1][k];
#endif
      }
    }
  }
  
  /* polar BCs */
  /* make sure b and u are antisymmetric at the poles */
  if(BL){
    for(i=-N1G;i<N1+N1G;i++) {
      for(k=-N3G;k<N3+N3G;k++) {
        for(j=N2;j<N2+N2G;j++) {
          prim[i][j][k][U2] *= -1. ;
          prim[i][j][k][B2] *= -1. ;
        }
      }
    }
  }

#endif
}

void bound_x3dn( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(N3!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(3, 0)) return;

  /* phi BCs */
  //inner phi-boundary
  for (i=-N1G; i<N1+N1G; i++)
  {
    for (j=-N2G; j<N2+N2G; j++)
    {
      for (kNg=-N3G; kNg<0; kNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[i][j][kNg][m] = prim[i][j][N3+kNg][m];
        pflag[i][j][kNg] = pflag[i][j][N3+kNg];
#elif(DIRICHLET==1)
        PLOOP prim[i][j][kNg][m] = pbound[i][j][kNg][m];
        pflag[i][j][kNg] = pflag[i][j][0];
#elif (OUTFLOW==1)
        PLOOP prim[i][j][kNg][m] = prim[i][j][0][m];
        pflag[i][j][kNg] = pflag[i][j][0];
#else
        //periodic by default
        PLOOP prim[i][j][kNg][m] = prim[i][j][N3+kNg][m];
        pflag[i][j][kNg] = pflag[i][j][N3+kNg];
#endif
      }
    }
  }
#endif

}

void bound_x3up( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(N3!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(3, 1)) return;

  /* phi BCs */
  //outer phi-boundary
  for (i=-N1G; i<N1+N1G; i++)
  {
    for (j=-N2G; j<N2+N2G; j++)
    {
      for (kNg=-N3G; kNg<0; kNg++)
      {
        for (kNg=0; kNg<N3G; kNg++)
        {
#if (PERIODIC==1)
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][kNg];
#elif(DIRICHLET==1)
          PLOOP prim[i][j][N3+kNg][m] = pbound[i][j][N3+kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][N3-1];
#elif (OUTFLOW==1)
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][N3-1][m];
          pflag[i][j][N3+kNg] = pflag[i][j][N3-1];
#else
          //periodic by default
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][kNg];
#endif
        }
      }
    }
  }
#endif
}

void bound_prim( double prim[][N2M][N3M][NPR] )
{
  int ispack, dim;
  
  //MPIMARK: could be optimized by individually passing corner zones:
  //         then, speed up by not doing comm dimension by dimension

  //x1-dim
  //packing, putting send and receive requests
  dim = 1;
  ispack = 1;
  bound_mpi_dim(dim, ispack, prim, pflag);
  //while waiting for MPI comm to complete, do physical boundaries
  bound_x1dn(prim);
  bound_x1up(prim);
  //waiting for comm to complete and unpacking
  ispack = 0;
  bound_mpi_dim(dim, ispack, prim, pflag);

  //x2-dim
  //packing, putting send and receive requests
  dim = 2;
  ispack = 1;
  bound_mpi_dim(dim, ispack, prim, pflag);
  //while waiting for MPI comm to complete, do physical boundaries
  bound_x2dn(prim);
  bound_x2up(prim);
  //waiting for comm to complete and unpacking
  ispack = 0;
  bound_mpi_dim(dim, ispack, prim, pflag);
  

  //x3-dim
  //waiting for comm to complete and unpacking
  dim = 3;
  ispack = 1;
  bound_mpi_dim(dim, ispack, prim, pflag);
  //while waiting for MPI comm to complete, do physical boundaries
  bound_x3dn(prim);
  bound_x3up(prim);
  //waiting for comm to complete and unpacking
  ispack = 0;
  bound_mpi_dim(dim, ispack, prim, pflag);
}

//packs (ispack=1) or unpacks (ispack=0) the cells to be communicated along
//dimension (dim), either upper (isup=1) or lower (isup=0) boundary
//returns the number of items packed (count)
int pack_prim(int ispack, int dim, int isup, double prim[][N2M][N3M][NPR], double *mpi_buf )
{
  int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;

  //if true, ensure it has value of unity for the below to work
  if(ispack) ispack=1;
  
  //do not do empty dimensions
  if(mpi_ntot[dim] == 1) return(0);
  
  //figure out the range of indices to transfer
  //x1: in x1-dim transfer only ghost cells immediately adjacent to active cells
  if(1==dim){
    jstart=0; jstop=N2-1;
    kstart=0; kstop=N3-1;
    if(isup) istart=N1-ispack*N1G; else istart=-N1G+ispack*N1G;
    istop=istart+N1G-1;
  }
  //x2: in x2-dim, additionally trasfer the ghost cells that have been communicated in x1-dim
  else if(2==dim){
    istart=-N1G; istop=N1+N1G-1;
    kstart=0; kstop=N3-1;
    if(isup) jstart=N2-ispack*N2G; else jstart=-N2G+ispack*N2G;
    jstop=jstart+N2G-1;
  }
  //x3, in x3-dim, additionally trasfer the ghost cells that have been communicated in x1-dim and x2-dim
  else if(3==dim){
    istart=-N1G; istop=N1+N1G-1;
    jstart=-N2G; jstop=N2+N2G-1;
    if(isup) kstart=N3-ispack*N3G; else kstart=-N3G+ispack*N3G;
    kstop=kstart+N3G-1;
  }
  
  //initialize the counter of the number of doubles (un)packed
  count = 0;
  
  //packing
  if(ispack){
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) PLOOP {
      mpi_buf[count++] = prim[i][j][k][m];
    }
  }
  ///unpacking
  else {
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) PLOOP {
      prim[i][j][k][m] = mpi_buf[count++];
    }
  }
  return(count);
}

//packs (ispack=1) or unpacks (ispack=0) the cells to be communicated along
//dimension (dim), either upper (isup=1) or lower (isup=0) boundary
//returns the number of items packed (count)
int pack_pflag(int ispack, int dim, int isup, int pflag[][N2M][N3M], int *mpi_buf )
{
  int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
  //number of ghost cells to copy: need only layer of thickness one for pflag
  int n1g = (N1G>0), n2g = (N2G>0), n3g = (N3G>0);
  
  //if true, ensure it has value of unity for the below to work
  if(ispack) ispack=1;
  
  //do not do empty dimensions
  if(mpi_ntot[dim] == 1) return(0);
  
  //figure out the range of indices to transfer
  //x1: in x1-dim transfer only ghost cells immediately adjacent to active cells
  if(1==dim){
    jstart=0; jstop=N2-1;
    kstart=0; kstop=N3-1;
    if(isup) istart=N1-ispack*n1g; else istart=-n1g+ispack*n1g;
    istop=istart+n1g-1;
  }
  //x2: in x2-dim, additionally trasfer the ghost cells that have been communicated in x1-dim
  else if(2==dim){
    istart=-n1g; istop=N1+n1g-1;
    kstart=0; kstop=N3-1;
    if(isup) jstart=N2-ispack*n2g; else jstart=-n2g+ispack*n2g;
    jstop=jstart+n2g-1;
  }
  //x3, in x3-dim, additionally trasfer the ghost cells that have been communicated in x1-dim and x2-dim
  else if(3==dim){
    istart=-n1g; istop=N1+n1g-1;
    jstart=-n2g; jstop=N2+n2g-1;
    if(isup) kstart=N3-ispack*n3g; else kstart=-n3g+ispack*n3g;
    kstop=kstart+n3g-1;
  }
  
  //initialize the counter of the number of doubles (un)packed
  count = 0;
  
  //packing
  if(ispack){
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      mpi_buf[count++] = pflag[i][j][k];
    }
  }
  ///unpacking
  else {
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      pflag[i][j][k] = mpi_buf[count++];
    }
  }
  return(count);
}


int is_physical_bc( int dim, int isup )
{
  
  //dimension is trivial => boundary is not physical
  if( 1 == mpi_ntot[dim] )
    return(0);
  
  //lower boundary is physical
  if( 0 == isup && 0 == mpi_coords[dim] && (1 == mpi_dims[dim] || 1 != mpi_periods[dim]) )
    return(1);
  
  //upper boundary is physical
  if( 1 == isup && mpi_dims[dim]-1 == mpi_coords[dim] && (1 == mpi_dims[dim] || 1 != mpi_periods[dim]) )
    return(1);
  
  //MPI boundary
  return(0);
}

//initiates send (isrecv = 0) or receive (isrecv = 1) operation in dimension dim (=0,1,2)
void bound_mpi_dim(int dim, int ispack, double prim[][N2M][N3M][NPR], int pflag[][N2M][N3M])
{
#ifdef MPI
#define DOMPIPFLAG (0)
    
  int count, count_pflag, tagsend, tagrecv, isup;
  
  //packing and sending
  if(1 == ispack) {
    for(isup=0;isup<=1;isup++) {
      //skip if on the physical boundary and the BC is not periodic
      if(is_physical_bc(dim,isup)) continue;

      //pack the data from prim[] into mpi send buffers
      count = pack_prim(ispack, dim, isup, prim, mpi_buf_send[dim][isup]);

#if(DOMPIPFLAG)
      //pack the data from pflag[] into mpi send buffers
      count_pflag = pack_pflag(ispack, dim, isup, pflag, mpi_buf_send_pflag[dim][isup]);
#endif

      //prims
      tagsend = 0+2*isup;
      tagrecv = 0+2*!isup;
      MPI_Isend(mpi_buf_send[dim][isup], //buffer that's being sent
                count,              //number of items sent
                MPI_DOUBLE,         //data type
                mpi_nbrs[dim][isup],//the rank of destination process
                tagsend,                //tag
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_send[dim][isup] //error
                );

      MPI_Irecv(mpi_buf_recv[dim][isup], //buffer that's being received
                count,              //number of items received (same as those sent)
                MPI_DOUBLE,         //data type
                mpi_nbrs[dim][isup],//the rank of source process
                tagrecv,                //tag (should be same as in the send process)
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_recv[dim][isup] //error
                );
#if(DOMPIPFLAG)
      //pflags
      tagsend = 1+2*isup;
      tagrecv = 1+2*!isup;
      MPI_Isend(mpi_buf_send_pflag[dim][isup], //buffer that's being sent
                count_pflag,              //number of items sent
                MPI_INT,            //data type
                mpi_nbrs[dim][isup],//the rank of destination process
                tagsend,                //tag
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_send_pflag[dim][isup] //error
                );
      
      MPI_Irecv(mpi_buf_recv_pflag[dim][isup], //buffer that's being received
                count_pflag,              //number of items received (same as those sent)
                MPI_INT,            //data type
                mpi_nbrs[dim][isup],//the rank of source process
                tagrecv,                //tag (should be same as in the send process)
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_recv_pflag[dim][isup] //error
                );
#endif
    }
  }
  //waiting, unpacking, and putting results back
  else {
    for(isup=0;isup<=1;isup++) {
      //skip if on the physical boundary and the BC is not periodic
      if(is_physical_bc(dim,isup)) continue;
      
      //wait for comminication to complete
      //prims
      MPI_Wait(&mpi_reqs_recv[dim][isup],&mpi_stat_recv[dim][isup]);
      MPI_Wait(&mpi_reqs_send[dim][isup],&mpi_stat_send[dim][isup]);
#if(DOMPIPFLAG)
      //pflags
      MPI_Wait(&mpi_reqs_recv_pflag[dim][isup],&mpi_stat_recv_pflag[dim][isup]);
      MPI_Wait(&mpi_reqs_send_pflag[dim][isup],&mpi_stat_send_pflag[dim][isup]);
#endif
      //unpack the data from mpi recv buffers into prim[]
      pack_prim(ispack, dim, isup, prim, mpi_buf_recv[dim][isup]);
#if(DOMPIPFLAG)
      //pack the data from mpi send buffers into  pflag[]
      pack_pflag(ispack, dim, isup, pflag, mpi_buf_recv_pflag[dim][isup]);
#endif
    }
  }
#endif
}

//do not allow "sucking" on the boundary:
//type = 0: do not allow ucon[1] > 0
//type = 1: do not allow ucon[1] < 0
//if a disallowed value detected, reset vcon[1] to zero
void inflow_check(double *pr, int ii, int jj, int kk, int type )
{
        struct of_geom geom ;
        double ucon[NDIM] ;
        int j,k ;
        double alpha,beta1,gamma,vsq ;

        get_geometry(ii,jj,kk,CENT,&geom) ;
        ucon_calc(pr, &geom, ucon) ;

        if( ((ucon[1] > 0.) && (type==0)) || ((ucon[1] < 0.) && (type==1)) ) { 
                /* find gamma and remove it from primitives */
	  if( gamma_calc(pr,&geom,&gamma) ) { 
	    fflush(stderr);
	    fprintf(stderr,"\ninflow_check(): gamma failure, (%d,%d,%d) \n",
                    ii+mpi_startn[1], jj+mpi_startn[2], kk+mpi_startn[3]);
	    fflush(stderr);
	    fail(FAIL_GAMMA);
	  }
	    pr[U1] /= gamma ;
	    pr[U2] /= gamma ;
	    pr[U3] /= gamma ;
	    alpha = 1./sqrt(-geom.gcon[0][0]) ;
	    beta1 = geom.gcon[0][1]*alpha*alpha ;

	    /* reset radial velocity so radial 4-velocity
	     * is zero */
	    pr[U1] = beta1/alpha ;

	    /* now find new gamma and put it back in */
	    vsq = 0. ;
	    SLOOP vsq += geom.gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
	    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
	    if( vsq >= 1. ) { 
	      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
	    }
	    gamma = 1./sqrt(1. - vsq) ;
	    pr[U1] *= gamma ;
	    pr[U2] *= gamma ;
	    pr[U3] *= gamma ;

	    /* done */
	  }
	  else
	    return ;

}

//YOAV: function that converts cartesian physical velocity to BL hpysical velocity
void cart_vel_to_BL(double x, double y, double z, double vel[NDIM]) {
  double w, dw_dt, r, theta, phi, v_x, v_y, v_z, v_r, v_t, v_p;

  v_x = vel[1];
  v_y = vel[2];
  v_z = vel[3];

  // w = pow(x,2) + pow(y,2) + pow(z,2) - pow(BL_A,2);

  // r = sqrt(0.5 * (w + sqrt(pow(w,2) + 4 * pow(BL_A,2) * pow(z,2))));
  // theta = acos(z / r);
  // phi = atan2(y, x);

  // dw_dt = 2 * (x * v_x + y * v_y + z * v_z);

  // v_r = (1 / (2 * r)) * (
  //                         (dw_dt / 2) + ((w * dw_dt + 4 * pow(BL_A,2) * z * v_z) /
  //                         (2 * sqrt(pow(w,2) + 4 * pow(BL_A,2) * pow(z,2))))
  //                       );
  // v_t = (-1 / sqrt(1 - pow(z / r,2))) * ((v_z * r - v_r * z) / (pow(r,2)));
  // v_p = (1 / (1 + pow(y / x,2))) * ((v_y * x - v_x * y) / (pow(x,2)));
  
// Assuming BL_A=0 and the its spherical
  v_r = (x*v_x+y*v_y+z*v_z)/sqrt(x*x+y*y+z*z);
  v_t = (v_x*y-v_y*x)/(x*x+y*y);
  v_p = (z*(x*v_x+y*v_y)-(x*x+y*y)*v_z)/((x*x+y*y+z*z)*sqrt(x*x+y*y));

  vel[1] = v_r;
  vel[2] = v_t;
  vel[3] = v_p;
}

//YOAV: disk creation for odd number of processors
void bound_disk_odd( double prim[][N2M][N3M][NPR]) {
  int i,k;
  double r;
  for (i=0; i<N1; i++)
  {
    for (k=0; k<N3; k++)
    { 
      // check if we are inside the disk
      get_phys_coord_r(i,N2/2-1,k,&r);
      if ((r<D_ROUT)&&(r>D_RIN)) {
        // we are inside the disk. inject wind
        rho_profile(prim,i,N2/2-1,k);
        rho_profile(prim,i,N2/2  ,k);

        beta_profile(prim,i,N2/2-1,k);
        beta_profile(prim,i,N2/2  ,k);
      }
    }
  }
}


//YOAV: disk creation for No MPI case (one processor)
void bound_disk_noMPI( double prim[][N2M][N3M][NPR]) {
  int i,k;
  double r;
  for (i=0; i<N1; i++)
  {
    for (k=0; k<N3; k++)
    { 
      get_phys_coord_r(i,N2/2-1,k,&r);
      if ((r<D_ROUT)&&(r>D_RIN)) {
        // we are inside the disk. inject wind
        rho_profile(prim,i,N2/2-1,k);
        rho_profile(prim,i,N2/2  ,k);

        beta_profile(prim,i,N2/2-1,k);
        beta_profile(prim,i,N2/2  ,k);
      }
    }
  }
}

//YOAV: disk creation of the bottom half of the processors
void bound_disk_even_bot( double prim[][N2M][N3M][NPR]) {
  int i,k;
  double r;
  for (i=0; i<N1; i++)
  {
    for (k=0; k<N3; k++)
    {
      // check if we are inside the disk
      get_phys_coord_r(i,N2-1,k,&r);
      if ((r<D_ROUT)&&(r>D_RIN)) {
        // we are inside the disk. inject wind
        rho_profile(prim,i,N2-1,k);
        beta_profile(prim,i,N2-1,k);
      }
    }
  }
}

//YOAV: disk creation of the top half of the processors
void bound_disk_even_top( double prim[][N2M][N3M][NPR]) {
  int i,k;
  double r;
  for (i=0; i<N1; i++)
  {
    for (k=0; k<N3; k++)
    {
      // check if we are inside the disk
      get_phys_coord_r(i,0,k,&r);
      if ((r<D_ROUT)&&(r>D_RIN)) {
        // we are inside the disk. inject wind
        rho_profile(prim,i,0,k);
        beta_profile(prim,i,0,k);
      }
    }
  }
}

//YOAV: call for wind boundry
void bound_wind( double prim[][N2M][N3M][NPR] ) {
  /*
  The whole function does:
  (1) Determine whether the number of proc is odd or even
  (2) If it's odd, create disk in the middle processor alone
  (3) If it's even, create disk in the connetion between the 2 middle processors
  (4) If not using multiprocessors, do the same as for mpi with odd number
  */
  int i, j, k;

#ifdef MPI
  if (mpi_dims[2]%2==1) {
    if (mpi_dims[2]/2 == mpi_coords[2]) {
      bound_disk_odd(prim);
    }
  } else {
      if (mpi_dims[2]/2-1 == mpi_coords[2]) { 
        bound_disk_even_bot(prim);
      } else if (mpi_dims[2]/2 == mpi_coords[2]) {
        bound_disk_even_top(prim);
      }
    }
#else
  bound_disk_noMPI(prim);
#endif
}

// YOAV: New version of the beta profile calculation
void beta_profile( double prim[][N2M][N3M][NPR], int ii, int jj, int kk) {

  double r, th, phi, x, y, z, vx, vy, vz;
  double z0, u_w0, theta_w, gamma, r_0;
  double A, B, C;
  double gam = 1.66;
  double beta[NDIM] = { 0 };
  double GridLoc = 0;
  struct of_geom geom;
  int i,j,k;
  double X[NDIM] = { 0 };

  // get metric and phys coords:
  coord(ii,jj,kk,CENT,X);
  bl_coord(X,&r,&th,&phi);
  blgset(ii,jj,kk,&geom);

  GridLoc = (M_PI_2-th)/fabs(M_PI_2-th); //-1 if bottom and 1 if top

  // fix for equatorial symmetry:
  th  = M_PI_2 - fabs(M_PI_2 - th);

  //calculate cartesian coords
  x = sqrt(r*r+BL_A*BL_A) * sin(th) * cos(phi);
  y = sqrt(r*r+BL_A*BL_A) * sin(th) * sin(phi);
  z = r * cos(th);
  
  //wind charactristics
  theta_w = -atan(x/(FOCAL_Z+z));
  r_0 = x/(1+z/FOCAL_Z);

  // calculate velocities
  A = 0.5;
  B = gam/(gam-1)*2*M_PI*fabs(r_0)*cos(theta_w)*BERN_W*log(D_ROUT/D_RIN)/LUMN_W*PRES_W;
  C = -(BERN_W+1)/fabs(r_0);
  // u_w0 = (-B+sqrt(B*B-4*A*C))/(2*A); // complete
  u_w0 = sqrt(-C/A); //p=0
  // u_w0 = 0.2/sqrt(fabs(r_0)); //for dipole!!

  // generate inject wind in cart coords
  vx = u_w0*sin(theta_w);
  vz = u_w0*cos(theta_w);

  // insert to beta array
  beta[1] = vx;
  beta[3] = vz;

  // fixing the sign of U2 as promised:
  beta[3] *= -GridLoc;

  // insert the new 4 velocity into prim array
  prim[ii][jj][kk][U1] = beta[1];
  prim[ii][jj][kk][U2] = beta[3];

  prim[ii][jj][kk][U1] *= sqrt(geom.gcon[1][1]);
  prim[ii][jj][kk][U2] *= sqrt(geom.gcon[2][2]);
  prim[ii][jj][kk][U3] *= sqrt(geom.gcon[3][3]);

  // and finally performing coord transform as seen at init_torus.c:
  // -> "convert from 4-vel in BL coords to relative 4-vel in code coords"
  coord_transform(prim[ii][jj][kk],ii,jj,kk);
  // fprintf(stderr, "%lf, %lf, %lf\n\n", prim[ii][jj][kk][U1], prim[ii][jj][kk][U2], prim[ii][jj][kk][U3]);
}

//YOAV: calculation of rho profile
void rho_profile( double prim[][N2M][N3M][NPR], int ii, int jj, int kk ) {
  double r, th, phi, x, y, z;
  double u_w0, theta_w, r_0;
  double A,B,C;
  double gam = 1.66;
  struct of_geom geom;
  int i,j,k;
  double X[NDIM] = { 0 };
  double rho_w;

  // get metric and phys coords:
  coord(ii,jj,kk,CENT,X);
  bl_coord(X,&r,&th,&phi);
  blgset(ii,jj,kk,&geom);

  // fix for equatorial symmetry:
  th  = M_PI_2 - fabs(M_PI_2 - th);

  //calculate cartesian coords
  x = sqrt(r*r+BL_A*BL_A) * sin(th) * cos(phi);
  y = sqrt(r*r+BL_A*BL_A) * sin(th) * sin(phi);
  z = r * cos(th);

  theta_w = -atan(x/(FOCAL_Z+z));
  r_0 = x/(1+z/FOCAL_Z);

  // calculate velocities
  A = 0.5;
  B = gam/(gam-1)*2*M_PI*fabs(r_0)*cos(theta_w)*BERN_W*log(D_ROUT/D_RIN)/LUMN_W*PRES_W;
  C = -(BERN_W+1)/fabs(r_0);
  // u_w0 = (-B+sqrt(B*B-4*A*C))/(2*A); // complete
  u_w0 = sqrt(-C/A); //p=0
  // u_w0 = 0.2/sqrt(fabs(r_0)); //for dipole!!

  // calculate rho profile
  rho_w = LUMN_W/(2*M_PI*u_w0*cos(theta_w)*BERN_W*log(D_ROUT/D_RIN)*fabs(r_0));

  // fprintf(stderr, "%1.15lf ,%1.15lf, %1.15lf\n", r, rho_w, u_w0);

  prim[ii][jj][kk][RHO] = rho_w;
  // prim[ii][jj][kk][UU]  = UUMINLIMIT;
}