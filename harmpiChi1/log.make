bounds.c: In function ‘pack_prim’:
bounds.c:498:14: warning: ‘istop’ may be used uninitialized in this function [-Wmaybe-uninitialized]
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
              ^~~~~
bounds.c:498:7: warning: ‘istart’ may be used uninitialized in this function [-Wmaybe-uninitialized]
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
       ^~~~~~
In file included from bounds.c:47:0:
bounds.c: In function ‘pack_pflag’:
decs.h:340:1: warning: ‘jstop’ may be used uninitialized in this function [-Wmaybe-uninitialized]
 for(j=jstart;j<=jstop;j++)\
 ^~~
bounds.c:552:27: note: ‘jstop’ was declared here
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
                           ^~~~~
bounds.c:552:20: warning: ‘jstart’ may be used uninitialized in this function [-Wmaybe-uninitialized]
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
                    ^~~~~~
In file included from bounds.c:47:0:
decs.h:339:1: warning: ‘istop’ may be used uninitialized in this function [-Wmaybe-uninitialized]
 for(i=istart;i<=istop;i++)\
 ^~~
bounds.c:552:14: note: ‘istop’ was declared here
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
              ^~~~~
bounds.c:552:7: warning: ‘istart’ may be used uninitialized in this function [-Wmaybe-uninitialized]
   int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
       ^~~~~~
coord.c:661:6: warning: ‘vofx_cylindrified’ defined but not used [-Wunused-function]
 void vofx_cylindrified( double *Xin, void (*vofx)(double*, double*), double *Vout )
      ^~~~~~~~~~~~~~~~~
diag.c: In function ‘diag_flux’:
diag.c:326:1: warning: multi-line comment [-Wcomment]
 //#pragma omp parallel for schedule(static,MY_MAX(N2*N3/nthreads,1)) collapse(2) \
 ^
diag.c: In function ‘diag’:
diag.c:196:7: warning: ‘kmax’ may be used uninitialized in this function [-Wmaybe-uninitialized]
       fprintf(stderr,"LOG      t=%g \t divbmax: %d %d %d %g\n",
       ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             t,imax,jmax,kmax,divbmax) ;
             ~~~~~~~~~~~~~~~~~~~~~~~~~
diag.c:196:7: warning: ‘jmax’ may be used uninitialized in this function [-Wmaybe-uninitialized]
diag.c:196:7: warning: ‘imax’ may be used uninitialized in this function [-Wmaybe-uninitialized]
metric.c: In function ‘gdet_func’:
metric.c:117:55: warning: iteration 4 invokes undefined behavior [-Waggressive-loop-optimizations]
   for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
                                         ~~~~~~~~~~~~~~^~~~~~~~~~~~
metric.c:117:3: note: within this loop
   for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
   ^~~
metric.c: In function ‘get_geometry’:
metric.c:273:21: warning: iteration 4 invokes undefined behavior [-Waggressive-loop-optimizations]
    geom->gcon[0][j] = gcon[ii][jj][kk][loc][0][j];
    ~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metric.c:272:2: note: within this loop
  for(j=0;j<=NDIM*NDIM-1;j++){
  ^~~
lu.c: In function ‘invert_matrix’:
lu.c:69:53: warning: iteration 4 invokes undefined behavior [-Waggressive-loop-optimizations]
   for( i = 0 ; i < NDIM*NDIM ; i++ ) {  Amtmp[0][i] = Am[0][i]; }
                                         ~~~~~~~~~~~~^~~~~~~~~~
lu.c:69:3: note: within this loop
   for( i = 0 ; i < NDIM*NDIM ; i++ ) {  Amtmp[0][i] = Am[0][i]; }
   ^~~
utoprim_1dfix1.c:776:13: warning: ‘func_1d_orig2’ defined but not used [-Wunused-function]
 static void func_1d_orig2(FTYPE x[], FTYPE dx[], FTYPE resid[],
             ^~~~~~~~~~~~~
utoprim_1dfix1.c:435:14: warning: ‘dvsq_dW’ defined but not used [-Wunused-function]
 static FTYPE dvsq_dW(FTYPE W)
              ^~~~~~~
utoprim_1dvsq2fix1.c:87:14: warning: ‘vsq_calc’ declared ‘static’ but never defined [-Wunused-function]
 static FTYPE vsq_calc(FTYPE W);
              ^~~~~~~~
