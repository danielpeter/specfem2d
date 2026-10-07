/*
!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================
*/


#ifdef USE_TEXTURES_FIELDS
realw_texture d_displ_tex;
realw_texture d_accel_tex;
// backward/reconstructed
realw_texture d_b_displ_tex;
realw_texture d_b_accel_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel(int x);


// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ<1>(int x) { return tex1Dfetch(d_displ_tex, x); }
template<> __device__ float texfetch_accel<1>(int x) { return tex1Dfetch(d_accel_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ<3>(int x) { return tex1Dfetch(d_b_displ_tex, x); }
template<> __device__ float texfetch_accel<3>(int x) { return tex1Dfetch(d_b_accel_tex, x); }
#endif

#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
__constant__ size_t d_hprime_xx_tex_offset;
realw_texture d_wxgll_xx_tex;
__constant__ size_t d_wxgll_xx_tex_offset;
#endif


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */


// loads displacement into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ(const int* tx, const int* iglob,
                                                          realw_const_p d_displ,
                                                          realw* sh_tempx,
                                                          realw* sh_tempz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
#ifdef USE_TEXTURES_FIELDS
  sh_tempx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*2);
  sh_tempz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*2 + 1);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_tempx[(*tx)] = d_displ[(*iglob)*2];
  sh_tempz[(*tx)] = d_displ[(*iglob)*2 + 1];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprime(const int* tx,
                                                           realw_const_p d_hprime_xx,
                                                           realw* sh_hprime_xx){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_hprime_xx[(*tx)] = tex1Dfetch(d_hprime_xx_tex,(*tx) + d_hprime_xx_tex_offset);
#else
  // hprime
  sh_hprime_xx[(*tx)] = d_hprime_xx[(*tx)];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_wxgll(const int* tx,
                                                           realw_const_p d_wxgll,
                                                           realw* sh_wxgll){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_wxgll[(*tx)] = tex1Dfetch(d_wxgll_xx_tex,(*tx) + d_wxgll_xx_tex_offset);
#else
  // hprime
  sh_wxgll[(*tx)] = d_wxgll[(*tx)];
#endif
}




/* ----------------------------------------------------------------------------------------------- */

// loads hprimewgll into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprimewgll(const int* tx,
                                                               realw_const_p d_hprimewgll_xx,
                                                               realw* sh_hprimewgll_xx) {

  // each thread reads its corresponding value
  // weighted hprime
//#ifdef USE_TEXTURES_CONSTANTS
  // hprime
//  sh_hprimewgll_xx[(*tx)] = tex1Dfetch(d_hprimewgll_xx_tex,(*tx));
//#else
  sh_hprimewgll_xx[(*tx)] = d_hprimewgll_xx[(*tx)];
//#endif
}

/* ----------------------------------------------------------------------------------------------- */



__device__  __forceinline__ void sum_hprime_xi(int I, int J,
                                              realw* tempxl,realw* tempzl,
                                              realw* sh_tempx,realw* sh_tempz, realw* sh_hprime) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+I];

    sumx += sh_tempx[J*NGLLX+l] * fac;
    sumz += sh_tempz[J*NGLLX+l] * fac;
  }

// counts:
// + NGLLX * ( 2 + 3*6 ) FLOP = 100 FLOP
//
// + 0 BYTE

  *tempxl = sumx;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */


__device__  __forceinline__ void sum_hprime_gamma(int I, int J,
                                                 realw* tempxl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprime) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+J];

    sumx += sh_tempx[l*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */



__device__  __forceinline__ void sum_hprimewgll_xi(int I, int J,
                                                   realw* tempxl,realw* tempzl,
                                                   realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];

    sumx += sh_tempx[J*NGLLX+l] * fac;
    sumz += sh_tempz[J*NGLLX+l] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}


/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_gamma(int I, int J,
                                                 realw* tempxl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[J*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];

    sumx += sh_tempx[l*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}


/* ----------------------------------------------------------------------------------------------- */
// computes dot product between tensor and derivatives

__device__  __forceinline__ void get_dot_product_wgll(realw jacobianl,
                                                      realw sigma_xx,realw sigma_xz,
                                                      realw sigma_zx,realw sigma_zz,
                                                      realw Dxl,realw Dzl,
                                                      realw* sh_tempx,realw* sh_tempz,
                                                      realw* sh_wxgll,
                                                      int idx,
                                                      int tx,
                                                      int p_sv){

  // fills shared memory arrays
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = sh_wxgll[idx] * jacobianl * (sigma_xx * Dxl + sigma_zx * Dzl); // sh_tempx1
      sh_tempz[tx] = sh_wxgll[idx] * jacobianl * (sigma_xz * Dxl + sigma_zz * Dzl); // sh_tempz1
    }else{
      // SH-case
      sh_tempx[tx] = sh_wxgll[idx] * jacobianl * (sigma_xx * Dxl + sigma_xz * Dzl); // sh_tempx1
      sh_tempz[tx] = 0.f;
    }
  }
  __syncthreads();
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for elastic domains

/* ----------------------------------------------------------------------------------------------- */

// note:
// kernel_2 is split into 2 kernels:
//  - a kernel without attenuation and for isotropic media: Kernel_2_noatt_iso_impl()
//  - a kernel without attenuation and for anisotropic media: Kernel_2_noatt_ani_impl()
//
// this should help with performance:
// the high number of registers needed for our kernels limits the occupancy; separation tries to reduce this.


// kernel without attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiz,
                        realw* d_gammax,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p wxgll,
                        realw* d_kappav,
                        realw* d_muv,
                        const int simulation_type,
                        const int p_sv,
                        const int PML,
                        const int* d_spec_to_pml,
                        realw* d_rhostore,
                        const int gravity,
                        const realw gravity_const){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .true. or .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;
  realw sum_terms1,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

// counts:
// 2 FLOP

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  //checks if element is outside the PML
  if (PML){
    if (d_spec_to_pml[working_element] > 0) return;
  }

// counts:
// + 4 FLOP
//
// + 1 float * 128 threads = 512 BYTE

  // limits thread ids to range [0,25-1]
  if (tx >= NGLL2 ) tx = tx - NGLL2 ;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

// counts:
// + 3 FLOP
//
// + 1 float * 128 threads = 512 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// + 2 * 1 float * 128 threads = 1024 BYTE

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load

//  xixl = d_xix[offset]; // first array with texture load
//  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
//  xizl = d_xiz[offset];

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

// counts:
// + 15 FLOP
//
// + 9 float * 128 threads = 4608 BYTE

  // local index
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

 // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  // note:
  // here, kappal and mul are taken from arrays kappastore and mustore,
  // while the CPU-routine takes values lambda and mu from poroelastcoef array
  //
  // conversion from kappa/mu to lambda/mu
  // AXISYM    : kappal = lambdal + TWO_THIRDS * mul
  // non-AXISYM: kappal = lambdal + mul

  // original
  //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  //lambdal = lambdalplus2mul - 2.0f * mul;

  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // compute the three components of the stress tensor sigma
  if (p_sv){
    // P_SV case
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_zx = sigma_xz;
  }else{
    // SH-case
    sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
    sigma_xz = mul * duxdzl;  // sigma_zy
  }

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // adds gravity contribution to stress (non-symmetric G-tensor term)
  // assumes a constant gravity acting only in (negative) z-direction: g = (0,-g)
  if (gravity && p_sv) {
    realw rhol_g = d_rhostore[offset] * gravity_const;
    // G-tensor contribution added: sigma_ij += rho * [ (s.g) delta_ij - s_i g_j ]
    sigma_xx -= rhol_g * sh_tempz[tx];    // sigma_xx = sigma_xx - rho * g * u_z
    sigma_xz += rhol_g * sh_tempx[tx];    // sigma_xz = sigma_xz + rho * g * u_x
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,xixl,xizl,sh_tempx,sh_tempz,sh_wxgll,J,tx,p_sv);
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,gammaxl,gammazl,sh_tempx,sh_tempz,sh_wxgll,I,tx,p_sv);
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }

// counts:
// + 8 FLOP
//
// + 2 float * 25 threads = 50 BYTE


// counts:
// -----------------
// total of: 790 FLOP per thread
//           ~ 128 * 790 = 101120 FLOP per block
//
//           11392 BYTE DRAM accesses per block
//
// arithmetic intensity: 101120 FLOP / 11392 BYTES ~ 8.9 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 883146240 FLOPS (Single) floating-point operations for 20736 elements
//          -> 42590 FLOP per block
// arithmetic intensity: 42590 FLOP / 11392 BYTES ~ 3.74 FLOP/BYTE
//
// roofline model: Kepler K20x
// ---------------------------
//   for a Kepler K20x card, the peak single-precision performance is about 3.95 TFlop/s.
//   global memory access has a bandwidth of ~ 250 GB/s.
//
//   memory bandwidth: 250 GB/s
//   single-precision peak performance: 3.95 TFlop/s -> corner arithmetic intensity = 3950./250. ~ 15.8 flop/byte
//
//   elastic kernel has an arithmetic intensity of: hand-counts   ~ 8.9 flop/byte
//                                                  nvprof-counts ~ 42590./11392. flop/byte = 3.74 flop/byte
//
//   -> we can only achieve about: (hand-counts)   56% of the peak performance
//                                 (nvprof-counts) 24% of the peak performance -> 935.0 GFlop/s
//
// roofline model: Tesla K20c (Kepler architecture: http://www.nvidia.com/content/tesla/pdf/Tesla-KSeries-Overview-LR.pdf)
// ---------------------------
//   memory bandwidth: 208 GB/s
//   single-precision peak performance: 3.52 TFlop/s -> corner arithmetic intensity = 3520 / 208 ~ 16.9 flop/byte
//
//   we can only achieve about: (hand-counts)   52% of the peak performance
//                              (nvprof-counts) 22% of the peak performance -> 779.0 GFlop/s - measured: 647.3 GFlop/s


} // kernel_2_noatt_iso_impl()


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - elastic isotropic compute forces kernel with PML (no attenuation)

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_PML_impl(const int nb_blocks_to_compute,
                            const int* d_ibool,
                            const int* d_phase_ispec_inner_elastic,
                            const int num_phase_ispec_elastic,
                            const int d_iphase,
                            realw_const_p d_displ,
                            realw_p d_accel,
                            realw* d_xix,realw* d_xiz,
                            realw* d_gammax,realw* d_gammaz,
                            realw_const_p d_hprime_xx,
                            realw_const_p d_hprimewgll_xx,
                            realw_const_p wxgll,
                            realw* d_kappav,
                            realw* d_muv,
                            const int simulation_type,
                            const int p_sv,
                            const int* d_spec_to_pml,
                            int NSPEC_PML_X,
                            int NSPEC_PML_Z,
                            realw deltat,
                            realw* PML_dux_dxl_old,
                            realw* PML_dux_dzl_old,
                            realw* PML_duz_dxl_old,
                            realw* PML_duz_dzl_old,
                            realw* d_displ_elastic_old,
                            realw* d_rmemory_dux_dx,
                            realw* d_rmemory_dux_dx2,
                            realw* d_rmemory_duz_dx,
                            realw* d_rmemory_duz_dx2,
                            realw* d_rmemory_dux_dz,
                            realw* d_rmemory_dux_dz2,
                            realw* d_rmemory_duz_dz,
                            realw* d_rmemory_duz_dz2,
                            realw* d_rmemory_displ_elastic,
                            realw* d_rmemory_displ_elastic2,
                            realw_p d_veloc,
                            const realw* d_rhostore,
                            realw* alphax_store,
                            realw* alphaz_store,
                            realw* betax_store,
                            realw* betaz_store,
                            const int gravity,
                            const realw gravity_const){

// elastic compute kernel without attenuation for isotropic elements with PML

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  int tx = threadIdx.x;

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  //realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;
  realw sum_terms1,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // PML
  int ispec_pml;
  int offset_pml,offset_local_pml;
  realw alpha1,beta1,alphax,betax,alphaz,betaz;
  realw c1,c2;
  realw r1,r2,r3,r4,r5,r6,r7,r8;
  realw r9_x,r9_z,r10_x,r10_z;
  realw rhol,rho_times_jacobianl;
  realw A1,A2,A3,A4;
  realw coef0_1,coef1_1,coef2_1;
  realw coef0_2,coef1_2,coef2_2;
  realw coef0_3,coef1_3,coef2_3;
  realw coef0_4,coef1_4,coef2_4;
  realw pml_contrib_x,pml_contrib_z;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  ispec_pml = d_spec_to_pml[working_element] - 1;

  // checks if element is inside the PML
  if (ispec_pml < 0) return;

  // limits thread ids to range [0,25-1]
  if (tx >= NGLL2 ) tx = tx - NGLL2 ;

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  kappal = d_kappav[offset];
  mul = d_muv[offset];
  rhol = d_rhostore[offset];

  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);
  rho_times_jacobianl = rhol * jacobianl;

  // local index
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  // synchronize all the threads
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // PML
  // local PML array index
  offset_pml = ispec_pml*NGLL2 + tx;  // ispec_pml elements in range [0,NSPEC_PML-1]
  offset_local_pml = (ispec_pml-(NSPEC_PML_X + NSPEC_PML_Z))*NGLL2 + tx; // local pml elements in range [0,NSPEC_PML_XZ-1]

  // coefficients
  alphax = alphax_store[offset_pml];
  betax  = betax_store[offset_pml];
  alphaz = alphaz_store[offset_pml];
  betaz  = betaz_store[offset_pml];

  // note: see pml_init.F90, compare to routine define_PML_coefficients()
  //       (starting around line 900);
  //       K_MIN_PML must be == 1.0 and K_MAX_PML == 1.0,
  //       thus, K_x == K_z == K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER == 1
  if (ispec_pml < NSPEC_PML_X){
    // in CPML_X_ONLY
    //
    // for CPML_X_ONLY:
    //   alpha1  == alpha_x
    //   alpha_z == 0
    //
    //   beta1  == beta_x  == alpha_x + d_x / K_x
    //                        with d_x == d0_x / damping_change_factor_acoustic * abscissa_normalized**NPOWER
    //                             K_x == K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
    //                     == alpha_x + (d0_x / damping_change_factor_acoustic * abscissa_normalized**NPOWER) /
    //                                  (K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER)
    //   beta_z == 0
    //
    alpha1 = alphax;
    beta1 = betax;
    alphaz = 0.f;
    betaz = 0.f;
  } else if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_Z_ONLY region
    //
    // for CPML_Z_ONLY:
    //   alpha1  == alpha_z == ALPHA_MAX * (1 - abscissa)
    //   alpha_x == 0
    //
    //   beta1  == beta_z  == alpha_z + d_z / K_z
    //                     == alpha_z + (2 * d0_z * abscissa**2)
    //   beta_x == 0
    alpha1 = alphaz;
    beta1  = betaz;
    alphax = 0.f;
    betax = 0.f;
  } else {
    // in CPML_XZ region
    alpha1 = alphaz;
    beta1  = betaz;
    //alphax = alphax;
    //betax  = betax;
  }

  // Update memory variables of derivatives
  //
  // coefficients
  // note: see file pml_compute_memory_variables.f90, routine pml_compute_memory_variables_elastic()
  //       assumes Newmark time scheme updates (line ~433)
  //         ! non-rotated, element aligns with x/y/z-coordinates
  //         rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
  //                                           coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
  //         ..
  //
  //       with coefficients coef0_zx_1,.. computed by routine lik_parameter_computation()
  //       and routine compute_coef_convolution() in file pml_compute.f90:
  //         with  c1 == exp( - 1/2 * alpha1 * dt )     alpha1 being either alpha_x or alpha_z
  //               c2 == exp( - 1/2 * beta1  * dt )     beta1  being either beta_x or beta_z
  //
  //         for CPML_X_ONLY region: alpha_x /= 0, alpha_z == 0, beta_x /= 0, beta_z == 0
  //           coef0_zx_1 == exp( - alpha_x * dt ) == c1**2
  //           coef0_zx_2 == exp( - beta_x  * dt ) == c2**2
  //
  //         for CPML_Z_ONLY region: alpha_x == 0, alpha_z /= 0, beta_x == 0, beta_z /= 0)
  //           coef0_zx_1 == exp( - alpha_z * dt ) == c1**2
  //           coef0_zx_2 == exp( - beta_z  * dt ) == c2**2
  //
  //         and
  //           coef1_zx_1 == (1 - c1)/alpha1    (or 1/2 * dt, if alpha1 <= 10^-5)
  //           coef1_zx_2 == (1 - c2)/beta1     (or 1/2 * dt, if beta1  <= 10^-5)
  //
  //           coef2_zx_1 == coef1_zx_1 * c1    (or coef1_zx_1, if alpha1 <= 10^-5)
  //           coef2_zx_2 == coef1_zx_2 * c2    (or coef1_zx_2, if beta1  <= 10^-5)
  //

  // for all PML regions
  c1 = __expf(-0.5f * deltat * alpha1);
  c2 = __expf(-0.5f * deltat * beta1);

  coef0_1 = c1 * c1;
  if (abs(alpha1) > 0.00001f){
    // coef1_zx_1 == (1 - c1)/alpha1
    // coef2_zx_1 == coef1 * c1
    coef1_1 = (1.f - c1) / alpha1;
    coef2_1 = coef1_1 * c1;
  } else {
    // coef1_zx_1 == 1/2 dt
    // coef2_zx_1 == coef1_zx_1
    coef1_1 = 0.5f * deltat;
    coef2_1 = coef1_1;
  }

  coef0_2 = c2 * c2;
  if (abs(beta1) > 0.00001f){
    // coef1_zx_2 == (1 - c2)/beta1
    // coef2_zx_2 == coef1 * c2
    coef1_2 = (1.f - c2) / beta1;
    coef2_2 = coef1_2 * c2;
  } else {
    // coef1_zx_2 == 1/2 dt
    // coef2_zx_2 == coef1_zx_2
    coef1_2 = 0.5f * deltat;
    coef2_2 = coef1_2;
  }

  if (ispec_pml >= (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_XZ region
    realw c3 = __expf(-0.5f * deltat * betax);
    realw c4 = __expf(-0.5f * deltat * alphax);

    coef0_3 = c3 * c3;
    if (abs(betax) > 0.00001f){
      // coef1_zx_2 == (1 - c3)/betax
      // coef2_zx_2 == coef1 * c3
      coef1_3 = (1.f - c3) / betax;
      coef2_3 = coef1_3 * c3;
    } else {
      // coef1_zx_2 == 1/2 dt
      // coef2_zx_2 == coef1_zx_2
      coef1_3 = 0.5f * deltat;
      coef2_3 = coef1_3;
    }

    coef0_4 = c4 * c4;
    if (abs(alphax) > 0.00001f){
      // coef1_zx_1 == (1 - c4)/alphax
      // coef2_zx_1 == coef1 * c4
      coef1_4 = (1.f - c4) / alphax;
      coef2_4 = coef1_4 * c4;
    } else {
      // coef1_zx_1 == 1/2 dt
      // coef2_zx_1 == coef1_zx_1
      coef1_4 = 0.5f * deltat;
      coef2_4 = coef1_4;
    }
  } else {
    coef0_3 = 0.f;
    coef1_3 = 0.f;
    coef2_3 = 0.f;

    coef0_4 = 0.f;
    coef1_4 = 0.f;
    coef2_4 = 0.f;
  }

  // memory variables update
  // see routine pml_compute_memory_variables_elastic() in file pml_compute_memory_variables.f90 (line ~433)
  //   ! non-rotated, element aligns with x/y/z-coordinates
  //   rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
  //                                          coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
  //   ..
  //
  if (ispec_pml < NSPEC_PML_X){
    // in CPML_X_ONLY region
    // rmemory dux_dx
    r1 = coef0_2 * d_rmemory_dux_dx[offset_pml] + coef1_2 * duxdxl + coef2_2 * PML_dux_dxl_old[offset_pml];
    // rmemory duz_dx
    r2 = coef0_2 * d_rmemory_duz_dx[offset_pml] + coef1_2 * duzdxl + coef2_2 * PML_duz_dxl_old[offset_pml];

    // rmemory dux_dz
    r3 = coef0_1 * d_rmemory_dux_dz[offset_pml] + coef1_1 * duxdzl + coef2_1 * PML_dux_dzl_old[offset_pml];
    // rmemory duz_dz
    r4 = coef0_1 * d_rmemory_duz_dz[offset_pml] + coef1_1 * duzdzl + coef2_1 * PML_duz_dzl_old[offset_pml];
  } else {
    // in CPML_Z_ONLY or in CPML_XZ region
    // rmemory dux_dx
    r1 = coef0_1 * d_rmemory_dux_dx[offset_pml] + coef1_1 * duxdxl + coef2_1 * PML_dux_dxl_old[offset_pml];
    // rmemory duz_dx
    r2 = coef0_1 * d_rmemory_duz_dx[offset_pml] + coef1_1 * duzdxl + coef2_1 * PML_duz_dxl_old[offset_pml];

    // rmemory dux_dz
    r3 = coef0_2 * d_rmemory_dux_dz[offset_pml] + coef1_2 * duxdzl + coef2_2 * PML_dux_dzl_old[offset_pml];
    // rmemory duz_dz
    r4 = coef0_2 * d_rmemory_duz_dz[offset_pml] + coef1_2 * duzdzl + coef2_2 * PML_duz_dzl_old[offset_pml];
  }
  d_rmemory_dux_dx[offset_pml] = r1;
  d_rmemory_duz_dx[offset_pml] = r2;
  d_rmemory_dux_dz[offset_pml] = r3;
  d_rmemory_duz_dz[offset_pml] = r4;

  if (ispec_pml >= (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_XZ region
    // rmemory dux_dx2
    r5 = coef0_3 * d_rmemory_dux_dx2[offset_local_pml] + coef1_3 * duxdxl + coef2_3 * PML_dux_dxl_old[offset_pml];
    // rmemory duz_dx2
    r6 = coef0_3 * d_rmemory_duz_dx2[offset_local_pml] + coef1_3 * duzdxl + coef2_3 * PML_duz_dxl_old[offset_pml];
    // rmemory dux_dz2
    r7 = coef0_4 * d_rmemory_dux_dz2[offset_local_pml] + coef1_4 * duxdzl + coef2_4 * PML_dux_dzl_old[offset_pml];
    // rmemory duz_dz2
    r8 = coef0_4 * d_rmemory_duz_dz2[offset_local_pml] + coef1_4 * duzdzl + coef2_4 * PML_duz_dzl_old[offset_pml];
    d_rmemory_dux_dx2[offset_local_pml] = r5;
    d_rmemory_duz_dx2[offset_local_pml] = r6;
    d_rmemory_dux_dz2[offset_local_pml] = r7;
    d_rmemory_duz_dz2[offset_local_pml] = r8;
  } else {
    r5 = 0.f;
    r6 = 0.f;
    r7 = 0.f;
    r8 = 0.f;
  } // ispec \in REGION_XZ

  // Update memory variables of displ
  r9_x = coef0_1 * d_rmemory_displ_elastic[2*offset_pml] + coef1_1 * sh_tempx[tx] + coef2_1 * d_displ_elastic_old[2*offset_pml];
  d_rmemory_displ_elastic[2*offset_pml] = r9_x;

  r9_z = coef0_1 * d_rmemory_displ_elastic[2*offset_pml+1] + coef1_1 * sh_tempz[tx] + coef2_1 * d_displ_elastic_old[2*offset_pml+1];
  d_rmemory_displ_elastic[2*offset_pml+1] = r9_z;

  if (ispec_pml >= (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_XZ region
    r10_x = coef0_4 * d_rmemory_displ_elastic2[2*offset_local_pml] + coef1_4 * sh_tempx[tx] + coef2_4 * d_displ_elastic_old[2*offset_pml];
    d_rmemory_displ_elastic2[2*offset_local_pml] = r10_x;
    r10_z = coef0_4 * d_rmemory_displ_elastic2[2*offset_local_pml+1] + coef1_4 * sh_tempz[tx] + coef2_4 * d_displ_elastic_old[2*offset_pml+1];
    d_rmemory_displ_elastic2[2*offset_local_pml+1] = r10_z;
  } else {
    r10_x = 0.f;
    r10_z = 0.f;
  } // ispec \in REGION_XZ

  // Update old derivatives
  PML_dux_dxl_old[offset_pml] = duxdxl;
  PML_dux_dzl_old[offset_pml] = duxdzl;
  PML_duz_dxl_old[offset_pml] = duzdxl;
  PML_duz_dzl_old[offset_pml] = duzdzl;

  d_displ_elastic_old[2*offset_pml] = sh_tempx[tx];
  d_displ_elastic_old[2*offset_pml+1] = sh_tempz[tx];

  // Compute contribution of the PML
  // note: compare to routine pml_compute_accel_contribution_elastic() (line ~323):
  //         accel_elastic_PML(1,i,j) = wxgll(i) * wzgll(j) * fac * &
  //            ( A1 * veloc_elastic(1,iglob) + A2 * dummy_loc(1,i,j) + &
  //              A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML))
  //         accel_elastic_PML(2,i,j) = wxgll(i) * wzgll(j) * fac * &
  //            ( A1 * veloc_elastic(2,iglob) + A2 * dummy_loc(2,i,j) + &
  //              A3 * rmemory_displ_elastic(1,2,i,j,ispec_PML) + A4 * rmemory_displ_elastic(2,2,i,j,ispec_PML))
  //
  //       with factor fac == rhol * jacobianl
  //
  if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_X_ONLY or in CPML_Z_ONLY region
    A1 = beta1 - alpha1;
    A2 = - alpha1 * A1;     // - alpha1 * (beta1 - alpha1)
    A3 = - alpha1 * A2;     // alpha1 * alpha1 * (beta1 - alpha1)
    A4 = 0.f;
  } else {
    // in CPML_XZ region
    realw fac1 = (alphax * alpha1 + alphax*alphax + 2.f * betax * beta1 - 2.f * alphax * (betax + beta1)) / (alpha1 - alphax);
    realw fac2 = (alphax * alpha1 + alpha1*alpha1 + 2.f * betax * beta1 - 2.f * alpha1 * (betax + beta1)) / (alphax - alpha1);

    A1 = 0.5f * (fac1 - alphax + fac2 - alpha1);
    A2 = 0.5f * (alphax*alphax - fac1 * alphax + alpha1*alpha1 - fac2 * alpha1);
    A3 = 0.5f * alpha1 * alpha1 * (fac2 - alpha1);
    A4 = 0.5f * alphax * alphax * (fac1 - alphax);
  }
  pml_contrib_x = sh_wxgll[J] * sh_wxgll[I] * rho_times_jacobianl * (A1 * d_veloc[iglob*2] + A2 * sh_tempx[tx] + A3 * r9_x + A4 * r10_x);
  pml_contrib_z = sh_wxgll[J] * sh_wxgll[I] * rho_times_jacobianl * (A1 * d_veloc[iglob*2+1] + A2 * sh_tempz[tx] + A3 * r9_z + A4 * r10_z);

  // Update derivatives
  // see routine pml_compute_memory_variables_elastic() in file pml_compute_memory_variables.f90 (line ~511):
  //    dux_dxl(i,j) = A5 * PML_dux_dxl(i,j) + A6 * rmemory_dux_dx(i,j,ispec_PML,1) + A7 * rmemory_dux_dx(i,j,ispec_PML,2)
  //    duz_dxl(i,j) = A5 * PML_duz_dxl(i,j) + A6 * rmemory_duz_dx(i,j,ispec_PML,1) + A7 * rmemory_duz_dx(i,j,ispec_PML,2)
  //    dux_dzl(i,j) = A8 * PML_dux_dzl(i,j) + A9 * rmemory_dux_dz(i,j,ispec_PML,1) + A10 * rmemory_dux_dz(i,j,ispec_PML,2)
  //    duz_dzl(i,j) = A8 * PML_duz_dzl(i,j) + A9 * rmemory_duz_dz(i,j,ispec_PML,1) + A10 * rmemory_duz_dz(i,j,ispec_PML,2)
  //
  //  with coefficients from routine lik_parameter_computation() in file pml_compute.f90 (line ~123):
  //  for X_ONLY_TEMP:
  //    A0 == A5  == Kx
  //    A1 == A6  == - Kx * (alpha_x - beta_x)
  //    A2 == A7  == 0
  //  for Z_ONLY_TEMP:
  //    A0 == A8  == 1/Kz
  //    A1 == A9  == 0
  //    A2 == A10 == - 1/Kz * (beta_z - alpha_z)
  //  for XZ_TEMP:
  //    A0 == Kx/Kz
  //    A1 == 1/2 * Kx/Kz * (gamma_x - alpha_x)
  //          where gamma_x == (alpha_x * beta_z + alpha_x**2 + 2 * beta_x * alpha_z - 2 * alpha_x * (beta_x + alpha_z)) / (beta_z - alpha_x)
  //    A2 == 1/2 * Kx/Kz * (gamma_z - beta_z)
  //          where gamma_z == (alpha_x * beta_z + beta_z**2 + 2 * beta_x * alpha_z - 2 * beta_z * (beta_x + alpha_z)) / (alpha_x - beta_z)
  //
  //  note: K_MIN_PML must be == 1.0 and K_MAX_PML == 1.0 for this implementation, thus Kx == Kz == 1
  //
  if (ispec_pml < NSPEC_PML_X){
    // in CPML_X_ONLY region
    // (alpha1 == alpha_x and beta1 == beta_x)
    //
    // dux_dx: uses coefficients for Z_ONLY_TEMP case
    //         with arguments alpha_z -> alpha_x == alpha1, beta_z -> beta_x == beta1 (index_ik 31)
    //         A5 == A0 == 1
    //         A6 == A2 == - (beta_z - alpha_z) == (alpha_z - beta_z)
    //         A7 == A1 == 0
    realw bar_A = (alpha1-beta1);
    duxdxl += bar_A * r1;
    duzdxl += bar_A * r2;
    // dux_dz: uses coefficients for X_ONLY_TEMP case
    //         with arguments alpha_x -> alpha_x == alpha1, beta_x -> beta_x == beta1 (index_ik 13)
    //         A8 == A0 == 1
    //         A9 == A1 == - (alpha_x - beta_x)
    //         A10 == A2 == 0
    duxdzl -= bar_A * r3;
    duzdzl -= bar_A * r4;
  } else if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)) {
    // in CPML_Z_ONLY region
    // (alpha1 == alpha_z and beta1 == beta_z)
    //
    // dux_dx: uses coefficients for X_ONLY_TEMP case
    //         with arguments alpha_x -> alpha_z == alpha1, beta_x -> beta_z == beta1 (index_ik 31)
    //         A5 == A0 == 1
    //         A6 == A1 == - (alpha_x - beta_x)
    //         A7 == A2 == 0
    realw bar_A = (alpha1-beta1);
    duxdxl -= bar_A * r1;
    duzdxl -= bar_A * r2;
    // dux_dz: uses coefficients for Z_ONLY_TEMP case
    //         with arguments alpha_z -> alpha_z == alpha1, beta_z -> beta_z == beta1 (index_ik 13)
    //         A8 == A0 == 1
    //         A9 == A2 == - (beta_z - alpha_z) == (alpha_z - beta_z)
    //         A10 == A1 == 0
    duxdzl += bar_A * r3;
    duzdzl += bar_A * r4;
  } else {
    // in CPML_XZ region
    // (alpha1 == alpha_z, beta1 == beta_z, and alphax == alpha_x, betax == beta_x)
    //
    // dux_dx: uses coefficients for XZ_TEMP case
    //         with arguments alpha_x -> alpha_z == alpha1, beta_x -> beta_z == beta1 (index_ik 31)
    //                        alpha_z -> alpha_x == alphax, beta_z -> beta_x == betax
    //         A5 == A0 == 1
    //         A6 == A1 == 1/2 * (gamma_x - alpha_x)
    //                     gamma_x == (alpha_x * beta_z + alpha_x**2 + 2 * beta_x * alpha_z - 2 * alpha_x * (beta_x + alpha_z)) / (beta_z - alpha_x)
    realw bar_A1 = 0.5f * ((alpha1 * betax + alpha1*alpha1 + 2.f * beta1 * alphax - 2.f * alpha1 * (beta1 + alphax)) / (betax - alpha1) - alpha1);
    //         A7 == A2 == 1/2 * (gamma_z - beta_z)
    //                     gamma_z == (alpha_x * beta_z + beta_z**2 + 2 * beta_x * alpha_z - 2 * beta_z * (beta_x + alpha_z)) / (alpha_x - beta_z)
    realw bar_A2 = 0.5f * ((alpha1 * betax + betax*betax + 2.f * beta1 * alphax - 2.f * betax * ( beta1 + alphax)) / (alpha1 - betax) - betax);
    duxdxl += bar_A1 * r1 + bar_A2 * r5;
    duzdxl += bar_A1 * r2 + bar_A2 * r6;

    // dux_dz: uses coefficients for XZ_TEMP case
    //         with arguments alpha_x -> alpha_x == alphax, beta_x -> beta_x == betax (index_ik 13)
    //                        alpha_z -> alpha_z == alpha1, beta_z -> beta_z == beta1
    //         A8 == A0 == 1
    //         A9 == A2 == 1/2 * (gamma_z - beta_z)
    //                     gamma_z == (alpha_x * beta_z + beta_z**2 + 2 * beta_x * alpha_z - 2 * beta_z * (beta_x + alpha_z)) / (alpha_x - beta_z)
    realw bar_A3 = 0.5f * ((alphax * beta1 + beta1*beta1 + 2.f * betax * alpha1 - 2.f * beta1 * (betax + alpha1)) / (alphax - beta1) - beta1);
    //         A10 == A1 == 1/2 * (gamma_x - alpha_x)
    //                     gamma_x == (alpha_x * beta_z + alpha_x**2 + 2 * beta_x * alpha_z - 2 * alpha_x * (beta_x + alpha_z)) / (beta_z - alpha_x)
    realw bar_A4 = 0.5f * ((alphax * beta1 + alphax*alphax + 2.f * betax * alpha1 - 2.f * alphax * (betax + alpha1)) / (beta1 - alphax) - alphax);
    duxdzl += bar_A3 * r3 + bar_A4 * r7;
    duzdzl += bar_A3 * r4 + bar_A4 * r8;
  }

  // precompute some sums to save CPU time
  //duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // stress calculations
  // isotropic case
  // compute elements with an elastic isotropic rheology

  // conversion from kappa/mu to lambda/mu
  // AXISYM    : kappal = lambdal + TWO_THIRDS * mul
  // non-AXISYM: kappal = lambdal + mul
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // see compute_forces_viscoelastic.F90 (line ~597):
  //    ! stress components:
  //    !   sigma_xx = (lambda+2mu) F^{-1}[s_z/s_x] * dux_dx + lambda duz_dz
  //    !   sigma_zz = (lambda+2mu) F^{-1}[s_x/s_z] * duz_dz + lambda dux_dx
  //    !
  //    !   sigma_zx = mu duz_dx + mu F^{-1}[s_x/s_z] * dux_dz
  //    !            = mu ( duz_dx + F^{-1}[s_x/s_z] * dux_dz)
  //    !   sigma_xz = mu F^{-1}[s_z/s_x] * duz_dx + mu dux_dz
  //    !            = mu ( dux_dz + F^{-1}[s_z/s_x] * duz_dx )
  //    !
  //    ! note that PML_dux_dxl,PML_dux_dzl,.. arrays contain the original, unmodified dux_dx,dux_dz,.. strain values.
  //    !
  //    sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
  //    sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
  //    sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
  //    sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
  //
  // compute the three components of the stress tensor sigma
  if (p_sv){
    // P_SV case
    // note: for elements w/out PML
    //       sigma_xx = lambdalplus2mul * duxdxl + lambdal * duzdzl;
    //       sigma_zz = lambdalplus2mul * duzdzl + lambdal * duxdxl;
    //       sigma_xz = mul * duzdxl_plus_duxdzl;
    //       sigma_zx = sigma_xz;
    sigma_xx = lambdalplus2mul * duxdxl + lambdal * PML_duz_dzl_old[offset_pml];
    sigma_zz = lambdalplus2mul * duzdzl + lambdal * PML_dux_dxl_old[offset_pml];
    sigma_zx = mul * (PML_duz_dxl_old[offset_pml] + duxdzl);
    sigma_xz = mul * (PML_dux_dzl_old[offset_pml] + duzdxl);
  }else{
    // SH-case
    sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
    sigma_xz = mul * duxdzl;  // sigma_zy
  }

  // adds gravity contribution to stress (non-symmetric G-tensor term)
  // assumes a constant gravity acting only in (negative) z-direction: g = (0,-g)
  if (gravity && p_sv) {
    realw rhol_g = d_rhostore[offset] * gravity_const;
    // G-tensor contribution added: sigma_ij += rho * [ (s.g) delta_ij - s_i g_j ]
    sigma_xx -= rhol_g * sh_tempz[tx];    // sigma_xx = sigma_xx - rho * g * u_z
    sigma_xz += rhol_g * sh_tempx[tx];    // sigma_xz = sigma_xz + rho * g * u_x
    // note: sigma_zx and sigma_zz remain unchanged
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,xixl,xizl,sh_tempx,sh_tempz,sh_wxgll,J,tx,p_sv);
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,gammaxl,gammazl,sh_tempx,sh_tempz,sh_wxgll,I,tx,p_sv);
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1 - pml_contrib_x);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3 - pml_contrib_z);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - elastic anisotropic compute forces kernel (no attenuation)

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_ani_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiz,
                        realw* d_gammax,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p wxgll,
                        realw_const_p d_kappav,
                        realw_const_p d_muv,
                        const int simulation_type,
                        const int p_sv,
                        const int* ispec_is_anisotropic,
                        realw* d_c11store,realw* d_c12store,realw* d_c13store,
                        realw* d_c15store,
                        realw* d_c23store,
                        realw* d_c25store,realw* d_c33store,
                        realw* d_c35store,
                        realw* d_c55store,
                        realw* d_rhostore,
                        const int gravity,
                        const realw gravity_const) {

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL2 ) tx = NGLL2-1;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;

  realw c11,c13,c15,c33,c35,c55;
  realw sum_terms1,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
                              // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // stress calculations
  if (ispec_is_anisotropic[working_element]){
    // full anisotropic case
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];

    // compute the three components of the stress tensor sigma (full anisotropy)
    if (p_sv){
      // P_SV case
      sigma_xx = c11*duxdxl + c13*duzdzl + c15*duzdxl_plus_duxdzl;
      sigma_zz = c13*duxdxl + c33*duzdzl + c35*duzdxl_plus_duxdzl;
      sigma_xz = c15*duxdxl + c35*duzdzl + c55*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = c55 * duxdxl;  // assumes c55 == mu, and still isotropic in both directions - no anisotropy implemented yet...
      sigma_xz = c55 * duxdzl;
    }
  }else{
    // isotropic case

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    // original
    //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    //lambdal = lambdalplus2mul - 2.0f * mul;

    // new
    lambdal = kappal - mul;
    lambdalplus2mul = kappal + mul;

    // compute the three components of the stress tensor sigma
    if (p_sv){
      // P_SV case
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
      sigma_xz = mul * duxdzl;  // sigma_zy
    }
  }

  // adds gravity contribution to stress (non-symmetric G-tensor term)
  // assumes a constant gravity acting only in (negative) z-direction: g = (0,-g)
  if (gravity && p_sv) {
    realw rhol_g = d_rhostore[offset] * gravity_const;
    // G-tensor contribution added: sigma_ij += rho * [ (s.g) delta_ij - s_i g_j ]
    sigma_xx -= rhol_g * sh_tempz[tx];    // sigma_xx = sigma_xx - rho * g * u_z
    sigma_xz += rhol_g * sh_tempx[tx];    // sigma_xz = sigma_xz + rho * g * u_x
    // note: sigma_zx and sigma_zz remain unchanged
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,xixl,xizl,sh_tempx,sh_tempz,sh_wxgll,J,tx,p_sv);
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,gammaxl,gammazl,sh_tempx,sh_tempz,sh_wxgll,I,tx,p_sv);
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  } // threadIdx.x

} // kernel_2_noatt_ani_impl()



/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for viscoelastic domains

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_iso_impl(const int nb_blocks_to_compute,
                      const int* d_ibool,
                      const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                      const int d_iphase,
                      realw_const_p d_displ,
                      realw_p d_accel,
                      realw* d_xix,realw* d_xiz,
                      realw* d_gammax,realw* d_gammaz,
                      realw_const_p d_hprime_xx,
                      realw_const_p d_hprimewgll_xx,
                      realw_const_p wxgll,
                      realw* d_kappav,
                      realw* d_muv,
                      const int simulation_type,
                      const int p_sv,
                      realw_const_p A_newmark_mu,realw_const_p B_newmark_mu,
                      realw_const_p A_newmark_kappa,realw_const_p B_newmark_kappa,
                      realw_p e1,realw_p e11,realw_p e13,
                      realw_p dux_dxl_old,realw_p duz_dzl_old,realw_p dux_dzl_plus_duz_dxl_old,
                      realw* d_rhostore,
                      const int gravity,
                      const realw gravity_const){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .true.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .true. or .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int iglob,offset,offset_align,i_sls;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl,duxdxl_plus_duzdzl;
  realw duxdxl_old,duzdzl_old,duxdzl_plus_duzdxl_old,duxdxl_plus_duzdzl_old;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;
  realw sum_terms1,sum_terms3;

  // attenuation
  realw e1_load[N_SLS],e11_load[N_SLS],e13_load[N_SLS];
  realw e1_sum,e11_sum,e13_sum,a_newmark,b_newmark;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // limits thread ids to range [0,25-1]
  if (tx >= NGLL2 ) tx = tx - NGLL2 ;

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;
  offset_align = working_element*NGLL2 + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  kappal = d_kappav[offset];
  mul = d_muv[offset];

  // attenuation
  for (i_sls=0;i_sls<N_SLS;i_sls++){
    e1_load[i_sls] = e1[N_SLS*offset_align+i_sls];
    e11_load[i_sls] = e11[N_SLS*offset_align+i_sls];
    e13_load[i_sls] = e13[N_SLS*offset_align+i_sls];
  }

  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

  // local index
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __syncthreads();

 // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // compute the three components of the stress tensor sigma
  if (p_sv){
    // P_SV case
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
  }else{
    // SH-case
    sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
    sigma_xz = mul * duxdzl;  // sigma_zy
  }

  // attenuation
  // get the contribution of attenuation and update the memory variables
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duxdxl_old = dux_dxl_old[offset_align];
  duzdzl_old = duz_dzl_old[offset_align];
  duxdxl_plus_duzdzl_old = duxdxl_old + duzdzl_old;
  duxdzl_plus_duzdxl_old = dux_dzl_plus_duz_dxl_old[offset_align];

  e1_sum = 0.f;
  e11_sum = 0.f;
  e13_sum = 0.f;
  for (i_sls=0;i_sls<N_SLS;i_sls++){
    // bulk attenuation
    a_newmark = A_newmark_kappa[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_kappa[N_SLS * offset_align + i_sls];

    e1_load[i_sls] = a_newmark * a_newmark * e1_load[i_sls] + b_newmark * (duxdxl_plus_duzdzl + a_newmark * (duxdxl_plus_duzdzl_old));
    e1_sum += e1_load[i_sls];
    e1[N_SLS*offset_align+i_sls] = e1_load[i_sls];

    // shear attenuation
    a_newmark = A_newmark_mu[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_mu[N_SLS * offset_align + i_sls];

    e11_load[i_sls] = a_newmark * a_newmark * e11_load[i_sls] + b_newmark * (duxdxl - 0.5f*duxdxl_plus_duzdzl + a_newmark * (duxdxl_old-0.5f*duxdxl_plus_duzdzl_old));
    e11_sum += e11_load[i_sls];
    e11[N_SLS*offset_align+i_sls] = e11_load[i_sls];

    e13_load[i_sls] = a_newmark * a_newmark * e13_load[i_sls] + b_newmark * (duzdxl_plus_duxdzl + a_newmark * duxdzl_plus_duzdxl_old);
    e13_sum += e13_load[i_sls];
    e13[N_SLS*offset_align+i_sls] = e13_load[i_sls];
  }

  // add the contribution of the attenuation
  if (p_sv){
    // P_SV case
    sigma_xx += (lambdalplus2mul-mul) * e1_sum + 2.0f * mul * e11_sum;
    sigma_xz += mul * e13_sum;
    sigma_zz += (lambdalplus2mul-mul) * e1_sum - 2.0f * mul * e11_sum;
    sigma_zx = sigma_xz;
  }else{
    // SH-case
    sigma_xx += 0.f;  // attenuation not implemented yet for SH
    sigma_xz += 0.f;
  }

  // adds gravity contribution to stress (non-symmetric G-tensor term)
  // assumes a constant gravity acting only in (negative) z-direction: g = (0,-g)
  if (gravity && p_sv) {
    realw rhol_g = d_rhostore[offset] * gravity_const;
    // G-tensor contribution added: sigma_ij += rho * [ (s.g) delta_ij - s_i g_j ]
    sigma_xx -= rhol_g * sh_tempz[tx];    // sigma_xx = sigma_xx - rho * g * u_z
    sigma_xz += rhol_g * sh_tempx[tx];    // sigma_xz = sigma_xz + rho * g * u_x
  }

  // saves the grad(displ) to use at the next iteration
  dux_dxl_old[offset_align] = duxdxl;
  duz_dzl_old[offset_align] = duzdzl;
  dux_dzl_plus_duz_dxl_old[offset_align] = duzdxl_plus_duxdzl;

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,xixl,xizl,sh_tempx,sh_tempz,sh_wxgll,J,tx,p_sv);
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,gammaxl,gammazl,sh_tempx,sh_tempz,sh_wxgll,I,tx,p_sv);
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }
} // kernel_2_att_iso_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_ani_impl(int nb_blocks_to_compute,
                      const int* d_ibool,
                      const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                      const int d_iphase,
                      realw_const_p d_displ,
                      realw_p d_accel,
                      realw* d_xix,realw* d_xiz,
                      realw* d_gammax,realw* d_gammaz,
                      realw_const_p d_hprime_xx,
                      realw_const_p d_hprimewgll_xx,
                      realw_const_p wxgll,
                      realw_const_p d_kappav,
                      realw_const_p d_muv,
                      const int simulation_type,
                      const int p_sv,
                      const int* ispec_is_anisotropic,
                      realw* d_c11store,realw* d_c12store,realw* d_c13store,
                      realw* d_c15store,
                      realw* d_c23store,
                      realw* d_c25store,realw* d_c33store,
                      realw* d_c35store,
                      realw* d_c55store,
                      realw_const_p A_newmark_mu,realw_const_p B_newmark_mu,
                      realw_const_p A_newmark_kappa,realw_const_p B_newmark_kappa,
                      realw_p e1,realw_p e11,realw_p e13,
                      realw_p dux_dxl_old,realw_p duz_dzl_old,realw_p dux_dzl_plus_duz_dxl_old,
                      realw* d_rhostore,
                      const int gravity,
                      const realw gravity_const) {

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .true.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL2 ) tx = NGLL2-1;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;
  realw c11,c13,c15,c33,c35,c55;
  realw sum_terms1,sum_terms3;

  // attenuation
  int offset_align;
  realw duzdxl_plus_duxdzl,duxdxl_plus_duzdzl;
  realw duxdxl_old,duzdzl_old;
  realw duxdzl_plus_duzdxl_old,duxdxl_plus_duzdzl_old;
  realw e1_load[N_SLS],e11_load[N_SLS],e13_load[N_SLS];
  realw e1_sum,e11_sum,e13_sum,a_newmark,b_newmark;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  // attenuation
  offset_align = working_element*NGLL2 + tx;
  for (int i_sls=0;i_sls<N_SLS;i_sls++){
    e1_load[i_sls] = e1[N_SLS*offset_align+i_sls];
    e11_load[i_sls] = e11[N_SLS*offset_align+i_sls];
    e13_load[i_sls] = e13[N_SLS*offset_align+i_sls];
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = d_xiz[offset]; // all subsequent without to avoid over-use of texture for coalescent access

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // compute elements with an elastic isotropic rheology
  // note: also needed for anisotropy with attenuation case
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  // original
  //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  //lambdal = lambdalplus2mul - 2.0f * mul;
  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // stress calculations
  if (ispec_is_anisotropic[working_element]){
    // full anisotropic case
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];

    // compute the three components of the stress tensor sigma (full anisotropy)
    if (p_sv){
      // P_SV case
      sigma_xx = c11*duxdxl + c13*duzdzl + c15*duzdxl_plus_duxdzl;
      sigma_zz = c13*duxdxl + c33*duzdzl + c35*duzdxl_plus_duxdzl;
      sigma_xz = c15*duxdxl + c35*duzdzl + c55*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = c55 * duxdxl;  // assumes c55 == mu, and still isotropic in both directions - no anisotropy implemented yet...
      sigma_xz = c55 * duxdzl;
    }
  }else{
    // isotropic case
    // compute the three components of the stress tensor sigma
    if (p_sv){
      // P_SV case
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
      sigma_xz = mul * duxdzl;  // sigma_zy
    }
  }

  // attenuation
  // get the contribution of attenuation and update the memory variables
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duxdxl_old = dux_dxl_old[offset_align];
  duzdzl_old = duz_dzl_old[offset_align];
  duxdxl_plus_duzdzl_old = duxdxl_old + duzdzl_old;
  duxdzl_plus_duzdxl_old = dux_dzl_plus_duz_dxl_old[offset_align];

  e1_sum = 0.f;
  e11_sum = 0.f;
  e13_sum = 0.f;
  for (int i_sls=0;i_sls<N_SLS;i_sls++){
    // bulk attenuation
    a_newmark = A_newmark_kappa[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_kappa[N_SLS * offset_align + i_sls];

    e1_load[i_sls] = a_newmark * a_newmark * e1_load[i_sls] + b_newmark * (duxdxl_plus_duzdzl + a_newmark * (duxdxl_plus_duzdzl_old));
    e1_sum += e1_load[i_sls];
    e1[N_SLS*offset_align+i_sls] = e1_load[i_sls];

    // shear attenuation
    a_newmark = A_newmark_mu[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_mu[N_SLS * offset_align + i_sls];

    e11_load[i_sls] = a_newmark * a_newmark * e11_load[i_sls] + b_newmark * (duxdxl - 0.5f*duxdxl_plus_duzdzl + a_newmark * (duxdxl_old-0.5f*duxdxl_plus_duzdzl_old));
    e11_sum += e11_load[i_sls];
    e11[N_SLS*offset_align+i_sls] = e11_load[i_sls];

    e13_load[i_sls] = a_newmark * a_newmark * e13_load[i_sls] + b_newmark * (duzdxl_plus_duxdzl + a_newmark * duxdzl_plus_duzdxl_old);
    e13_sum += e13_load[i_sls];
    e13[N_SLS*offset_align+i_sls] = e13_load[i_sls];
  }

  // add the contribution of the attenuation
  if (p_sv){
    // P_SV case
    sigma_xx += (lambdalplus2mul-mul) * e1_sum + 2.0f * mul * e11_sum;
    sigma_zz += (lambdalplus2mul-mul) * e1_sum - 2.0f * mul * e11_sum;
    sigma_xz += mul * e13_sum;
    sigma_zx = sigma_xz;
  }else{
    // SH-case
    sigma_xx += 0.f;  // attenuation not implemented yet for SH
    sigma_xz += 0.f;
  }

  // adds gravity contribution to stress (non-symmetric G-tensor term)
  // assumes a constant gravity acting only in (negative) z-direction: g = (0,-g)
  if (gravity && p_sv) {
    realw rhol_g = d_rhostore[offset] * gravity_const;
    // G-tensor contribution added: sigma_ij += rho * [ (s.g) delta_ij - s_i g_j ]
    sigma_xx -= rhol_g * sh_tempz[tx];    // sigma_xx = sigma_xx - rho * g * u_z
    sigma_xz += rhol_g * sh_tempx[tx];    // sigma_xz = sigma_xz + rho * g * u_x
    // note: sigma_zx and sigma_zz remain unchanged
  }

  // saves the grad(displ) to use at the next iteration
  dux_dxl_old[offset_align] = duxdxl;
  duz_dzl_old[offset_align] = duzdzl;
  dux_dzl_plus_duz_dxl_old[offset_align] = duzdxl_plus_duxdzl;

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,xixl,xizl,sh_tempx,sh_tempz,sh_wxgll,J,tx,p_sv);
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product_wgll(jacobianl,sigma_xx,sigma_xz,sigma_zx,sigma_zz,gammaxl,gammazl,sh_tempx,sh_tempz,sh_wxgll,I,tx,p_sv);
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  } // threadIdx.x
} // kernel_2_att_ani_impl()


/* ----------------------------------------------------------------------------------------------- */

// note: we used templating to be able to call the same kernel_2 twice for both,
//       forward and backward wavefields. that is, calling it by
//          Kernel_2_noatt_iso_impl<1>
//       and
//          Kernel_2_noatt_iso_impl<3>
//       the templating helped to use textures for forward/backward fields.
//
//       most of this has become obsolete, textures are hardly needed for speedup anymore
//       and the Kernel_2 has become more and more specialized for different cases to
//       reduce register pressure and increase occupancy for better performance.
//       thus, in future we might re-evaluate and remove this template-feature.
//
// "forced" template instantiation
// see: https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
//      https://stackoverflow.com/questions/31705764/cuda-c-using-a-template-function-which-calls-a-template-kernel
//
// for compute_forces_viscoelastic_cuda.cu:
// Kernel_2_noatt_iso_impl<1> needs an explicit instantiation here to be able to link against it from a different .cu file, ..

// isotropic, no attenuation
template __global__ void Kernel_2_noatt_iso_impl<1>(const int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,
                                                    realw*,realw*,const int,const int,const int,const int*,
                                                    realw*,const int,const realw);

template __global__ void Kernel_2_noatt_iso_impl<3>(const int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,
                                                    realw*,realw*,const int,const int,const int,const int*,
                                                    realw*,const int,const realw);

// isotropic, no attenuation, w/ PML
template __global__ void Kernel_2_noatt_iso_PML_impl<1>(const int,const int*,const int*,const int,const int,
                                                        realw_const_p,realw_p,
                                                        realw*,realw*,realw*,realw*,
                                                        realw_const_p,realw_const_p,realw_const_p,
                                                        realw*,realw*,const int,const int,const int*,int,int,realw,
                                                        realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                        realw*,realw*,realw*,realw*,realw*,
                                                        realw_p,const realw*,
                                                        realw*,realw*,realw*,realw*,
                                                        const int,const realw);

template __global__ void Kernel_2_noatt_iso_PML_impl<3>(const int,const int*,const int*,const int,const int,
                                                        realw_const_p,realw_p,
                                                        realw*,realw*,realw*,realw*,
                                                        realw_const_p,realw_const_p,realw_const_p,
                                                        realw*,realw*,const int,const int,const int*,int,int,realw,
                                                        realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                        realw*,realw*,realw*,realw*,realw*,
                                                        realw_p,const realw*,
                                                        realw*,realw*,realw*,realw*,
                                                        const int,const realw);

// anisotropic, no attenuation
template __global__ void Kernel_2_noatt_ani_impl<1>(int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                    const int,const int,const int*,
                                                    realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                    realw*,const int,const realw);

template __global__ void Kernel_2_noatt_ani_impl<3>(int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                    const int,const int,const int*,
                                                    realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                    realw*,const int,const realw);

// isotropic, w/ attenuation
template __global__ void Kernel_2_att_iso_impl<1>(const int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,
                                                  realw*,realw*,const int,const int,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p,
                                                  realw*,const int,const realw);

template __global__ void Kernel_2_att_iso_impl<3>(const int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,
                                                  realw*,realw*,const int,const int,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p,
                                                  realw*,const int,const realw);

// anisotropic, w/ attenuation
template __global__ void Kernel_2_att_ani_impl<1>(int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  const int,const int,const int*,
                                                  realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                  realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p,
                                                  realw*,const int,const realw);

template __global__ void Kernel_2_att_ani_impl<3>(int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  const int,const int,const int*,
                                                  realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                  realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p,
                                                  realw*,const int,const realw);
