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
realw_texture d_potential_tex;
realw_texture d_potential_dot_dot_tex;
//backward/reconstructed
realw_texture d_b_potential_tex;
realw_texture d_b_potential_dot_dot_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential_dot_dot(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_potential<1>(int x) { return tex1Dfetch(d_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<1>(int x) { return tex1Dfetch(d_potential_dot_dot_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_potential<3>(int x) { return tex1Dfetch(d_b_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<3>(int x) { return tex1Dfetch(d_b_potential_dot_dot_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
// already defined in compute_forces_viscoelastic_cuda.cu
extern realw_texture d_hprime_xx_tex;
//extern realw_texture d_hprimewgll_xx_tex;
extern realw_texture d_wxgll_xx_tex;
#endif


// note on performance optimizations:
//
//   performance tests done:
//   - registers: we were trying to reduce the number of registers, as this is the main limiter for the
//                occupancy of the kernel. however, there is only little difference in register pressure for one "general" kernel
//                or multiple "spezialized" kernels. reducing registers is mainly achieved through the launch_bonds() directive.
//   - branching: we were trying to reduce code branches, such as the if-active check in earlier code versions.
//                reducing the branching helps the compiler to better optimize the executable.
//   - memory accesses: the global memory accesses are avoiding texture reads for coalescent arrays, as this is
//                still faster. thus we were using no __ldg() loads or __restricted__ pointer usage,
//                as those implicitly lead the compiler to use texture reads.
//   - arithmetic intensity: ratio of floating-point operations vs. memory accesses is still low for our kernels.
//                tests with using a loop over elements to re-use the constant arrays (like hprime, wgllwgll,..) and thus
//                increasing the arithmetic intensity failed because the number of registers increased as well.
//                this increased register pressure reduced the occupancy and slowed down the kernel performance.
//   - hiding memory latency: to minimize waiting times to retrieve a memory value from global memory, we put
//                some more calculations into the same code block before calling syncthreads(). this should help the
//                compiler to move independent calculations to wherever it can overlap it with memory access operations.
//                note, especially the if (gravity)-block locations are very sensitive
//                for optimal register usage and compiler optimizations
//

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - acoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_impl(const int nb_blocks_to_compute,
                       const int* d_ibool,
                       const int* d_phase_ispec_inner_acoustic,
                       const int num_phase_ispec_acoustic,
                       const int d_iphase,
                       realw_const_p d_potential_acoustic,
                       realw_p d_potential_dot_dot_acoustic,
                       realw_const_p d_b_potential_acoustic,
                       realw_p d_b_potential_dot_dot_acoustic,
                       const int nb_field,
                       const realw* d_xix, const realw* d_xiz,
                       const realw* d_gammax,const realw* d_gammaz,
                       realw_const_p d_hprime_xx,
                       realw_const_p d_hprimewgll_xx,
                       realw_const_p d_wxgll,
                       const realw* d_rhostore,
                       const realw* d_kappastore,
                       int PML,
                       int* d_spec_to_pml,
                       const int gravity,
                       const realw gravity_const){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^2 = 25 active threads, plus 7 inactive/ghost threads,
  //       because we used memory padding from NGLL^2 = 25 to 32 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J;
  int ispec,iglob,offset;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;
  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl_times_jacobianl,jacobianl;
  realw fac_z;
  realw sum_terms;

  __shared__ realw s_dummy_loc[2*NGLL2];

  __shared__ realw s_temp1[NGLL2];
  __shared__ realw s_temp3[NGLL2];

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

//         counts floating-point operations (FLOP) per thread
//         counts global memory accesses in bytes (BYTES) per block
// 2 FLOP
//
// 0 BYTES

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  ispec = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)] - 1; // array indexing starts at 0

  //checks if element is outside the PML
  if (PML){
    if (d_spec_to_pml[ispec] > 0) return;
  }

  // local padded index
  offset = ispec * NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;

// counts:
// + 7 FLOP
//
// + 2 float * 32 threads = 256 BYTE

#ifdef USE_TEXTURES_FIELDS
  s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
  if (nb_field==2) s_dummy_loc[NGLL2+tx]=texfetch_potential<3>(iglob);
#else
  // changing iglob indexing to match fortran row changes fast style
  s_dummy_loc[tx] = d_potential_acoustic[iglob];
  if (nb_field==2) s_dummy_loc[NGLL2+tx]=d_b_potential_acoustic[iglob];
#endif


// counts:
// + 0 FLOP
//
// + 1 float * 25 threads = 100 BYTE

  // local index
  J = (tx/NGLLX);
  I = (tx-J*NGLLX);

// counts:
// + 3 FLOP
//
// + 0 BYTES

  // note: loads mesh values here to give compiler possibility to overlap memory fetches with some computations;
  //       arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads (arrays accesses are coalescent, thus no need for texture reads)
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] );
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);
  rho_invl_times_jacobianl = jacobianl / d_rhostore[offset];

// counts:
// + 5 FLOP
//
// + 5 float * 32 threads = 160 BYTE

  // loads hprime into shared memory

#ifdef USE_TEXTURES_CONSTANTS
  sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
  sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  // loads hprimewgll into shared memory
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];

  if (threadIdx.x < NGLLX){
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    // changing iglob indexing to match fortran row changes fast style
    sh_wxgll[tx] = d_wxgll[tx];
#endif
  }


// counts:
// + 0 FLOP
//
// + 2 * 1 float * 25 threads = 200 BYTE

  for (int k=0 ; k < nb_field ; k++) {

    // synchronize all the threads (one thread for each of the NGLL grid points of the
    // current spectral element) because we need the whole element to be ready in order
    // to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    // computes first matrix product
    temp1l = 0.f;
    temp3l = 0.f;

    for (int l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      // 1. cut-plane along xi-direction
      temp1l += s_dummy_loc[NGLL2*k+J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
      // 3. cut-plane along gamma-direction
      temp3l += s_dummy_loc[NGLL2*k+l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
    }

// counts:
// + NGLLX * 2 * 6 FLOP = 60 FLOP
//
// + 0 BYTE

    // compute derivatives of ux, uy and uz with respect to x, y and z
    // derivatives of potential
    dpotentialdxl = xixl*temp1l + gammaxl*temp3l;
    dpotentialdzl = xizl*temp1l + gammazl*temp3l;

// counts:
// + 2 * 3 FLOP = 6 FLOP
//
// + 0 BYTE

    // form the dot product with the test vector
    s_temp1[tx] = sh_wxgll[J]*rho_invl_times_jacobianl * (dpotentialdxl*xixl  + dpotentialdzl*xizl)  ;
    s_temp3[tx] = sh_wxgll[I]*rho_invl_times_jacobianl * (dpotentialdxl*gammaxl + dpotentialdzl*gammazl)  ;

// counts:
// + 2 * 6 FLOP = 12 FLOP
//
// + 2 BYTE

    // synchronize all the threads (one thread for each of the NGLL grid points of the
    // current spectral element) because we need the whole element to be ready in order
    // to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    sum_terms = 0.f;
    for (int l=0;l<NGLLX;l++) {
      //assumes hprimewgll_xx = hprimewgll_zz
      sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
    }

// counts:
// + NGLLX * 11 FLOP = 55 FLOP
//
// + 0 BYTE

    // adds gravity contribution
    // gravity volume term: 1/kappa grad(chi) . g = - g/kappa \partial_z chi
    // assumes a constant gravity acting only in (negative) z-direction
    if (gravity) {
      fac_z = - gravity_const / d_kappastore[ispec*NGLL2 + tx];
      sum_terms += jacobianl * sh_wxgll[I] * sh_wxgll[J] * fac_z * dpotentialdzl;
    }

    // assembles potential array
    if (k==0) {
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    } else {
      atomicAdd(&d_b_potential_dot_dot_acoustic[iglob],sum_terms);
    }
// counts:
// + 1 FLOP
//
// + 1 float * 25 threads = 100 BYTE

// -----------------
// total of: 149 FLOP per thread
//           ~ 32 * 149 = 4768 FLOP per block
//
//           818 BYTE DRAM accesses per block
//
//           -> arithmetic intensity: 4768 FLOP / 818 BYTES ~ 5.83 FLOP/BYTE (hand-count)
  } // nb_field loop
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - acoustic compute forces kernel with PML

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_PML_impl(const int nb_blocks_to_compute,
                           const int* d_ibool,
                           const int* d_phase_ispec_inner_acoustic,
                           const int num_phase_ispec_acoustic,
                           const int d_iphase,
                           realw_const_p d_potential_acoustic,
                           realw_p d_potential_dot_dot_acoustic,
                           const realw* d_xix, const realw* d_xiz,
                           const realw* d_gammax,const realw* d_gammaz,
                           realw_const_p d_hprime_xx,
                           realw_const_p d_hprimewgll_xx,
                           realw_const_p d_wxgll,
                           const realw* d_rhostore,
                           int* d_spec_to_pml,
                           int NSPEC_PML_X,
                           int NSPEC_PML_Z,
                           realw deltat,
                           realw* PML_dpotentialdxl_old,
                           realw* PML_dpotentialdzl_old,
                           realw* d_potential_old,
                           realw* d_rmemory_acoustic_dux_dx,
                           realw* d_rmemory_acoustic_dux_dz,
                           realw* d_rmemory_acoustic_dux_dx2,
                           realw* d_rmemory_acoustic_dux_dz2,
                           realw* d_rmemory_pot_acoustic,
                           realw* d_rmemory_pot_acoustic2,
                           realw_p potential_dot,
                           realw* d_kappastore,
                           realw* alphax_store,
                           realw* alphaz_store,
                           realw* betax_store,
                           realw* betaz_store,
                           const int gravity,
                           const realw gravity_const){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  int tx = threadIdx.x;

  int I,J;
  int ispec,iglob,offset;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;
  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl_times_jacobianl,jacobianl;
  realw fac_z;
  realw gravity_term;
  realw sum_terms;

  __shared__ realw s_dummy_loc[NGLL2];

  __shared__ realw s_temp1[NGLL2];
  __shared__ realw s_temp3[NGLL2];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // PML
  int ispec_pml;
  int offset_pml,offset_local_pml;
  realw alpha1,beta1,alphax,betax,alphaz,betaz;
  realw c1,c2;
  realw r1,r2,r3,r4,r5,r6;
  realw rhol,kappal;
  realw fac,A1,A2,A3,A4;
  realw coef0_1,coef1_1,coef2_1;
  realw coef0_2,coef1_2,coef2_2;
  realw coef0_3,coef1_3,coef2_3;
  realw coef0_4,coef1_4,coef2_4;
  realw pml_contrib;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  ispec = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)] - 1; // array indexing starts at 0
  ispec_pml = d_spec_to_pml[ispec] - 1;

  // checks if element is inside the PML
  if (ispec_pml < 0) return;

  // local padded index
  offset = ispec * NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;

#ifdef USE_TEXTURES_FIELDS
  s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
  // changing iglob indexing to match fortran row changes fast style
  s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif

  // local index
  J = (tx/NGLLX);
  I = (tx-J*NGLLX);

  // calculates jacobian
  xixl = get_global_cr( &d_xix[offset] );
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  rhol = d_rhostore[offset];
  kappal = d_kappastore[ispec*NGLL2 + tx]; // non-padded

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);
  rho_invl_times_jacobianl = jacobianl / rhol;

  // loads hprime into shared memory
#ifdef USE_TEXTURES_CONSTANTS
  sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
  sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  // loads hprimewgll into shared memory
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];

  if (threadIdx.x < NGLLX){
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    // changing iglob indexing to match fortran row changes fast style
    sh_wxgll[tx] = d_wxgll[tx];
#endif
  }

  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp3l = 0.f;

  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    // 1. cut-plane along xi-direction
    temp1l += s_dummy_loc[J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
    // 3. cut-plane along gamma-direction
    temp3l += s_dummy_loc[l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
  }

  // compute derivatives of ux, uy and uz with respect to x and z
  // derivatives of potential
  dpotentialdxl = xixl*temp1l + gammaxl*temp3l;
  dpotentialdzl = xizl*temp1l + gammazl*temp3l;

  // adds gravity contribution
  // gravity volume term: 1/kappa grad(chi) . g = - g/kappa \partial_z chi
  // assumes a constant gravity acting only in (negative) z-direction;
  // note: uses the potential gradient before PML convolution (dpotentialdzl is modified below)
  gravity_term = 0.f;
  if (gravity) {
    fac_z = - gravity_const / kappal;
    gravity_term = jacobianl * sh_wxgll[I] * sh_wxgll[J] * fac_z * dpotentialdzl;
  }

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
  // note: see file pml_compute_memory_variables.f90, routine pml_compute_memory_variables_acoustic()
  //       assumes Newmark time scheme updates (line ~180)
  //         rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
  //                                                    coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
  //         rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
  //                                                    coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
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
  // see routine pml_compute_memory_variables_acoustic() in file pml_compute_memory_variables.f90 (line ~180)
  //   rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
  //                                              coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
  //   ..
  //
  //   rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
  //                                              coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
  //   ..
  //
  if (ispec_pml < NSPEC_PML_X){
    // in CPML_X_ONLY region
    // rmemory dux_dx
    r1 = coef0_2 * d_rmemory_acoustic_dux_dx[offset_pml] + coef1_2 * dpotentialdxl + coef2_2 * PML_dpotentialdxl_old[offset_pml];
    // rmemory dux_dz
    r2 = coef0_1 * d_rmemory_acoustic_dux_dz[offset_pml] + coef1_1 * dpotentialdzl + coef2_1 * PML_dpotentialdzl_old[offset_pml];
  } else {
    // in CPML_Z_ONLY or in CPML_XZ region
    // rmemory dux_dx
    r1 = coef0_1 * d_rmemory_acoustic_dux_dx[offset_pml] + coef1_1 * dpotentialdxl + coef2_1 * PML_dpotentialdxl_old[offset_pml];
    // rmemory dux_dz
    r2 = coef0_2 * d_rmemory_acoustic_dux_dz[offset_pml] + coef1_2 * dpotentialdzl + coef2_2 * PML_dpotentialdzl_old[offset_pml];
  }
  d_rmemory_acoustic_dux_dx[offset_pml] = r1;
  d_rmemory_acoustic_dux_dz[offset_pml] = r2;

  if (ispec_pml >= (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_XZ region
    // rmemory dux_dx2
    r3 = coef0_3 * d_rmemory_acoustic_dux_dx2[offset_local_pml] + coef1_3 * dpotentialdxl + coef2_3 * PML_dpotentialdxl_old[offset_pml];
    // rmemory dux_dz2
    r4 = coef0_4 * d_rmemory_acoustic_dux_dz2[offset_local_pml] + coef1_4 * dpotentialdzl + coef2_4 * PML_dpotentialdzl_old[offset_pml];
    d_rmemory_acoustic_dux_dx2[offset_local_pml] = r3;
    d_rmemory_acoustic_dux_dz2[offset_local_pml] = r4;
  } else {
    r3 = 0.f;
    r4 = 0.f;
  } // ispec \in REGION_XZ

  // Update memory variables of potential
  r5 = coef0_1 * d_rmemory_pot_acoustic[offset_pml] + coef1_1 * s_dummy_loc[tx] + coef2_1 * d_potential_old[offset_pml];
  d_rmemory_pot_acoustic[offset_pml] = r5 ;

  if (ispec_pml >= (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_XZ region
    r6 = coef0_4 * d_rmemory_pot_acoustic2[offset_local_pml] + coef1_4 * s_dummy_loc[tx] + coef2_4 * d_potential_old[offset_pml];
    d_rmemory_pot_acoustic2[offset_local_pml] = r6;
  } else {
    r6 = 0.f;
  } // ispec \in REGION_XZ

  // Update old potential
  PML_dpotentialdxl_old[offset_pml] = dpotentialdxl;
  PML_dpotentialdzl_old[offset_pml] = dpotentialdzl;
  d_potential_old[offset_pml] = s_dummy_loc[tx];

  // Compute contribution of the PML to second derivative of potential
  // note: compare to routine pml_compute_accel_contribution_acoustic():
  //          potential_dot_dot_acoustic_PML(i,j)= wxgll(i) * wzgll(j) * fac * &
  //                                               (A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
  //                                                A3 * rmemory_potential_acoustic(1,i,j,ispec_PML) + &
  //                                                A4 * rmemory_potential_acoustic(2,i,j,ispec_PML))
  if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)){
    // in CPML_X_ONLY or in CPML_Z_ONLY region
    A1 = beta1 - alpha1;
    A2 = - alpha1 * A1;     // - alpha1 * (beta1 - alpha1)
    A3 = - alpha1 * A2;     // alpha1 * alpha1 * (beta1 - alpha1)
    A4 = 0.f;
    //pml_contrib = sh_wxgll[J] * sh_wxgll[I] * fac * (A1 * potential_dot[iglob] + A2 * s_dummy_loc[tx] + A3 * r5);
  } else {
    // in CPML_XZ region
    realw fac1 = (alphax * alpha1 + alphax*alphax + 2.f * betax * beta1 - 2.f * alphax * (betax + beta1)) / (alpha1 - alphax);
    realw fac2 = (alphax * alpha1 + alpha1*alpha1 + 2.f * betax * beta1 - 2.f * alpha1 * (betax + beta1)) / (alphax - alpha1);

    A1 = 0.5f * (fac1 - alphax + fac2 - alpha1);
    A2 = 0.5f * (alphax*alphax - fac1 * alphax + alpha1*alpha1 - fac2 * alpha1);
    A3 = 0.5f * alpha1 * alpha1 * (fac2 - alpha1);
    A4 = 0.5f * alphax * alphax * (fac1 - alphax);
    //pml_contrib = sh_wxgll[J] * sh_wxgll[I] * fac * (A1 * potential_dot[iglob] + A2 * s_dummy_loc[tx] + A3 * r5 + A4 * r6);
  }
  fac = jacobianl / kappal;    // jacobianl / kappal
  pml_contrib = sh_wxgll[J] * sh_wxgll[I] * fac * (A1 * potential_dot[iglob] + A2 * s_dummy_loc[tx] + A3 * r5 + A4 * r6);

  // Update derivatives
  // see routine pml_compute_memory_variables_acoustic() in file pml_compute_memory_variables.f90 (line ~218):
  //    dux_dxl(i,j) = A5 * PML_dux_dxl(i,j) + A6 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
  //                                           A7 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2)
  //    dux_dzl(i,j) = A8 * PML_dux_dzl(i,j) + A9 *  rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
  //                                           A10 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2)
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
    dpotentialdxl += bar_A * r1;
    // dux_dz: uses coefficients for X_ONLY_TEMP case
    //         with arguments alpha_x -> alpha_x == alpha1, beta_x -> beta_x == beta1 (index_ik 13)
    //         A8 == A0 == 1
    //         A9 == A1 == - (alpha_x - beta_x)
    //         A10 == A2 == 0
    dpotentialdzl -= bar_A * r2;
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
    dpotentialdxl -= bar_A * r1;
    // dux_dz: uses coefficients for Z_ONLY_TEMP case
    //         with arguments alpha_z -> alpha_z == alpha1, beta_z -> beta_z == beta1 (index_ik 13)
    //         A8 == A0 == 1
    //         A9 == A2 == - (beta_z - alpha_z) == (alpha_z - beta_z)
    //         A10 == A1 == 0
    dpotentialdzl += bar_A * r2;
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
    dpotentialdxl += bar_A1 * r1 + bar_A2 * r3;
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
    dpotentialdzl += bar_A3 * r2 + bar_A4 * r4;
  }
  //__syncthreads();  // not needed... still everything thread local

  // form the dot product with the test vector
  s_temp1[tx] = sh_wxgll[J] * rho_invl_times_jacobianl * (dpotentialdxl * xixl    + dpotentialdzl * xizl);
  s_temp3[tx] = sh_wxgll[I] * rho_invl_times_jacobianl * (dpotentialdxl * gammaxl + dpotentialdzl * gammazl);

  __syncthreads();

  sum_terms = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_zz
    sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
  }

  // adds gravity contribution
  if (gravity) {
    sum_terms += gravity_term;
  }

  // assembles potential array (including gravity contribution)
  atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms - pml_contrib);
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - viscoacoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_viscoacoustic_impl(const int nb_blocks_to_compute,
                            const int* d_ibool,
                            const int* d_phase_ispec_inner_acoustic,
                            const int num_phase_ispec_acoustic,
                            const int d_iphase,
                            realw_const_p d_potential_acoustic,
                            realw_p d_potential_dot_dot_acoustic,
                            const realw* d_xix, const realw* d_xiz,
                            const realw* d_gammax,const realw* d_gammaz,
                            realw_const_p d_hprime_xx,
                            realw_const_p d_hprimewgll_xx,
                            realw_const_p d_wxgll,
                            const realw* d_rhostore,
                            realw_p d_e1_acous,
                            const realw* d_A_newmark,
                            const realw* d_B_newmark,
                            realw_p d_sum_forces_old,
                            const realw* d_kappastore,
                            const int gravity,
                            const realw gravity_const){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int I,J;
  int ispec,iglob,offset,offset_align,i_sls;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;
  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl_times_jacobianl,jacobianl;
  realw fac_z;
  realw sum_terms;
  realw sum_forces_old,forces_attenuation,a_newmark;
  realw e1_acous_load[N_SLS];

  __shared__ realw s_dummy_loc[NGLL2];
  __shared__ realw s_temp1[NGLL2];
  __shared__ realw s_temp3[NGLL2];
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  if (bx >= nb_blocks_to_compute ) return;

  ispec = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)] - 1;

  offset = ispec*NGLL2_PADDED + tx;
  offset_align = ispec*NGLL2 + tx;

  iglob = d_ibool[offset] - 1;

#ifdef USE_TEXTURES_FIELDS
  s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
  s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif

  // local index
  J = (tx/NGLLX);
  I = (tx-J*NGLLX);

  xixl = get_global_cr( &d_xix[offset] );
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);
  rho_invl_times_jacobianl = jacobianl / d_rhostore[offset];

  for (i_sls=0;i_sls<N_SLS;i_sls++)  e1_acous_load[i_sls] = d_e1_acous[N_SLS*offset_align+i_sls];

#ifdef USE_TEXTURES_CONSTANTS
  sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
  sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  // loads hprimewgll into shared memory
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];

  if (threadIdx.x < NGLLX){
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    sh_wxgll[tx] = d_wxgll[tx];
#endif
  }

  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp3l = 0.f;

  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    // 1. cut-plane along xi-direction
    temp1l += s_dummy_loc[J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
    // 3. cut-plane along gamma-direction
    temp3l += s_dummy_loc[l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
  }

  dpotentialdxl = xixl*temp1l +  gammaxl*temp3l;
  dpotentialdzl = xizl*temp1l +  gammazl*temp3l;
  s_temp1[tx] = sh_wxgll[J] * rho_invl_times_jacobianl * (dpotentialdxl*xixl  + dpotentialdzl*xizl)  ;
  s_temp3[tx] = sh_wxgll[I] * rho_invl_times_jacobianl * (dpotentialdxl*gammaxl + dpotentialdzl*gammazl)  ;

  __syncthreads();

  sum_terms = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_zz
    sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
  }

  sum_forces_old = d_sum_forces_old[offset_align];
  forces_attenuation = 0.f;

  for (i_sls=0;i_sls<N_SLS;i_sls++){
    a_newmark = d_A_newmark[N_SLS * offset_align + i_sls];
    e1_acous_load[i_sls] = a_newmark * a_newmark * e1_acous_load[i_sls] + d_B_newmark[N_SLS * offset_align + i_sls] * (sum_terms + a_newmark * sum_forces_old);
    forces_attenuation += e1_acous_load[i_sls];
    d_e1_acous[N_SLS*offset_align+i_sls] = e1_acous_load[i_sls];
  }

  d_sum_forces_old[offset_align] = sum_terms;
  sum_terms += forces_attenuation;

  // adds gravity contribution
  // gravity volume term: 1/kappa grad(chi) . g = - g/kappa partial_z chi
  // assumes a constant gravity acting only in (negative) z-direction
  if (gravity) {
    fac_z = - gravity_const / d_kappastore[offset_align];
    sum_terms += jacobianl * sh_wxgll[I] * sh_wxgll[J] * fac_z * dpotentialdzl;
  }

  atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
}

/* ----------------------------------------------------------------------------------------------- */

// note: we used templating to be able to call the same kernel_2 twice for both,
//       forward and backward wavefields. that is, calling it by
//          Kernel_2_acoustic_impl<1>
//       and
//          Kernel_2_acoustic_impl<3>
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
// for compute_forces_acoustic_cuda.cu:
// Kernel_2_acoustic_impl<1> needs an explicit instantiation here to be able to link against it from a different .cu file, ..

template __global__ void Kernel_2_acoustic_impl<1>(const int,const int*,const int*,const int,const int,
                                                   realw_const_p,realw_p,realw_const_p,realw_p,
                                                   const int,const realw*,const realw*,const realw*,const realw*,
                                                   realw_const_p,realw_const_p,realw_const_p,const realw*,const realw*,int,int*,
                                                   const int,const realw);

template __global__ void Kernel_2_acoustic_impl<3>(const int,const int*,const int*,const int,const int,
                                                   realw_const_p,realw_p,realw_const_p,realw_p,
                                                   const int,const realw*,const realw*,const realw*,const realw*,
                                                   realw_const_p,realw_const_p,realw_const_p,const realw*,const realw*,int,int*,
                                                   const int,const realw);

template __global__ void Kernel_2_acoustic_PML_impl<1>(const int,const int*,const int*,const int,const int,
                                                       realw_const_p,realw_p,
                                                       const realw*, const realw*,const realw*,const realw*,
                                                       realw_const_p,realw_const_p,realw_const_p,
                                                       const realw*,int*,int,int,realw,
                                                       realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                       realw_p,realw*,realw*,realw*,realw*,realw*,
                                                       const int,const realw);

template __global__ void Kernel_2_acoustic_PML_impl<3>(const int,const int*,const int*,const int,const int,
                                                       realw_const_p,realw_p,
                                                       const realw*, const realw*,const realw*,const realw*,
                                                       realw_const_p,realw_const_p,realw_const_p,
                                                       const realw*,int*,int,int,realw,
                                                       realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                       realw_p,realw*,realw*,realw*,realw*,realw*,
                                                       const int,const realw);

template __global__ void Kernel_2_viscoacoustic_impl<1>(const int,const int*,const int*,const int,const int,
                                                        realw_const_p,realw_p,
                                                        const realw*, const realw*,const realw*,const realw*,
                                                        realw_const_p,realw_const_p,realw_const_p,
                                                        const realw*,realw_p,const realw*,const realw*,realw_p,
                                                        const realw*,const int,const realw);

template __global__ void Kernel_2_viscoacoustic_impl<3>(const int,const int*,const int*,const int,const int,
                                                        realw_const_p,realw_p,
                                                        const realw*, const realw*,const realw*,const realw*,
                                                        realw_const_p,realw_const_p,realw_const_p,
                                                        const realw*,realw_p,const realw*,const realw*,realw_p,
                                                        const realw*,const int,const realw);

