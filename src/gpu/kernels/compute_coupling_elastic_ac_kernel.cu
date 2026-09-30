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


__global__ void compute_coupling_elastic_ac_kernel(realw* potential_dot_dot_acoustic,
                                                   realw* accel,
                                                   int num_coupling_ac_el_faces,
                                                   int* coupling_ac_el_ispec,
                                                   int* coupling_ac_el_ijk,
                                                   realw* coupling_ac_el_normal,
                                                   realw* coupling_ac_el_jacobian1Dw,
                                                   int* d_ibool,
                                                   const realw* d_xiz,
                                                   const realw* d_gammaz,
                                                   realw_const_p d_hprime_xx,
                                                   int simulation_type,
                                                   int backward_simulation,
                                                   const int PML,
                                                   const int gravity,
                                                   const realw gravity_const,
                                                   const int* d_spec_to_pml,
                                                   const int NSPEC_PML_X,
                                                   const int NSPEC_PML_Z,
                                                   const realw deltat,
                                                   realw_const_p potential_acoustic,
                                                   realw_const_p potential_acoustic_grav,
                                                   realw_const_p potential_dot_acoustic,
                                                   realw_const_p d_potential_old,
                                                   realw* d_rmemory_sfb_potential_ddot_acoustic,
                                                   realw_const_p alphax_store,
                                                   realw_const_p alphaz_store,
                                                   realw_const_p betax_store,
                                                   realw_const_p betaz_store) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,iglob,ispec;
  realw pressure;
  realw nx,nz;
  realw jacobianw;
  realw dux_dxi,dux_dgamma,dpotentialdzl;

  if (iface < num_coupling_ac_el_faces){
    // gets iglob index
    // note: iglob is the same for accel(..) and potential_dot_dot(..)
    //       we take iglob from the acoustic element side, since for PML the indexing requires the acoustic element.

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = coupling_ac_el_ispec[iface*2] - 1;                        // ispec from acoustic element

    i = coupling_ac_el_ijk[INDEX4(NDIM,2,NGLLX,0,0,igll,iface)] - 1;  // (i,j) from acoustic element
    j = coupling_ac_el_ijk[INDEX4(NDIM,2,NGLLX,1,0,igll,iface)] - 1;

    iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

    // gets associated normal on GLL point
    // note: normal points away from acoustic element
    nx = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,0,igll,iface)]; // (1,igll,iface)
    nz = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,1,igll,iface)]; // (2,igll,iface)

    // gets associated, weighted jacobian
    jacobianw = coupling_ac_el_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

    // uses potential chi such that displacement s = 1/rho grad(chi)
    // pressure p = - kappa ( div( s )) then becomes: p = - dot_dot_chi
    pressure = - potential_dot_dot_acoustic[iglob];

    // gravity: uses potential chi such that displacement s = 1/rho grad(chi),
    //          pressure becomes: p = - kappa ( div( s ) ) = - dot_dot_chi + grad(chi) * g
    //
    // g only acting in negative z-direction, i.e. g = - g z_hat and term: + grad(chi) * g = - g \partial_z \chi
    //
    // see compute_coupling_viscoelastic_ac.f90,
    // for details, see: doc/how_to/README_acoustic_gravity_free_surface.md
    //
    // note: skipped for the adjoint fields of SIMULATION_TYPE == 3, where the adjoint definition
    //       pressure^dagger = potential^dagger applies (as in the CPU implementation),
    //       and overwritten below for PML edges (the CPU PML overwrite drops the gravity term as well)
    if (gravity && !(simulation_type == 3 && backward_simulation == 0)){
      // computes term \partial_z \chi on the acoustic element
      dux_dxi = 0.f;
      dux_dgamma = 0.f;
      for (int l = 0; l < NGLLX; l++){
        // derivative along xi: dchi/dxi(i,j) = sum_k potential(k,j) * hprime_xx(i,k)
        dux_dxi += potential_acoustic_grav[d_ibool[INDEX3_PADDED(NGLLX,NGLLX,l,j,ispec)] - 1]
                   * d_hprime_xx[i + l*NGLLX];
        // derivative along gamma: dchi/dgamma(i,j) = sum_k potential(i,k) * hprime_zz(j,k)
        // note: assumes hprime_zz == hprime_xx
        dux_dgamma += potential_acoustic_grav[d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,l,ispec)] - 1]
                      * d_hprime_xx[j + l*NGLLX];
      }
      // mapping derivatives on acoustic element
      int offset = INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec);
      dpotentialdzl = dux_dxi * d_xiz[offset] + dux_dgamma * d_gammaz[offset];

      // perturbed fluid pressure under gravity: p* = - \kappa \nabla \cdot \mathbf{s} = - \partial_t^2 chi + \nabla chi \cdot g
      // with g only acting in negative z-direction g = (0, -g): \nabla chi \cdot g = - g \partial_z \chi
      pressure -= gravity_const * dpotentialdzl;
    }

    if (simulation_type == 3 && backward_simulation == 0){
      // handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
      // adjoint definition: pressure^\dagger = potential^\dagger
      pressure = - pressure;
    }

    // PML
    if (PML) {
      // PML element index
      // takes acoustic side since we access d_potential_old with offset_pml for local PML elements
      int ispec_pml = d_spec_to_pml[ispec] - 1;
      // checks if element is inside the PML
      if (ispec_pml >= 0) {
        realw alpha1,beta1;
        realw coef0_1,coef1_1,coef2_1;
        realw A1,A2,A3;

        // to match offset for local (i,j) index to thread index (tx == i + j * NGLLX from I = (tx-J*NGLLX))
        int tx = i + j * NGLLX;

        // local PML array index (tx == igll)
        int offset_pml = ispec_pml * NGLL2 + tx;  // ispec_pml elements in range [0,NSPEC_PML-1]

        // coefficients
        if (ispec_pml < NSPEC_PML_X){
          // in CPML_X_ONLY
          //   alpha1  == alpha_x
          //   alpha_z == 0
          //
          //   beta1  == beta_x  == alpha_x + d_x / K_x
          //   beta_z == 0
          //
          alpha1 = alphax_store[offset_pml];
          beta1 = betax_store[offset_pml];
        } else if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)){
          // in CPML_Z_ONLY region
          //   alpha1  == alpha_z
          //   alpha_x == 0
          //
          //   beta1  == beta_z  == alpha_z + d_z / K_z
          //   beta_x == 0
          alpha1 = alphaz_store[offset_pml];
          beta1  = betaz_store[offset_pml];
        } else {
          // should not occur - there should be no coupling edge interfaces in the CPML_XZ regions at the corners
          alpha1 = 0.f;
          beta1  = 0.f;
        }

        // for all PML regions
        realw c1 = __expf(-0.5f * deltat * alpha1);

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

        // memory variables update
        // see compute_coupling_viscoelastic_ac.f90 (line ~127):
        //   ! Newmark
        //   rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) = &
        //                coef0_1 * rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) + &
        //                coef1_1 * potential_acoustic(iglob) + coef2_1 * potential_acoustic_old(iglob)
        //
        realw rm = coef0_1 * d_rmemory_sfb_potential_ddot_acoustic[INDEX2(NGLLX,igll,iface)]
                    + coef1_1 * potential_acoustic[iglob] + coef2_1 * d_potential_old[offset_pml];

        d_rmemory_sfb_potential_ddot_acoustic[INDEX2(NGLLX,igll,iface)] = rm;

        // potential update
        // see compute_coupling_viscoelastic_ac.f90 (line ~144):
        //   pressure = - (A0 * potential_dot_dot_acoustic(iglob) + A1 * potential_dot_acoustic(iglob) + &
        //                A2 * potential_acoustic(iglob) + A3 * rmemory_sfb_potential_ddot_acoustic(1,i,j,inum))
        //
        // with coefficients A0,A1,A2,A3 from routine l_parameter_computation(..),
        // note that we require K_x == K_z == 1:
        // for CPML_X_ONLY
        //    A_0 == kappa_x == 1
        //    A_1 == A_0 * (beta_x - alpha_x)
        //    A_2 == - A_0 * alpha_x * (beta_x - alpha_x)
        //    A_3 == A_0 * alpha_x**2 * (beta_x - alpha_x)
        //    A_4 == 0
        // for CPML_Z_ONLY_TEMP
        //    A_0 == kappa_z == 1
        //    A_1 == A_0 * (beta_z - alpha_z)
        //    A_2 == - A_0 * alpha_z * (beta_z - alpha_z)
        //    A_3 == 0
        //    A_4 == A_0 * alpha_z**2 * (beta_z - alpha_z)
        if (ispec_pml < (NSPEC_PML_X + NSPEC_PML_Z)){
          // in CPML_X_ONLY or in CPML_Z_ONLY region
          //A0 = 1.0f;
          A1 = beta1 - alpha1;
          A2 = - alpha1 * A1;     // - alpha1 * (beta1 - alpha1)
          A3 = - alpha1 * A2;     // alpha1 * alpha1 * (beta1 - alpha1)
          //A4 = 0.f;

          /*
          // note: there might be some difference with the CPU update of routine compute_coupling_viscoelastic_ac(), i.e,
          //       for CPML_X_ONLY the coefficient A4 == 0 and for CPML_Z_ONLY the coefficient A3 == 0
          // thus, to match CPU implementation
          if (ispec_pml < NSPEC_PML_X){
            // in CPML_X_ONLY region
            A3 = - alpha1 * A2;
          } else {
            // in CPML_Z_ONLY region
            A3 = 0.f;
          }
          */
        } else {
          // in CPML_XZ region - should not occur
          //A0 = 1.0f; // keeps pressure as is
          A1 = 0.f;
          A2 = 0.f;
          A3 = 0.f;
        }
        // overwrites pressure value
        pressure = - (potential_dot_dot_acoustic[iglob] + A1 * potential_dot_acoustic[iglob] + A2 * potential_acoustic[iglob] + A3 * rm);
      }
    } // PML

    // continuity of displacement and pressure on global point
    //
    // note: Newmark time scheme together with definition of scalar potential:
    //          pressure = - chi_dot_dot
    //          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
    //          pressure at time step [t + delta_t]
    //          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
    //          it means you have to calculate and update the acoustic pressure first before
    //          calculating this term...
    atomicAdd(&accel[iglob*2],jacobianw * nx * pressure);
    atomicAdd(&accel[iglob*2+1],jacobianw * nz * pressure);
  }
}

