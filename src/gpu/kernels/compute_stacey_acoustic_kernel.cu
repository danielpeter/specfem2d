/*
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
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


__global__ void compute_stacey_acoustic_kernel(realw* potential_dot_acoustic,
                                               realw* potential_dot_dot_acoustic,
                                               int* abs_boundary_ispec,
                                               int* abs_boundary_ij,
                                               realw* abs_boundary_jacobian1Dw,
                                               int* d_ibool,
                                               realw* rhostore,
                                               realw* kappastore,
                                               int* ispec_is_acoustic,
                                               int read_abs,
                                               int write_abs,
                                               int UNDO_ATTENUATION_AND_OR_PML,
                                               int compute_wavefield1,
                                               int compute_wavefield2,
                                               int num_abs_boundary_faces,
                                               realw* b_potential_dot_acoustic,
                                               realw* b_potential_dot_dot_acoustic,
                                               realw* b_absorb_potential_left,
                                               realw* b_absorb_potential_right,
                                               realw* b_absorb_potential_top,
                                               realw* b_absorb_potential_bottom,
                                               int* ib_left,
                                               int* ib_right,
                                               int* ib_top,
                                               int* ib_bottom,
                                               int* edge_abs) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,iglob,ispec,num_local;
  realw rhol,kappal,cpl;
  realw jacobianw;
  realw vel;

  if (iface >= num_abs_boundary_faces) return;

  // "-1" from index values to convert from Fortran-> C indexing
  ispec = abs_boundary_ispec[iface]-1;
  if ( ! ispec_is_acoustic[ispec]) return;

  i = abs_boundary_ij[INDEX3(NDIM,NGLLX,0,igll,iface)]-1;
  j = abs_boundary_ij[INDEX3(NDIM,NGLLX,1,igll,iface)]-1;

  //daniel todo: check if we can simplify.
  //       in fortran routine, we set i == NGLLX+1 or j == NGLLX+1
  //       to indicate points which duplicate contributions and can be left out
  //
  //check if the point must be computed
  if (i==NGLLX || j==NGLLX) return;

  iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

  // determines bulk sound speed
  rhol = rhostore[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)];
  kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];
  cpl = sqrt( kappal / rhol );

  // gets associated, weighted jacobian
  jacobianw = abs_boundary_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

  // uses a potential definition of: s = 1/rho grad(chi)
  if (compute_wavefield1) {
    vel = potential_dot_acoustic[iglob] / rhol;
    // Sommerfeld condition
    atomicAdd(&potential_dot_dot_acoustic[iglob],-vel*jacobianw/cpl);
  }

  // adjoint simulations
  if (compute_wavefield2) {
    // we distinguish between undo_attenuation or classical, because undo recomputes it meanwhile classical just reads it
    if (UNDO_ATTENUATION_AND_OR_PML){
      vel = b_potential_dot_acoustic[iglob] / rhol;
      atomicAdd(&b_potential_dot_dot_acoustic[iglob],-vel*jacobianw/cpl);
    }else{
      if (edge_abs[iface] == 1)     { num_local = ib_bottom[iface] - 1;
                                      atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                                -b_absorb_potential_bottom[INDEX2(NGLLX,igll,num_local)]);}
      else if (edge_abs[iface] == 2){ num_local = ib_right[iface] - 1;
                                      atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                                -b_absorb_potential_right[INDEX2(NGLLX,igll,num_local)]);}
      else if (edge_abs[iface] == 3){ num_local = ib_top[iface] - 1;
                                      atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                                -b_absorb_potential_top[INDEX2(NGLLX,igll,num_local)]);}
      else if (edge_abs[iface] == 4){ num_local = ib_left[iface] - 1;
                                      atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                                -b_absorb_potential_left[INDEX2(NGLLX,igll,num_local)]);}
    }
    if (write_abs) {
      // saves boundary values
      if (edge_abs[iface] == 1)      { num_local = ib_bottom[iface] - 1;
                                       b_absorb_potential_bottom[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
      else if (edge_abs[iface] == 2) { num_local = ib_right[iface] - 1;
                                       b_absorb_potential_right[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
      else if (edge_abs[iface] == 3) { num_local = ib_top[iface] - 1;
                                       b_absorb_potential_top[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
      else if (edge_abs[iface] == 4) { num_local = ib_left[iface] - 1;
                                       b_absorb_potential_left[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
    }
  } //if compute_wavefield2
}
