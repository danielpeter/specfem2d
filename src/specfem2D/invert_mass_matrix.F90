!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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

  subroutine invert_mass_matrix_init()

!  builds the global mass matrix

  use constants, only: IMAIN,CUSTOM_REAL,NGLLX,NGLLZ,ONE,TWO,TWO_THIRDS,FOUR_THIRDS, &
                       CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ, &
                       IEDGE1,IEDGE2,IEDGE3,IEDGE4,pi

  use specfem_par, only: myrank,any_elastic,any_acoustic,any_poroelastic,any_electromagnetic, &
                         rmass_inverse_elastic, &
                         rmass_inverse_acoustic,rmass_inverse_e1,ATTENUATION_VISCOACOUSTIC,phi_nu1,N_SLS, &
                         time_stepping_scheme,rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic, &
                         rmass_inverse_electromagnetic, &
                         nspec,ibool,wxgll,wzgll,jacobian, &
                         ispec_is_elastic,ispec_is_acoustic,ispec_is_poroelastic,ispec_is_electromagnetic, &
                         phistore,tortstore,rhoarraystore, &
                         rhostore,kappastore, &
                         spermittivitystore,sconductivitystore, &
                         deltatover2,kmato, &
                         gammaz, &
                         AXISYM,is_on_the_axis,coord,wxglj,xiglj, &
                         time_stepping_scheme,STACEY_ABSORBING_CONDITIONS, &
                         ATTENUATION_PERMITTIVITY,ATTENUATION_CONDUCTIVITY, &
                         Qe11_electromagnetic,Qe33_electromagnetic, &
                         Qs11_electromagnetic,Qs33_electromagnetic,f0_electromagnetic


  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML,region_CPML,spec_to_PML, &
                         K_x_store,K_z_store,d_x_store,d_z_store

  use shared_parameters, only: GRAVITY

  implicit none

  ! local parameter
  integer :: ispec,i,j,iglob

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: xxi
  real(kind=CUSTOM_REAL) :: phinu1

  double precision :: rhol,kappal_relaxed
  double precision :: rho_s,rho_f,rho_bar,phi,tort

  integer :: ispec_PML
  logical :: this_element_has_PML

  integer :: i_sls

 ! electromagnetic parameters
  double precision :: permlxx,permlzz,condlxx,condlzz
  double precision :: taus_x,taue_x,taud_x
  double precision :: taus_z,taue_z,taud_z
  double precision, dimension(2):: eperml,perm,cond

  if (myrank == 0) then
    write(IMAIN,*) "  initializing mass matrices"
    call flush_IMAIN()
  endif

  ! initialize mass matrix
  if (any_elastic) then
    rmass_inverse_elastic(:,:) = 0._CUSTOM_REAL
  endif

  if (any_poroelastic) then
    rmass_s_inverse_poroelastic(:) = 0._CUSTOM_REAL
    rmass_w_inverse_poroelastic(:) = 0._CUSTOM_REAL
  endif

  if (any_acoustic) then
    rmass_inverse_acoustic(:) = 0._CUSTOM_REAL
    if (ATTENUATION_VISCOACOUSTIC) rmass_inverse_e1(:,:) = 0._CUSTOM_REAL
  endif

  if (any_electromagnetic) then
    rmass_inverse_electromagnetic(:,:) = 0._CUSTOM_REAL
  endif

  ! computes mass matrix for each element (poroelastic/elastic/acoustic/electromagnetic)
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)

        ! gets density model (elastic or acoustic)
        rhol = rhostore(i,j,ispec)
        kappal_relaxed = kappastore(i,j,ispec)

        if (ispec_is_poroelastic(ispec)) then
          ! material is poroelastic

          !!! PML NOT WORKING YET !!!
          this_element_has_PML = .false.
          if (PML_BOUNDARY_CONDITIONS) then
            if (ispec_is_PML(ispec)) call stop_the_code('PML not implemented yet for poroelastic case')
          endif

          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          rho_bar = (1.d0-phi)*rho_s + phi*rho_f

          ! for the solid mass matrix
          rmass_s_inverse_poroelastic(iglob) = rmass_s_inverse_poroelastic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rho_bar - phi*rho_f/tort)
          ! for the fluid mass matrix
          rmass_w_inverse_poroelastic(iglob) = rmass_w_inverse_poroelastic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rho_bar*rho_f*tort - phi*rho_f*rho_f)/(rho_bar*phi)

        else if (ispec_is_elastic(ispec)) then
          ! for elastic medium

          this_element_has_PML = .false.
          if (PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) this_element_has_PML = .true.

          if (this_element_has_PML) then
            ! PML
            ispec_PML = spec_to_PML(ispec)

            ! mass matrix contribution for PML
            if (region_CPML(ispec) == CPML_X_ONLY) then
              if (AXISYM) then  ! This PML can't be on the axis
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML))
              else ! not axisym
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML))
              endif
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)

            else if (region_CPML(ispec) == CPML_XZ) then
              if (AXISYM) then  ! This corner can't be on the axis
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                     * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
              else ! not axisym
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
              endif
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)

            else if (region_CPML(ispec) == CPML_Z_ONLY) then
              if (AXISYM) then
                if (is_on_the_axis(ispec)) then
                  if (i == 1) then ! First GLJ point
                    xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                    rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                       + xxi*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
                  else
                    rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                       + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
                  endif
                else ! not on the axis
                  rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
                endif
              else ! not axisym
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
              endif
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
            endif

            ! adds additional Newmark term (dt/2 * C) to mass matrix for stabilization
            if (time_stepping_scheme == 1) then
              ! Newmark
              if (region_CPML(ispec) == CPML_X_ONLY) then
                if (AXISYM) then  ! This PML can't be on the axis
                   rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                        + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (d_x_store(i,j,ispec_PML) * deltatover2)
                else ! not axisym
                   rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                        + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (d_x_store(i,j,ispec_PML) * deltatover2)
                endif
                rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)

              else if (region_CPML(ispec) == CPML_XZ) then
                if (AXISYM) then  ! This corner can't be on the axis
                   rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob)  &
                        + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                        * ((d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) + &
                            d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltatover2)
                else ! not axisym
                   rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                        + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                        * ((d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) + &
                            d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltatover2)
                endif
                rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)

              else if (region_CPML(ispec) == CPML_Z_ONLY) then
                if (AXISYM) then
                  if (is_on_the_axis(ispec)) then
                    if (i == 1) then ! First GLJ point
                      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                      rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                         + xxi*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (d_z_store(i,j,ispec_PML)* deltatover2)
                    else
                      rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                         + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                         * (d_z_store(i,j,ispec_PML)* deltatover2)
                    endif
                  else ! not on the axis
                    rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                         + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (d_z_store(i,j,ispec_PML)* deltatover2)
                  endif
                else ! not axisym
                  rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                       + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (d_z_store(i,j,ispec_PML)* deltatover2)
                endif
                rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
              endif
            endif ! Newmark

          else
            ! elastic no PML
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (i == 1) then
                  ! First GLJ point
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                      + xxi*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
                else
                  rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                      + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
                endif
              else
                ! not on the axis
                rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                    + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
              endif
            else
              ! not axisym
              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) &
                      + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
            endif
            rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
          endif

        else if (ispec_is_acoustic(ispec)) then
          ! for acoustic medium

          this_element_has_PML = .false.
          if (PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) this_element_has_PML = .true.

          if (this_element_has_PML) then
            ! PML
            ispec_PML = spec_to_PML(ispec)

            ! mass matrix contribution for PML
            if (region_CPML(ispec) == CPML_X_ONLY) then
              if (AXISYM) then   !! ABAB: This PML cannot be on the axis: it is a right PML
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                     + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                       * (K_x_store(i,j,ispec_PML))
              else ! not axisym
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                     + wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                       * (K_x_store(i,j,ispec_PML))
              endif

            else if (region_CPML(ispec) == CPML_XZ) then
              if (AXISYM) then   !! ABAB: This corner cannot be on the axis
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                     + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                       * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
              else ! not axisym
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                  + wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                    * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
              endif

            else if (region_CPML(ispec) == CPML_Z_ONLY) then
              if (AXISYM) then
                if (is_on_the_axis(ispec)) then
                  if (i == 1) then ! First GLJ point
                    xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                    rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                          + xxi*wxglj(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                            * (K_z_store(i,j,ispec_PML))
                  else
                    rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                          + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)/kappal_relaxed*jacobian(i,j,ispec) &
                            * (K_z_store(i,j,ispec_PML))
                  endif
                else ! not on the axis
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                         + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                           * (K_z_store(i,j,ispec_PML))
                endif
              else ! not axisym
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                     + wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                       * (K_z_store(i,j,ispec_PML))

              endif
            endif

            ! adds additional Newmark term (dt/2 * C) to mass matrix for stabilization
            if (time_stepping_scheme == 1) then
              if (region_CPML(ispec) == CPML_X_ONLY) then
                if (AXISYM) then   !! ABAB: This PML cannot be on the axis: it is a right PML
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                       + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                         * (d_x_store(i,j,ispec_PML) * deltatover2)
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                       + wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                         * (d_x_store(i,j,ispec_PML) * deltatover2)
                endif

              else if (region_CPML(ispec) == CPML_XZ) then
                if (AXISYM) then   !! ABAB: This corner cannot be on the axis
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                       + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                         * ((d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) &
                           + d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltatover2)
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                       + wxgll(i)*wzgll(j)/ kappal_relaxed*jacobian(i,j,ispec) &
                         * ((d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) &
                           + d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltatover2)
                endif

              else if (region_CPML(ispec) == CPML_Z_ONLY) then
                if (AXISYM) then
                  if (is_on_the_axis(ispec)) then
                    if (i == 1) then ! First GLJ point
                      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                      rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                         + xxi*wxglj(i)*wzgll(j)/kappal_relaxed*jacobian(i,j,ispec) &
                           * (d_z_store(i,j,ispec_PML)* deltatover2)
                    else
                      rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                         + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)/kappal_relaxed*jacobian(i,j,ispec) &
                           * (d_z_store(i,j,ispec_PML)* deltatover2)
                    endif
                  else ! not on the axis
                    rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                         + coord(1,iglob)*wxgll(i)*wzgll(j)/kappal_relaxed*jacobian(i,j,ispec) &
                           * (d_z_store(i,j,ispec_PML)* deltatover2)
                  endif
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                       + wxgll(i)*wzgll(j)/kappal_relaxed*jacobian(i,j,ispec) &
                         * (d_z_store(i,j,ispec_PML)* deltatover2)
                endif
              endif
            endif ! Newmark

          else

            ! acoustic no PML
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (i == 1) then
                  ! First GLJ point
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                      + xxi*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / kappal_relaxed

                  if (ATTENUATION_VISCOACOUSTIC) then
                    ! loop over relaxation mechanisms
                    do i_sls = 1,N_SLS
                      phinu1 = 1.0
                      if (time_stepping_scheme > 1) phinu1 = phi_nu1(i,j,ispec,i_sls)
                      rmass_inverse_e1(iglob,i_sls) = rmass_inverse_e1(iglob,i_sls) &
                         + xxi*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / phinu1
                    enddo
                  endif

                else
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                      + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / kappal_relaxed

                  if (ATTENUATION_VISCOACOUSTIC) then
                    ! loop over relaxation mechanisms
                    do i_sls = 1,N_SLS
                      phinu1 = 1.
                      if (time_stepping_scheme > 1) phinu1 = phi_nu1(i,j,ispec,i_sls)
                      rmass_inverse_e1(iglob,i_sls) = rmass_inverse_e1(iglob,i_sls) &
                         + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / phinu1
                    enddo
                  endif

                endif
              else
                ! not on the axis
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                    + coord(1,iglob)*wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal_relaxed

                if (ATTENUATION_VISCOACOUSTIC) then
                  ! loop over relaxation mechanisms
                  do i_sls = 1,N_SLS
                    phinu1 = 1.0
                    if (time_stepping_scheme > 1) phinu1 = phi_nu1(i,j,ispec,i_sls)
                    rmass_inverse_e1(iglob,i_sls) = rmass_inverse_e1(iglob,i_sls) &
                       + coord(1,iglob)*wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / phinu1
                  enddo
                endif

              endif
            else
              ! not axisym
              rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                   + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal_relaxed

              if (ATTENUATION_VISCOACOUSTIC) then
                ! loop over relaxation mechanisms
                do i_sls = 1,N_SLS
                  phinu1 = 1.0
                  if (time_stepping_scheme > 1) phinu1 = phi_nu1(i,j,ispec,i_sls)
                  rmass_inverse_e1(iglob,i_sls) = rmass_inverse_e1(iglob,i_sls) &
                       + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / phinu1
                enddo
              endif

            endif

          endif

        else if (ispec_is_electromagnetic(ispec)) then
         ! electromagnetic medium

          !!! PML NOT WORKING YET !!!
          this_element_has_PML = .false.
          if (PML_BOUNDARY_CONDITIONS) then
            if (ispec_is_PML(ispec)) stop 'PML not implemented yet for electromagnetic case'
          endif

          permlxx = spermittivitystore(1,i,j,ispec)  !e11
          permlzz = spermittivitystore(2,i,j,ispec)  !e33
          condlxx = sconductivitystore(1,i,j,ispec)  !sig11
          condlzz = sconductivitystore(2,i,j,ispec)  !sig33

          ! permittivity (Zener model)
          taud_x = (sqrt(Qe11_electromagnetic(kmato(ispec))**2+1.d0) +1.d0)/ &
                (2.d0*pi*f0_electromagnetic*Qe11_electromagnetic(kmato(ispec)))
          taud_z = (sqrt(Qe33_electromagnetic(kmato(ispec))**2+1.d0) +1.d0)/ &
                (2.d0*pi*f0_electromagnetic*Qe33_electromagnetic(kmato(ispec)))
          taue_x = (sqrt(Qe11_electromagnetic(kmato(ispec))**2+1.d0) -1.d0)/ &
                (2.d0*pi*f0_electromagnetic*Qe11_electromagnetic(kmato(ispec)))
          taue_z = (sqrt(Qe33_electromagnetic(kmato(ispec))**2+1.d0) -1.d0)/ &
                (2.d0*pi*f0_electromagnetic*Qe33_electromagnetic(kmato(ispec)))

          ! conductivity (Kelvin-Voigt model)
          taus_x = Qs11_electromagnetic(kmato(ispec)) / (2.d0*pi*f0_electromagnetic)
          taus_z = Qs33_electromagnetic(kmato(ispec)) / (2.d0*pi*f0_electromagnetic)

          if (ATTENUATION_PERMITTIVITY .and. Qe11_electromagnetic(kmato(ispec)) < 90.d0 ) then
          ! high frequency, with memory variables
          perm(1) = permlxx * taue_x/taud_x
          perm(2) = permlzz * taue_z/taud_z
          else
          perm(1) = permlxx
          perm(2) = permlzz
          endif

          if (ATTENUATION_CONDUCTIVITY) then
          cond(1) = condlxx * taus_x
          cond(2) = condlzz * taus_z
          else
          cond(1) = 0.d0
          cond(2) = 0.d0
          endif

          ! effective permeabilities
          eperml(:) = perm(:) + cond(:)

          rmass_inverse_electromagnetic(1,iglob) = rmass_inverse_electromagnetic(1,iglob) &
                  + wxgll(i)*wzgll(j)*eperml(1)*jacobian(i,j,ispec)
          rmass_inverse_electromagnetic(2,iglob) = rmass_inverse_electromagnetic(2,iglob) &
                  + wxgll(i)*wzgll(j)*eperml(2)*jacobian(i,j,ispec)

       else
          call stop_the_code('Invalid element type found in routine invert_mass_matrix_init()')
       endif
      enddo
    enddo
  enddo ! of do ispec = 1,nspec

  !
  !--- DK and Zhinan Xie: add C Delta_t / 2 contribution to the mass matrix
  !--- DK and Zhinan Xie: in the case of Clayton-Engquist absorbing boundaries;
  !--- DK and Zhinan Xie: see for instance the book of Hughes (1987) chapter 9.
  !--- DK and Zhinan Xie: IMPORTANT: note that this implies that we must have two different mass matrices,
  !--- DK and Zhinan Xie: one per component of the wave field i.e. one per spatial dimension.
  !--- DK and Zhinan Xie: This was also suggested by Jean-Paul Ampuero in 2003.
  !

  ! Stacey contribution for elastic medium
  if (STACEY_ABSORBING_CONDITIONS .and. time_stepping_scheme == 1) then
    ! Newmark scheme contribution to mass matrix
    call invert_mass_matrix_init_Stacey()
  endif

  ! adds acoustic gravity free surface contribution to mass matrix
  if (any_acoustic .and. GRAVITY) then
    call prepare_mass_matrix_acoustic_gravity()
  endif

  end subroutine invert_mass_matrix_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine invert_mass_matrix_init_Stacey()

  use constants, only: IMAIN,CUSTOM_REAL,NGLLX,NGLLZ,ONE,TWO,TWO_THIRDS,FOUR_THIRDS, &
                       IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: myrank,any_elastic,any_acoustic,any_electromagnetic, &
                         rmass_inverse_elastic, &
                         rmass_inverse_acoustic, &
                         rmass_inverse_electromagnetic, &
                         time_stepping_scheme, &
                         ibool,wxgll,wzgll,jacobian, &
                         ispec_is_elastic,ispec_is_acoustic, &
                         rho_vpstore,rho_vsstore, &
                         spermittivitystore,sconductivitystore,inv_magpermeabilitystore, &
                         num_abs_boundary_faces,abs_boundary_ispec, &
                         deltatover2,kmato, &
                         codeabs,codeabs_corner, &
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                         xix,xiz,gammaz,gammax, &
                         time_stepping_scheme,P_SV,STACEY_ABSORBING_CONDITIONS,ATTENUATION_PERMITTIVITY,ATTENUATION_CONDUCTIVITY, &
                         Qe11_electromagnetic,Qe33_electromagnetic,Qs11_electromagnetic,Qs33_electromagnetic,f0_electromagnetic

  ! axisymmetry
  use specfem_par, only: AXISYM,is_on_the_axis,coord,xiglj,wxglj

  implicit none

  ! local parameter
  integer :: ispecabs,ibegin,iend,jbegin,jend,ispec,i,j,iglob

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: rho_vp,rho_vs
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,tx,ty,tz
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: r_xiplus1

  ! electromagnetic parameters
  double precision :: permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability
  double precision :: cpxsquare,cpzsquare
  double precision, dimension(2):: econdl,eperml
  real(kind=CUSTOM_REAL) :: cpxl,cpzl

  ! check if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (time_stepping_scheme /= 1) return  ! only for Newmark time scheme

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  adding Stacey contributions to mass matrix"
    call flush_IMAIN()
  endif

  ! elastic medium
  if (any_elastic) then

    do ispecabs = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(ispecabs)

      if (ispec_is_elastic(ispec)) then

        !--- left absorbing boundary
        if (codeabs(IEDGE4,ispecabs)) then
          i = 1
          do j = 1,NGLLZ
            iglob = ibool(i,j,ispec)

            rho_vp = rho_vpstore(i,j,ispec)
            rho_vs = rho_vsstore(i,j,ispec)

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D

            weight = jacobian1D * wzgll(j)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) weight = weight * coord(1,iglob)

            ! Clayton-Engquist condition if elastic
            if (P_SV) then
              ! P_SV-case
              vx = deltatover2
              vz = deltatover2

              vn = nx * vx + nz * vz

              tx = rho_vp * vn * nx + rho_vs * (vx - vn * nx)
              tz = rho_vp * vn * nz + rho_vs * (vz - vn * nz)

              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + tx*weight
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(2,iglob) + tz*weight
            else
              ! SH-case
              vy = deltatover2
              ty = rho_vs * vy
              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + ty*weight
              ! ficticous for SH case, but to be save when inverting
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
            endif
          enddo
        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if (codeabs(IEDGE2,ispecabs)) then
          i = NGLLX
          do j = 1,NGLLZ
            iglob = ibool(i,j,ispec)

            ! for analytical initial plane wave for Bielak's conditions
            ! left or right edge, horizontal normal vector

            rho_vp = rho_vpstore(i,j,ispec)
            rho_vs = rho_vsstore(i,j,ispec)

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D

            weight = jacobian1D * wzgll(j)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) weight = weight * coord(1,iglob)

            ! Clayton-Engquist condition if elastic
            if (P_SV) then
              ! P_SV-case
              vx = deltatover2
              vz = deltatover2

              vn = nx * vx + nz * vz

              tx = rho_vp * vn * nx + rho_vs * (vx - vn * nx)
              tz = rho_vp * vn * nz + rho_vs * (vz - vn * nz)

              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + tx*weight
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(2,iglob) + tz*weight
            else
              ! SH-case
              vy = deltatover2
              ty = rho_vs * vy
              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + ty*weight
              ! ficticous for SH case, but to be save when inverting
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
            endif
          enddo
        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if (codeabs(IEDGE1,ispecabs)) then
          j = 1
          do i = 1,NGLLX
            ! exclude corners to make sure there is no contradiction on the normal
            ! for Stacey absorbing conditions but not for incident plane waves;
            ! thus subtract nothing i.e. zero in that case
            if ((codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX)) cycle

            iglob = ibool(i,j,ispec)

            rho_vp = rho_vpstore(i,j,ispec)
            rho_vs = rho_vsstore(i,j,ispec)

            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D

            weight = jacobian1D * wxgll(i)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) then
              ! using GLJ quadrature along the radius for elements on the axis
              if (is_on_the_axis(ispec)) then
                ! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
                if (i == 1) then
                  ! first GLJ point: limit of r/(1+xi) at the axis
                  r_xiplus1 = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                else
                  r_xiplus1 = coord(1,iglob) / (xiglj(i) + ONE)
                endif
                weight = jacobian1D * wxglj(i) * r_xiplus1
              else
                weight = weight * coord(1,iglob)
              endif
            endif

            ! Clayton-Engquist condition if elastic
            if (P_SV) then
              ! P_SV-case
              vx = deltatover2
              vz = deltatover2

              vn = nx * vx + nz * vz

              tx = rho_vp * vn * nx + rho_vs * (vx - vn * nx)
              tz = rho_vp * vn * nz + rho_vs * (vz - vn * nz)

              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + tx*weight
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(2,iglob) + tz*weight
            else
              ! SH-case
              vy = deltatover2
              ty = rho_vs * vy
              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + ty*weight
              ! ficticous for SH case, but to be save when inverting
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
            endif
          enddo
        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if (codeabs(IEDGE3,ispecabs)) then
          j = NGLLZ
          do i = 1,NGLLX
            ! exclude corners to make sure there is no contradiction on the normal
            ! for Stacey absorbing conditions but not for incident plane waves;
            ! thus subtract nothing i.e. zero in that case
            if ((codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX)) cycle

            iglob = ibool(i,j,ispec)

            rho_vp = rho_vpstore(i,j,ispec)
            rho_vs = rho_vsstore(i,j,ispec)

            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D

            weight = jacobian1D * wxgll(i)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) then
              ! using GLJ quadrature along the radius for elements on the axis
              if (is_on_the_axis(ispec)) then
                ! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
                if (i == 1) then
                  ! first GLJ point: limit of r/(1+xi) at the axis
                  r_xiplus1 = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                else
                  r_xiplus1 = coord(1,iglob) / (xiglj(i) + ONE)
                endif
                weight = jacobian1D * wxglj(i) * r_xiplus1
              else
                weight = weight * coord(1,iglob)
              endif
            endif

            ! Clayton-Engquist condition if elastic
            if (P_SV) then
              ! P_SV-case
              vx = deltatover2
              vz = deltatover2

              vn = nx * vx + nz * vz

              tx = rho_vp * vn * nx + rho_vs * (vx - vn * nx)
              tz = rho_vp * vn * nz + rho_vs * (vz - vn * nz)

              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + tx*weight
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(2,iglob) + tz*weight
            else
              ! SH-case
              vy = deltatover2
              ty = rho_vs * vy
              rmass_inverse_elastic(1,iglob) = rmass_inverse_elastic(1,iglob) + ty*weight
              ! ficticous for SH case, but to be save when inverting
              rmass_inverse_elastic(2,iglob) = rmass_inverse_elastic(1,iglob)
            endif
          enddo
        endif  !  end of top absorbing boundary
      endif ! ispec_is_elastic
    enddo

  endif ! any_elastic

  ! acoustic elements
  if (any_acoustic) then

    do ispecabs = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(ispecabs)

      ! Sommerfeld condition if acoustic
      if (ispec_is_acoustic(ispec)) then

        !--- left absorbing boundary
        if (codeabs(IEDGE4,ispecabs)) then
          i = 1
          jbegin = ibegin_edge4(ispecabs)
          jend = iend_edge4(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            rho_vp = rho_vpstore(i,j,ispec)

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)

            weight = jacobian1D * wzgll(j)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) weight = weight * coord(1,iglob)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + deltatover2*weight/rho_vp
          enddo
        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if (codeabs(IEDGE2,ispecabs)) then
          i = NGLLX
          jbegin = ibegin_edge2(ispecabs)
          jend = iend_edge2(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            rho_vp = rho_vpstore(i,j,ispec)

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)

            weight = jacobian1D * wzgll(j)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral
            if (AXISYM) weight = weight * coord(1,iglob)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + deltatover2*weight/rho_vp
          enddo
        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if (codeabs(IEDGE1,ispecabs)) then
          j = 1
          ibegin = ibegin_edge1(ispecabs)
          iend = iend_edge1(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if (codeabs_corner(1,ispecabs)) ibegin = 2
          if (codeabs_corner(2,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            rho_vp = rho_vpstore(i,j,ispec)

            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)

            weight = jacobian1D * wxgll(i)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral,
            if (AXISYM) then
              ! using GLJ quadrature along the radius for elements on the axis
              if (is_on_the_axis(ispec)) then
                if (i == 1) then
                  ! first GLJ point: limit of r/(1+xi) at the axis
                  r_xiplus1 = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                else
                  r_xiplus1 = coord(1,iglob) / (xiglj(i) + ONE)
                endif
                weight = jacobian1D * wxglj(i) * r_xiplus1
              else
                weight = weight * coord(1,iglob)
              endif
            endif

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + deltatover2*weight/rho_vp
          enddo
        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if (codeabs(IEDGE3,ispecabs)) then
          j = NGLLZ
          ibegin = ibegin_edge3(ispecabs)
          iend = iend_edge3(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if (codeabs_corner(3,ispecabs)) ibegin = 2
          if (codeabs_corner(4,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            rho_vp = rho_vpstore(i,j,ispec)

            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)

            weight = jacobian1D * wxgll(i)
            ! axisymmetry: adds cylindrical metric factor r to the boundary integral,
            if (AXISYM) then
              ! using GLJ quadrature along the radius for elements on the axis
              if (is_on_the_axis(ispec)) then
                if (i == 1) then
                  ! first GLJ point: limit of r/(1+xi) at the axis
                  r_xiplus1 = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                else
                  r_xiplus1 = coord(1,iglob) / (xiglj(i) + ONE)
                endif
                weight = jacobian1D * wxglj(i) * r_xiplus1
              else
                weight = weight * coord(1,iglob)
              endif
            endif

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + deltatover2*weight/rho_vp
          enddo
        endif  !  end of top absorbing boundary
      endif ! ispec_is_acoustic
    enddo

  endif ! any_acoustic

  ! electromagnetic elements
  if (any_electromagnetic) then

    do ispecabs = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(ispecabs)

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1

        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)), &
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)

        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

        ! normal pointing left
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        vx = deltatover2
        vz = deltatover2
        vn = nx*vx+nz*vz
        tx = eperml(1) * cpxl *vn*nx
        tz = eperml(2) * cpzl *vn*nz
        weight = jacobian1D * wzgll(j)
        rmass_inverse_electromagnetic(1,iglob) = rmass_inverse_electromagnetic(1,iglob) + tx*weight
        rmass_inverse_electromagnetic(2,iglob) = rmass_inverse_electromagnetic(2,iglob) + tz*weight
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1

        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)), &
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)

        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

        ! normal pointing right
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        vx = deltatover2
        vz = deltatover2
        vn = nx*vx+nz*vz
        tx = eperml(1) * cpxl *vn*nx
        tz = eperml(2) * cpzl *vn*nz
        weight = jacobian1D * wzgll(j)
        rmass_inverse_electromagnetic(1,iglob) = rmass_inverse_electromagnetic(1,iglob) + tx*weight
        rmass_inverse_electromagnetic(2,iglob) = rmass_inverse_electromagnetic(2,iglob) + tz*weight
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
      do i = 1,NGLLX
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1

        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)), &
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)

        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        vx = deltatover2
        vz = deltatover2
        vn = nx*vx+nz*vz
        tx = eperml(1) * cpxl *vn*nx
        tz = eperml(2) * cpzl *vn*nz
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(1,ispecabs) .and. i == 1) .or.  (codeabs_corner(2,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif
        weight = jacobian1D * wxgll(i)
        rmass_inverse_electromagnetic(1,iglob) = rmass_inverse_electromagnetic(1,iglob) + tx*weight
        rmass_inverse_electromagnetic(2,iglob) = rmass_inverse_electromagnetic(2,iglob) + tz*weight
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
      do i = 1,NGLLX
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1

        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)), &
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)

        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        vx = deltatover2
        vz = deltatover2
        vn = nx*vx+nz*vz
        tx = eperml(1) * cpxl *vn*nx
        tz = eperml(2) * cpzl *vn*nz
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(3,ispecabs) .and. i == 1) .or.  (codeabs_corner(4,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif
        weight = jacobian1D * wxgll(i)
        rmass_inverse_electromagnetic(1,iglob) = rmass_inverse_electromagnetic(1,iglob) + tx*weight
        rmass_inverse_electromagnetic(2,iglob) = rmass_inverse_electromagnetic(2,iglob) + tz*weight
      enddo
    endif  !  end of top absorbing boundary

    enddo

  endif ! any_electromagnetic

  end subroutine invert_mass_matrix_init_Stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_mass_matrix_acoustic_gravity()

! adds the linearized free-surface gravity contribution to the acoustic mass matrix

  use constants, only: IMAIN,CUSTOM_REAL,NGLLX,NGLLZ,ONE,ZERO,TINYVAL
  use specfem_par, only: myrank,any_acoustic,ibool,nelem_acoustic_surface,acoustic_surface, &
                         rhostore,jacobian,xix,xiz,gammax,gammaz,wxgll,wzgll, &
                         rmass_inverse_acoustic,ispec_is_acoustic, &
                         AXISYM,is_on_the_axis,coord,wxglj,xiglj
  use shared_parameters, only: GRAVITY,GRAVITY_CONST

  implicit none

  ! local parameters
  integer :: ispec_acoustic_surface,ispec,ixmin,ixmax,izmin,izmax,i,j,iglob
  double precision :: gravityl,rhol,weight,M_g
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,r_xiplus1

  ! checks if anything to do
  if (.not. any_acoustic) return
  if (.not. GRAVITY) return
  if (nelem_acoustic_surface == 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  adding acoustic gravity free surface to mass matrix"
    call flush_IMAIN()
  endif

  ! check non-zero gravity
  gravityl = dabs(GRAVITY_CONST)
  if (gravityl < TINYVAL) then
    write(IMAIN,*) "  warning: using zero gravity"
    call flush_IMAIN()
  endif

  ! loop over acoustic free surface elements
  do ispec_acoustic_surface = 1, nelem_acoustic_surface
    ispec = acoustic_surface(1,ispec_acoustic_surface)

    ! acoustic domains only
    if (.not. ispec_is_acoustic(ispec)) cycle

    ixmin = acoustic_surface(2,ispec_acoustic_surface)
    ixmax = acoustic_surface(3,ispec_acoustic_surface)
    izmin = acoustic_surface(4,ispec_acoustic_surface)
    izmax = acoustic_surface(5,ispec_acoustic_surface)

    ! skip single point nodes
    if (ixmin == ixmax .and. izmin == izmax) cycle

    ! note: linearized dynamic free surface has a boundary term:
    !         int_Gamma chi_ddot w / (rho g) dGamma
    !       which is added as a positive contribution to the acoustic mass matrix:
    !         M = M_V + M_g
    !       for GLL nodes on the free surface.
    !       in a weak form, M_g term is: weight * jacobian1D / (rho g)
    !
    ! for details, see: doc/how_to/README_acoustic_gravity_free_surface.md

    if (izmin == izmax) then
      ! horizontal edge (ITOP or IBOTTOM)
      j = izmin
      do i = 1, NGLLX
        iglob = ibool(i,j,ispec)

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)

        weight = dble(jacobian1D * wxgll(i))
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            if (i == 1) then
              r_xiplus1 = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            else
              r_xiplus1 = coord(1,iglob) / (xiglj(i) + ONE)
            endif
            weight = dble(jacobian1D * wxglj(i) * r_xiplus1)
          else
            weight = weight * coord(1,iglob)
          endif
        endif

        rhol = dble(rhostore(i,j,ispec))
        if (rhol <= 0.d0) call exit_MPI(myrank,'zero density on acoustic free surface')

        ! add M_g term: weight_gll * jacobian1D / (rho g)
        if (gravityl < TINYVAL) then
          ! zero gravity
          M_g = 0.d0
        else
          M_g = weight / (rhol * gravityl)
        endif
        rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + real(M_g, kind=CUSTOM_REAL)
      enddo

    else if (ixmin == ixmax) then
      ! vertical edge (ILEFT or IRIGHT)
      i = ixmin
      do j = 1, NGLLZ
        iglob = ibool(i,j,ispec)

        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)

        weight = dble(jacobian1D * wzgll(j))
        if (AXISYM) weight = weight * coord(1,iglob)

        rhol = dble(rhostore(i,j,ispec))
        if (rhol <= 0.d0) call exit_MPI(myrank,'zero density on acoustic free surface')

        ! add M_g term: weight_gll * jacobian1D / (rho g)
        if (gravityl < TINYVAL) then
          ! zero gravity
          M_g = 0.d0
        else
          M_g = weight / (rhol * gravityl)
        endif
        rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + real(M_g, kind=CUSTOM_REAL)
      enddo

    endif

  enddo

  end subroutine prepare_mass_matrix_acoustic_gravity

!
!-------------------------------------------------------------------------------------------------
!

  subroutine invert_mass_matrix()

! inverts the global mass matrix

  use specfem_par, only: myrank,any_elastic,any_acoustic,any_poroelastic,any_electromagnetic, &
                                rmass_inverse_elastic, &
                                rmass_inverse_acoustic, &
                                rmass_inverse_e1,ATTENUATION_VISCOACOUSTIC, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic,rmass_inverse_electromagnetic
  implicit none
  include 'constants.h'

  if (myrank == 0) then
    write(IMAIN,*) "  inverting mass matrices"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
! (this can happen when some degrees of freedom have been removed from some of the global arrays)
  if (any_elastic) then
    where(rmass_inverse_elastic <= 0._CUSTOM_REAL) rmass_inverse_elastic = 1._CUSTOM_REAL
  endif

  if (any_poroelastic) then
    where(rmass_s_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_s_inverse_poroelastic = 1._CUSTOM_REAL
    where(rmass_w_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_w_inverse_poroelastic = 1._CUSTOM_REAL
  endif

  if (any_acoustic) then
    where(rmass_inverse_acoustic <= 0._CUSTOM_REAL) rmass_inverse_acoustic = 1._CUSTOM_REAL
    if (ATTENUATION_VISCOACOUSTIC) where(abs(rmass_inverse_e1) <= 0._CUSTOM_REAL) rmass_inverse_e1 = 1._CUSTOM_REAL
  endif

  if (any_electromagnetic) then
    where(rmass_inverse_electromagnetic <= 0._CUSTOM_REAL) rmass_inverse_electromagnetic = 1._CUSTOM_REAL
  endif

! compute the inverse of the mass matrix
  if (any_elastic) then
    rmass_inverse_elastic(:,:) = 1._CUSTOM_REAL / rmass_inverse_elastic(:,:)
  endif

  if (any_poroelastic) then
    rmass_s_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_s_inverse_poroelastic(:)
    rmass_w_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_w_inverse_poroelastic(:)
  endif

  if (any_acoustic) then
    rmass_inverse_acoustic(:) = 1._CUSTOM_REAL / rmass_inverse_acoustic(:)
    if (ATTENUATION_VISCOACOUSTIC) rmass_inverse_e1(:,:) = 1._CUSTOM_REAL / rmass_inverse_e1(:,:)
  endif

  if (any_electromagnetic) then
    rmass_inverse_electromagnetic(:,:) = 1._CUSTOM_REAL / rmass_inverse_electromagnetic(:,:)
  endif

  end subroutine invert_mass_matrix

