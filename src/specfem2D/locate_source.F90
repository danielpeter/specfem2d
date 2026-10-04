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

!----
!---- locate_source finds the correct position of the (point force/momen-tensor) source
!----

  subroutine locate_sources()

  use constants, only: NDIM,NGLLX,NGLLZ,IMAIN,HUGEVAL,TINYVAL,NUM_ITER,USE_BEST_LOCATION_FOR_SOURCE, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC,IDOMAIN_ELECTROMAGNETIC

  use specfem_par, only: initialfield, &
                         coord,ibool,nglob,nspec, &
                         xigll,zigll,npgeo, &
                         NPROC,myrank,coorg,knods,NGNOD, &
                         nelem_acoustic_surface,acoustic_surface, &
                         ispec_is_elastic,ispec_is_poroelastic,ispec_is_electromagnetic, &
                         GRAVITY

  use specfem_par, only: NSOURCES,source_type, &
                         x_source,z_source, &
                         ispec_selected_source,islice_selected_source,iglob_source, &
                         xi_source,gamma_source

#ifndef WITH_MPI
  use constants, only: OUTPUT_FILES
#endif

  implicit none
  ! local parameters
  integer :: i_source
  double precision :: x_target,z_target
  integer :: ispec_acoustic_surface,ispec
  integer :: ixmin, ixmax, izmin, izmax

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating sources'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! locates sources
  do i_source = 1,NSOURCES
    ! double-checks source type for simulations w/out initialfield
    if (.not. initialfield) then
      if (source_type(i_source) /= 1 .and. source_type(i_source) /= 2) then
        print *,'Error: Source ',i_source,' has invalid source type ',source_type(i_source)
        print *,'       must be either 1 (force) or 2 (moment tensor)'
        call exit_MPI(myrank,'incorrect source type')
      endif
    endif

    ! target position
    x_target = x_source(i_source)
    z_target = z_source(i_source)

    ! gets source positioning
    ! note: iglob_source is not really needed for moment-tensor sources, but left as argument since it's already allocated
    !       and might help for future routines...
    call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                       x_target,z_target, &
                       ispec_selected_source(i_source),islice_selected_source(i_source), &
                       NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,NGNOD,npgeo, &
                       iglob_source(i_source),source_type(i_source),i_source)

    ! checks force source positions for acoustic domains
    if (source_type(i_source) == 1 .and. .not. GRAVITY) then
      ! only slice which holds source
      if (myrank == islice_selected_source(i_source)) then
        ! checks that acoustic source is not exactly on the free surface because pressure is zero there
        do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          ! checks if source element
          if (ispec == ispec_selected_source(i_source)) then
            if (.not. ispec_is_elastic(ispec) .and. &
                .not. ispec_is_poroelastic(ispec) .and. &
                .not. ispec_is_electromagnetic(ispec)) then
              ! acoustic element
              ixmin = acoustic_surface(2,ispec_acoustic_surface)
              ixmax = acoustic_surface(3,ispec_acoustic_surface)
              izmin = acoustic_surface(4,ispec_acoustic_surface)
              izmax = acoustic_surface(5,ispec_acoustic_surface)
              if ((izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == NGLLX .and. &
                  gamma_source(i_source) < -0.99d0) .or. &
                  (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == NGLLX .and. &
                  gamma_source(i_source) > 0.99d0) .or. &
                  (izmin == 1 .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
                  xi_source(i_source) < -0.99d0) .or. &
                  (izmin == 1 .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                  xi_source(i_source) > 0.99d0) .or. &
                  (izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == 1 .and. &
                  gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) < -0.99d0) .or. &
                  (izmin == 1 .and. izmax == 1 .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                  gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) > 0.99d0) .or. &
                  (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
                  gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) < -0.99d0) .or. &
                  (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                  gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) > 0.99d0)) then
                call exit_MPI(myrank,'an acoustic source cannot be located exactly '// &
                              'on the free surface because pressure is zero there')
              endif
            endif
          endif
        enddo
      endif
    endif

  enddo  ! NSOURCES

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == islice_selected_rec(irec)) then" to display the right sources
!! DK DK and receivers carried by each mesh slice, and not fictitious values coming from other slices
#ifndef WITH_MPI
  if (myrank == 0) then
     ! write actual source locations to file
     ! note that these may differ from input values, especially if source_surf = .true. in SOURCE
     ! note that the exact source locations are determined from (ispec,xi,gamma) values
     open(unit=14,file=trim(OUTPUT_FILES)//'for_information_SOURCE_actually_used',status='unknown')
     do i_source = 1,NSOURCES
        write(14,*) x_source(i_source), z_source(i_source)
     enddo
     close(14)
  endif
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine locate_sources

!
!----------------------------------------------------------------------------------
!

  subroutine locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_target,z_target, &
                           ispec_selected_source,islice_selected_source, &
                           NPROC,myrank, &
                           xi_source,gamma_source,coorg,knods,NGNOD,npgeo, &
                           iglob_source,source_type,i_source)

  use constants, only: NDIM,NGLLX,NGLLZ,IMAIN,HUGEVAL,TINYVAL,NUM_ITER,USE_BEST_LOCATION_FOR_SOURCE, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC,IDOMAIN_ELECTROMAGNETIC

  use specfem_par, only: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic,ispec_is_electromagnetic, &
    SOURCE_IS_MOVING

  use specfem_par, only: AXISYM,is_on_the_axis,xiglj

  implicit none

  integer,intent(in) :: nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec),intent(in) :: ibool

  ! array containing coordinates of the points
  double precision,intent(in) :: coord(NDIM,nglob)

  ! Gauss-Lobatto-Legendre points of integration
  double precision,intent(in) :: xigll(NGLLX)
  double precision,intent(in) :: zigll(NGLLZ)

  double precision,intent(in) :: x_target,z_target

  integer,intent(in) :: NGNOD,npgeo

  integer,intent(in) :: knods(NGNOD,nspec)
  double precision,intent(in) :: coorg(NDIM,npgeo)

  ! source information
  integer,intent(out) :: ispec_selected_source,islice_selected_source,iglob_source

  integer,intent(in)  :: NPROC, myrank
  double precision,intent(out) :: xi_source,gamma_source

  integer,intent(in) :: source_type,i_source

  !also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
  logical, dimension(nglob) :: flag_topological
  integer :: number_of_mesh_elements_for_the_initial_guess
  integer, dimension(:), allocatable :: array_of_all_elements_of_ispec_selected_source

  ! local parameters
  integer :: i,j,ispec,iglob
  integer :: iter_loop,ix_initial_guess,iz_initial_guess,number_of_iterations
  integer :: ispec_found,islice_found,iglob_found,idomain_found
  integer :: imin,imax,jmin,jmax

  double precision :: dist_squared
  double precision :: xi,gamma,dx,dz,dxi,dgamma
  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian
  double precision :: distmin_squared,dist_glob_squared
  double precision :: final_distance,final_distance_this_element
  double precision :: xi_found,gamma_found

  integer :: is_proc_source
  integer, dimension(1:NPROC)  :: allgather_is_proc_source
  integer, dimension(1)  :: locate_is_proc_source
  logical :: is_force_source

  ! initializes slice and element which hold source
  is_proc_source = 0
  ispec_found = 0
  islice_found = 0
  iglob_found = 0 ! closest global point

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! search range
  imin = 1
  imax = NGLLX
  jmin = 1
  jmax = NGLLZ

  ! loops over all elements
  do ispec = 1,nspec
    ! loop only on points inside the element
    ! exclude edges to ensure this point is not shared with other elements
    do j = jmin,jmax
      do i = imin,imax
        iglob = ibool(i,j,ispec)

        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        dist_squared = (x_target - dble(coord(1,iglob)))**2 + (z_target - dble(coord(2,iglob)))**2

        ! keep this point if it is closer to the source
        if (dist_squared < distmin_squared) then
          iglob_found = iglob
          distmin_squared = dist_squared
          ispec_found = ispec
          ix_initial_guess = i
          iz_initial_guess = j
        endif
      enddo
    enddo
  enddo

  ! also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
  flag_topological(:) = .false.

  ! mark the four corners of the initial guess element
  flag_topological(ibool(1,1,ispec_found)) = .true.
  flag_topological(ibool(NGLLX,1,ispec_found)) = .true.
  flag_topological(ibool(NGLLX,NGLLZ,ispec_found)) = .true.
  flag_topological(ibool(1,NGLLZ,ispec_found)) = .true.

  ! loop on all the elements to count how many are shared with the initial guess
  number_of_mesh_elements_for_the_initial_guess = 1
  do ispec = 1,nspec
    if (ispec == ispec_found) cycle  ! skip the initial guess element

    ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
    do j = 1,NGLLZ,NGLLZ-1
      do i = 1,NGLLX,NGLLX-1
        if (flag_topological(ibool(i,j,ispec))) then
          ! this element is in contact with the initial guess
          number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
          ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
          goto 700  ! Exit the loop over the four corners: we know this element is adjacent
        endif
      enddo
    enddo

    700 continue
  enddo

  ! now that we know the number of adjacent elements, we can allocate the list of elements and create it
  allocate(array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess))
  array_of_all_elements_of_ispec_selected_source(:) = 0

  ! first store the initial guess itself
  number_of_mesh_elements_for_the_initial_guess = 1
  array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess) = ispec_found

  ! then store all the others
  do ispec = 1,nspec
    if (ispec == ispec_found) cycle

    ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
    do j = 1,NGLLZ,NGLLZ-1
      do i = 1,NGLLX,NGLLX-1
        if (flag_topological(ibool(i,j,ispec))) then
          ! this element is in contact with the initial guess
          number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
          array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess) = ispec
          ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
          goto 800
        endif
      enddo
    enddo

    800 continue
  enddo

  ! initialize best position found
  ispec_found = 0
  xi_found = 0.d0
  gamma_found = 0.d0
  final_distance = HUGEVAL
  idomain_found = 0

  do i = 1,number_of_mesh_elements_for_the_initial_guess

    ! set initial guess in the middle of the element, since we computed the true one only for the true initial guess
    ! the nonlinear process below will converge anyway
    if (i > 1) then
      ix_initial_guess = int(NGLLX / 2.0)
      iz_initial_guess = int(NGLLZ / 2.0)
    endif

    ispec = array_of_all_elements_of_ispec_selected_source(i)

    ! ****************************************
    ! find the best (xi,gamma) for each source
    ! ****************************************

    ! use initial guess in xi and gamma
    xi = xigll(ix_initial_guess)
    gamma = zigll(iz_initial_guess)
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        xi = xiglj(ix_initial_guess)
      endif
    endif

    ! iterate to solve the nonlinear system
    if (USE_BEST_LOCATION_FOR_SOURCE) then
      number_of_iterations = NUM_ITER
    else
      number_of_iterations = 0 ! this means that the loop below will not be executed, i.e. we will not iterate
    endif

    do iter_loop = 1,number_of_iterations

      ! recompute jacobian for the new point
      call recompute_jacobian_with_negative_stop(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                                                 coorg,knods,ispec,NGNOD,nspec,npgeo,.true.)

      ! compute distance to target location
      dx = - (x - x_target)
      dz = - (z - z_target)

      ! compute increments
      dxi  = xix*dx + xiz*dz
      dgamma = gammax*dx + gammaz*dz

      ! update values
      xi = xi + dxi
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a source outside the mesh for instance)
      ! we can go slightly outside the [1,1] segment since with finite elements
      ! the polynomial solution is defined everywhere
      ! this can be useful for convergence of itertive scheme with distorted elements
      if (xi > 1.01d0) xi = 1.01d0
      if (xi < -1.01d0) xi = -1.01d0
      if (gamma > 1.01d0) gamma = 1.01d0
      if (gamma < -1.01d0) gamma = -1.01d0

    ! end of nonlinear iterations
    enddo

    ! compute final coordinates of point found
    call recompute_jacobian_with_negative_stop(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                                               coorg,knods,ispec,NGNOD,nspec,npgeo,.true.)

    ! compute final distance between asked and found
    final_distance_this_element = sqrt((x_target - x)**2 + (z_target - z)**2)

    ! if we have found an element that gives a shorter distance
    if (final_distance_this_element < final_distance) then
      ! store element number found
      ispec_found = ispec

      ! store xi,gamma of point found
      xi_found = xi
      gamma_found = gamma

      ! store final distance between asked and found
      final_distance = final_distance_this_element

      ! determines domain for outputting element type
      if (ispec_is_acoustic(ispec)) then
        idomain_found = IDOMAIN_ACOUSTIC
      else if (ispec_is_elastic(ispec)) then
        idomain_found = IDOMAIN_ELASTIC
      else if (ispec_is_poroelastic(ispec)) then
        idomain_found = IDOMAIN_POROELASTIC
      else if (ispec_is_electromagnetic(ispec)) then
        idomain_found = IDOMAIN_ELECTROMAGNETIC
      else
        call stop_the_code('Invalid element type in locating source found!')
      endif
    endif
  enddo

  ! global minimum distance computed over all processes
  call min_all_all_dp(final_distance, dist_glob_squared)

  ! check if this process contains the source
  if (abs(sqrt(dist_glob_squared) - sqrt(final_distance)) < TINYVAL ) is_proc_source = 1

  ! main collects info
  call gather_all_singlei(is_proc_source,allgather_is_proc_source,NPROC)

  if (myrank == 0) then
    ! select slice with maximum rank which contains source
    locate_is_proc_source = maxloc(allgather_is_proc_source) - 1
    islice_found = locate_is_proc_source(1)
  endif

  ! selects slice which holds source
  call bcast_all_singlei(islice_found)

#ifdef WITH_MPI
  ! for MPI version, gather information from all the nodes
  if (islice_found /= 0) then
    ! source is in another slice than the main process
    if (myrank == islice_found) then
      ! send information from process holding source
      call send_singlei(ispec_found,0,0)
      call send_singlei(idomain_found,0,1)
      call send_singledp(xi_found,0,2)
      call send_singledp(gamma_found,0,3)
      call send_singledp(final_distance,0,4)
      call send_singlei(iglob_found,0,5)
    else if (myrank == 0) then
      ! main collects
      call recv_singlei(ispec_found,islice_found,0)
      call recv_singlei(idomain_found,islice_found,1)
      call recv_singledp(xi_found,islice_found,2)
      call recv_singledp(gamma_found,islice_found,3)
      call recv_singledp(final_distance,islice_found,4)
      call recv_singlei(iglob_found,islice_found,5)
    endif
  endif
#endif

  ! store found position
  ! note: iglob_source is not really needed for moment-tensor sources, but left as argument since it's already allocated
  !       and might help for future routines...
  islice_selected_source = islice_found
  ispec_selected_source = ispec_found
  iglob_source = iglob_found
  xi_source = xi_found
  gamma_source = gamma_found

  ! sets source type flag (must be either 1 == force or 2 == moment tensor
  if (source_type == 1) then
    is_force_source = .true.
  else
    is_force_source = .false.
  endif

  ! user output
  if (myrank == 0) then
    if ((.not. is_force_source) .or. (is_force_source .and. (.not. SOURCE_IS_MOVING))) then
      ! source info
      write(IMAIN,*) 'Source # ',i_source
      write(IMAIN,*) '  source located in slice ',islice_selected_source
      write(IMAIN,*) '                 in element ',ispec_selected_source
      if (idomain_found == IDOMAIN_ACOUSTIC) then
        write(IMAIN,*) '                 in acoustic domain'
      else if (idomain_found == IDOMAIN_ELASTIC) then
        write(IMAIN,*) '                 in elastic domain'
      else if (idomain_found == IDOMAIN_POROELASTIC) then
        write(IMAIN,*) '                 in poroelastic domain'
      else if (idomain_found == IDOMAIN_ELECTROMAGNETIC) then
        write(IMAIN,*) '                 in electromagnetic domain'
      else
        write(IMAIN,*) '                 in unknown domain'
      endif
      write(IMAIN,*)
      if (is_force_source) then
        write(IMAIN,*) '  using force point source:'
      else
        write(IMAIN,*) '  using moment-tensor source:'
      endif
      write(IMAIN,*) '    xi coordinate of source in that element   :',xi_source
      write(IMAIN,*) '    gamma coordinate of source in that element:',gamma_source
      write(IMAIN,*)
      write(IMAIN,*) '    original x            : ',sngl(x_target)
      write(IMAIN,*) '    original z            : ',sngl(z_target)
      write(IMAIN,*) '    closest estimate found: ',sngl(final_distance),' m away'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! check position
    if (final_distance == HUGEVAL) call exit_MPI(myrank,'Error locating source')
  endif

  ! free temporary array
  deallocate(array_of_all_elements_of_ispec_selected_source)

  end subroutine locate_source

