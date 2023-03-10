! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts
!
!> \brief MPAS barotropic ocean SSPRK3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!
!-----------------------------------------------------------------------

module ocn_time_integration_ssprk3 

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_mesh
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_equation_of_state
   use ocn_time_average_coupled
   use ocn_time_varying_forcing 

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   public :: ocn_time_integrator_ssprk3

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_ssprk3
!
!> \brief MPAS barotropic ocean LTS time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!>  This routine integrates one timestep (dt) using an LTS time integrator
!>  with a splitting of the fast and slow tendency terms
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_ssprk3(domain,dt)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    ! a local time stepping scheme with spltting of the fast and slow tendencies
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !-----------------------------------------------------------------
    ! Input variables
    !-----------------------------------------------------------------

    real (kind=RKIND), intent(in) :: &
       dt      !< [in] time step (sec) to move forward

    !-----------------------------------------------------------------
    ! Input/output variables
    !-----------------------------------------------------------------

    type (domain_type), intent(inout) :: &
       domain  !< [inout] model state to advance forward

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    integer :: &
       iCell, iEdge, iRegion, k, ic, ie, im, err, & !iterators
       M, & ! M = dtCoarse / dtFine
       nRegions !number of interface regions (two)

    type (block_type), pointer :: &
       block ! structure with subdomain data

    type (mpas_pool_type), pointer :: &
       tendPool, & ! structure holding tendencies
       statePool, & ! structure holding state variables
       meshPool, & ! structure holding mesh variables
       forcingPool, & ! structure holding forcing variables
       scratchPool, & ! structure holding temporary variables
       tracersPool, &  ! structure holding tracers variables
       verticalMeshPool

    ! LTS Pools
    type (mpas_pool_type), pointer :: &
       LTSPool  ! structure holding LTS variables

    ! Tend Array Pointers
    real (kind=RKIND), dimension(:,:), pointer ::  &
       normalVelocityTend, & ! normal velocity fast tendency
       layerThicknessTend  ! layer thickness tendency

    ! State Array Pointers
    real (kind=RKIND), dimension(:), pointer :: &
       sshCur, &
       sshNew
    real (kind=RKIND), dimension(:,:), pointer :: &
       normalVelocityCur, & ! normal velocity at time n
       normalVelocityNew, & ! normal velocity at time n+1
       layerThicknessCur, & ! layer thickness at time n
       layerThicknessNew  ! layer thickness at time n+1

    ! LTS objects
    integer, dimension(:,:), pointer :: &
       nCellsInLTSRegion, & ! number of cells in a given LTS region
       nEdgesInLTSRegion ! number of edges in a given LTS region
    integer, dimension(:,:,:), pointer :: & 
       cellsInLTSRegion, & ! list of cells in a given LTS region
       edgesInLTSRegion ! list of edges in a given LTS region

    integer :: ssprk_step

    real (kind=RKIND), dimension(3) :: weights_old, weights_new, weights_tend

    call mpas_timer_start("lts time-step prep")

    err = 0

    nRegions = 2

    block => domain % blocklist

    !--- Retrieve model state, pools
    !--- No longer support sub-blocks so this retrieves the only block

    call mpas_pool_get_subpool(block%structs, 'mesh', & 
                                               meshPool)
    call mpas_pool_get_subpool(block%structs, 'state', &
                                               statePool)
    call mpas_pool_get_subpool(block%structs, 'forcing', &
                                               forcingPool)
    call mpas_pool_get_subpool(block%structs, 'LTS', &
                                               LTSPool)
    call mpas_pool_get_subpool(block%structs, 'tend', &
                                               tendPool)
    call mpas_pool_get_subpool(block%structs, 'scratch', &
                                               scratchPool)
    call mpas_pool_get_subpool(statePool, 'tracers', &
                                               tracersPool)
    call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

    !--- Retrieve state variables at necessary time levels

    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocityCur, 1)
    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocityNew, 2)

    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessCur, 1)
    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessNew, 2)

    call mpas_pool_get_array(statePool, 'ssh', &
                                         sshCur, 1)
    call mpas_pool_get_array(statePool, 'ssh', &
                                         sshNew, 2)

    !--- Retrieve tendency variables

    call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                        normalVelocityTend)
    call mpas_pool_get_array(tendPool, 'layerThickness', &
                                        layerThicknessTend)

    !--- Retrieve LTS arrays
    call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', &
                                       cellsInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', &
                                       nCellsInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', & 
                                       edgesInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', &
                                       nEdgesInLTSRegion)


    !--- Init variables
    do iEdge = 1, nEdgesAll
       do k = 1, maxLevelEdgeTop(iEdge)
          normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
       end do
    end do

    do iCell = 1, nCellsAll
       sshNew(iCell) = sshCur(iCell)
       do k = 1, maxLevelCell(iCell)
          layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
       end do
    end do

    call mpas_timer_stop("lts time-step prep")

    call mpas_timer_start("lts main loop")

    weights_old(1) = 1.0_RKIND
    weights_old(2) = 0.75_RKIND !use this for SSPRK3
    !weights_old(2) = 0.5   !use this for SSPRK2
    weights_old(3) = 1.0_RKIND/3.0_RKIND

    weights_new(1) = 0.0_RKIND
    weights_new(2) = 0.25_RKIND !use this for SSPRK3
    !weights_new(2) = 0.5   !use this for SSPRK2
    weights_new(3) = 2.0_RKIND/3.0_RKIND

    weights_tend(1) = dt * 1.0_RKIND
    weights_tend(2) = dt * 0.25_RKIND !use this for SSPRK3
    !weights_tend(2) = dt * 0.5   !use this for SSPRK2
    weights_tend(3) = dt * 2.0_RKIND/3.0_RKIND


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BEGIN SSPRK3 SCHEME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do ssprk_step = 1, 3

      if (config_disable_thick_all_tend .and. config_disable_vel_all_tend .and. config_disable_tr_all_tend) then
          exit ! don't compute in loop meant to update velocity, thickness, and tracers
      end if

      !if(ssprk_step == 1) then

      ! do ic = 1, nCellsAll
      !      do k = 1, maxLevelCell(ic)
      !     print*,  layerThicknessCur(k,ic)
      !     end do
      ! end do
      !print*,'this is the tendency call --------------------------------'
      !    end if

      call ocn_vert_transport_velocity_top(verticalMeshPool, layerThicknessNew, layerThickEdgeFlux, normalVelocityNew, sshNew, dt, vertAleTransportTop, err)
      call ocn_tend_vel(domain,tendPool, statePool, forcingPool, 2, domain % dminfo, dt)
      call ocn_vert_transport_velocity_top(verticalMeshPool, layerThicknessNew, layerThickEdgeFlux, normalTransportVelocity, sshNew, dt, vertAleTransportTop, err)
      call ocn_tend_thick(tendPool, forcingPool)

      call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
      call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

      do iEdge = 1, nEdgesAll
             normalVelocityNew(:,iEdge) = weights_old(ssprk_step) * normalVelocityCur(:,iEdge) &
                                        + weights_new(ssprk_step) * normalVelocityNew(:,iEdge) &
                                        + weights_tend(ssprk_step) * normalVelocityTend(:,iEdge)
      end do

      do iCell = 1, nCellsAll
         do k = 1, maxLevelCell(iCell)
            layerThicknessNew(k,iCell) = weights_old(ssprk_step) * layerThicknessCur(k,iCell) &
                                         + weights_new(ssprk_step) * layerThicknessNew(k,iCell) &
                                         + weights_tend(ssprk_step) * layerThicknessTend(k,iCell)
         end do
      end do

    !  if(ssprk_step == 1) then

    ! interface layers
    !     do iRegion =1,nRegions
    !      do ic = 1, nCellsInLTSRegion(iRegion,2)
    !      iCell = cellsInLTSRegion(iRegion,2,ic)
    !       do k = 1, maxLevelCell(iCell)
    !        print*, layerThicknessNew(k,iCell)
    !        end do
    !       end do
    !      end do     

    ! do iRegion =1,nRegions
    !    do ie = 1, nEdgesInLTSRegion(iRegion,2)
    !     iEdge = edgesInLTSRegion(iRegion,2,ie)
    !     do k = 1, nVertLevels
    !       print*, normalVelocityNew(k,iEdge)
    !     end do
    ! end do
    ! end do

    ! coarse
    !    do ic = 1, nCellsInLTSRegion(2,1)
    !     iCell = cellsInLTSRegion(2,1,ic)
    !     do k = 1, maxLevelCell(iCell)
    !     print*, layerThicknessNew(k,iCell)
    !     end do
    !    end do

    !do ie = 1, nEdgesInLTSRegion(2,1)
    !   iEdge = edgesInLTSRegion(2,1,ie)
    !   do k = 1, nVertLevels
    !       print*, normalVelocityNew(k,iEdge)
    !     end do
    ! end do

        ! fine
    !    do ic = 1, nCellsInLTSRegion(1,1)
    !     iCell = cellsInLTSRegion(1,1,ic)
    !     do k = 1, maxLevelCell(iCell)
    !      print*, layerThicknessNew(k,iCell)
    !     end do
    !    end do

    !    do ic = 1, nCellsInLTSRegion(1,3)
    !     iCell = cellsInLTSRegion(1,3,ic)
    !     do k = 1, maxLevelCell(iCell)
    !      print*, layerThicknessNew(k,iCell)
    !     end do
    !    end do

    !do ie = 1, nEdgesInLTSRegion(1,1)
    !   iEdge = edgesInLTSRegion(1,1,ie)
    !   do k = 1, nVertLevels
           !print*, normalVelocityNew(k,iEdge)
           !print*, normalVelocityCur(k,iEdge)
    !       if(ie==1) then
    !       print*, normalVelocityTend(k,iEdge)
    !       print*, iEdge
    !       print*, k
    !       end if
    !     end do
    ! end do

    !do ie = 1, nEdgesInLTSRegion(1,3)
    !   iEdge = edgesInLTSRegion(1,3,ie)
    !   do k = 1, nVertLevels
           !print*, normalVelocityNew(k,iEdge)
           !print*, normalVelocityCur(k,iEdge)
           !print*, normalVelocityTend(k,iEdge)
    !     end do
    ! end do

    !  end if    

      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool,  scratchPool, tracersPool, 2)

      do iEdge = 1, nEdgesAll
         normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
      end do

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END SSPRK3 SCHEME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpas_timer_stop("lts main loop")
      
    !
    !  A little clean up at the end
    !

    call mpas_timer_start("lts cleanup phase")

    if (config_prescribe_velocity) then
       do iEdge = 1, nEdgesAll
          normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
       end do
    end if

    if (config_prescribe_thickness) then
       do iCell = 1, nCellsAll
          layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
       end do
    end if

    do iEdge = 1, nEdgesAll
       normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
    end do

    call mpas_reconstruct(meshPool,  normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

    call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)
    call mpas_threading_barrier()

    do iCell = 1, nCellsAll
       surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
       surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

       SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
       SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
    end do

    call ocn_time_average_coupled_accumulate(statePool, forcingPool, 2)

    call mpas_timer_stop("lts cleanup phase")


  end subroutine ocn_time_integrator_ssprk3!}}}

end module ocn_time_integration_ssprk3
