PROGRAM POD
  USE def_type_mesh
  USE initialization
  USE my_util
  USE input_data
  USE symmetric_field
  ! USE arpack_mhd
  USE fourier_to_real_for_vtu
  USE user_data
  USE post_processing_debug
  USE sft_parallele
  USE boundary
  USE vtk_viz
  USE restart
  USE fem_tn_axi
  USE verbose
! CN 15/03/2024
  USE basis_change
#include "petsc/finclude/petsc.h"
  USE petsc
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn
  TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un
  !===Maxwell fields==============================================================
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: Hn, Bn, phin, vel
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: sigma_field, mu_H_field
  !===Temperature field===========================================================
  TYPE(mesh_type), POINTER                        :: temp_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: temperature
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: vol_heat_capacity_field
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: temperature_diffusivity_field
  !===Concentration field===========================================================
  TYPE(mesh_type), POINTER                        :: conc_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: concentration
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: concentration_diffusivity_field
  !===Level_set===================================================================
  REAL(KIND = 8), POINTER, DIMENSION(:, :, :, :)  :: level_set
  !===Density=====================================================================
  REAL(KIND = 8), POINTER, DIMENSION(:, :, :)     :: density
  !===LES=========================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: visc_LES
  !===Fourier modes===============================================================
  INTEGER                                         :: m_max_c
  INTEGER,      POINTER,      DIMENSION(:)        :: list_mode
  !===Time iterations=============================================================
  REAL(KIND=8)                                    :: time
  INTEGER                                         :: it
  !===Timing======================================================================
  REAL(KIND=8)                                    :: tps, tploc, tploc_max=0.d0

  !===Declare PETSC===============================================================
  PetscErrorCode :: ierr
  PetscMPIInt    :: rank
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d, comm_one_d_ns, comm_one_d_temp
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_conc

  !===pp_mesh w_c creation========================================================
  REAL(KIND=8), DIMENSION(:,:), POINTER, SAVE :: aij_p1p2, aij_p2p3
  REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
  LOGICAL :: once_p1p2=.TRUE., once_p2p3=.TRUE.
  INTEGER :: l, n

  !===Start PETSC and MPI (mandatory)=============================================
  CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

  !===User reads his/her own data=================================================
  CALL read_user_data('data')

  !===Initialize SFEMANS (mandatory)==============================================
  CALL initial(vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, conc_mesh,&
      interface_H_phi, interface_H_mu, list_mode, &
      un, pn, Hn, Bn, phin, vel, &
      vol_heat_capacity_field, temperature_diffusivity_field, &
      concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
      comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
      concentration, level_set, density, &
      der_un, visc_LES)

  
!================================================================================
!                                START TEST
!================================================================================
 
CALL run_pod_galerkin()

IF (rank == 0) THEN
  WRITE(*,*) "============================"
  WRITE(*,*) "CODE RAN SUCCESSFULLY"
  WRITE(*,*) "============================"
END IF
!================================================================================
!                                END TEST
!================================================================================
    
  CALL error_petsc('End of pod galerking')
  CONTAINS

END PROGRAM POD
  