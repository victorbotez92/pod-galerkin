MODULE pod_galerkin
   USE def_type_mesh
   USE symmetric_field
   USE input_data
   USE user_data
#include "petsc/finclude/petsc.h"
   USE petsc
   IMPLICIT NONE
   ! PUBLIC :: compute_pod_galerking


   TYPE dyn_real_ksp
      KSP, DIMENSION(:), POINTER :: DKSP ! POINTER attribute
   END TYPE dyn_real_ksp     ! instead of ALLOCATABLE

   TYPE dyn_real_mat
      Mat, DIMENSION(:), POINTER :: DMat ! POINTER attribute
   END TYPE dyn_real_mat     ! instead of ALLOCATABLE



    
CONTAINS

! SUBROUTINE Duchon_Robert_MHD(communicator_ns, communicator_mhd, vv_3_LA, vv_1_LA, LA_H_3, LA_H_1, vvrtz_per, vvz_per, &
!                                     pp_mesh, vv_mesh, H_mesh, jj_v_to_H, un, Bn, Hn, sigma, mu, list_mode, number_coeff_DR, it, time, interface_H_mu)
   !------------------------------------------------------------------------------------------------
   !-------------VARIABLES AND MODULES--------------------------------------------------------------
   !------------------------------------------------------------------------------------------------

! ------MODULES----------------------------------------------------------------------------
!-------
   SUBROUTINE gauss_to_nodes(mesh, LA_1, one, zero, V_gauss, V_out, cc_ksp, cb_1, cb_2, x_1, x_1_ghost, my_par)
      USE solve_petsc
      USE fem_rhs_axi
      USE st_matrix

      TYPE(mesh_type), INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), INTENT(IN) :: LA_1
      REAL(KIND = 8), DIMENSION(:), POINTER, INTENT(IN) :: one, zero
      REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V_gauss
      KSP, INTENT(IN) :: cc_ksp
      Vec, INTENT(IN) :: cb_1, cb_2, x_1, x_1_ghost
      TYPE(solver_param), INTENT(IN) :: my_par

      REAL(KIND = 8), DIMENSION(mesh%np, SIZE(V_gauss, 2), SIZE(V_gauss, 3)) :: V_out

      PetscErrorCode :: ierr
      INTEGER :: i, k

      DO i = 1, SIZE(V_gauss, 3)
         DO k = 1, SIZE(V_gauss, 2) / 2

            !===Compute rhs
            CALL qs_00_gauss(mesh, LA_1, one, zero, V_gauss(:, 2 * k - 1, i), cb_1)
            CALL qs_00_gauss(mesh, LA_1, one, zero, V_gauss(:, 2 * k, i), cb_2)

            !===RHS periodicity !todo add periodicity
            !IF (vvz_per%n_bord/=0) THEN
            !    CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, cb_DR_1, vv_1_LA)
            !    CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, cb_DR_2, vv_1_LA)
            !END IF

            !===Solve gauss to node equation
            !Solve system cc_c
            CALL solver(cc_ksp, cb_1, x_1, reinit = .FALSE., verbose = my_par%verbose)
            CALL VecGhostUpdateBegin(x_1, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL VecGhostUpdateEnd(x_1, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL extract(x_1_ghost, 1, 1, LA_1, V_out(:, 2 * k - 1, i))

            !Solve system cc_s
            CALL solver(cc_ksp, cb_2, x_1, reinit = .FALSE., verbose = my_par%verbose)
            CALL VecGhostUpdateBegin(x_1, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL VecGhostUpdateEnd(x_1, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL extract(x_1_ghost, 1, 1, LA_1, V_out(:, 2 * k, i))
         END DO
      END DO

   END SUBROUTINE gauss_to_nodes

   SUBROUTINE compute_L_ij_and_Q_ijk(communicator_ns, communicator_mhd, vv_1_LA, vv_3_LA, LA_H_1, LA_H_3, vvrtz_per, vvz_per, &
                                     pp_mesh, vv_mesh, H_mesh, jj_v_to_H, list_mode)
   
      USE chaine_caractere
      USE fem_M_axi
      USE solve_petsc
      USE st_matrix
      USE fem_rhs_axi
      USE periodic
      USE Gauss_points
      USE sft_parallele
      USE boundary
      USE Dir_nodes_petsc
      !  USE user_mean_function
      USE pdf
      USE fft_DR
      !  USE restart_DR
      USE fourier_to_real_for_vtu
      !  USE projections
      USE tn_axi
      !  USE dt_mu


      

      !----------VARIABLES-----------------------------------------------------------------------------


      !  REAL(KIND = 8), INTENT(IN) :: time
      !  INTEGER, INTENT(IN) :: it
      INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
      !  INTEGER, INTENT(IN) :: number_coeff_DR
      TYPE(mesh_type), INTENT(IN) :: pp_mesh, vv_mesh, H_mesh
      INTEGER, DIMENSION(:), INTENT(IN) :: jj_v_to_H
      TYPE(periodic_type), INTENT(IN) :: vvrtz_per, vvz_per

      !  REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)), INTENT(IN) :: un
      !  REAL(KIND = 8), DIMENSION(H_mesh%np, 6, SIZE(list_mode)), INTENT(IN) :: Hn
      !  REAL(KIND = 8), DIMENSION(H_mesh%np, 6, SIZE(list_mode)), INTENT(IN) :: Bn
      !  REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: sigma, mu
      TYPE(interface_type), TARGET :: interface_H_mu

      REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE, SAVE :: vv_rr_gauss
      REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE, SAVE :: H_rr_gauss

      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: sort_vv_rr_gauss
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: sort_H_rr_gauss

      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: div_u_cube, grad_u_square
      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: Bn_vv, u_cross_b
      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: one_over_mu_H
      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: rot_B_H, j_cross_b_H, j_H, grad_one_over_mu

      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: u_gauss
      REAL(KIND = 8), DIMENSION(:, :, :, :), ALLOCATABLE, SAVE :: Gu_gauss

      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: Vs
      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 2) :: mus
      REAL(KIND = 8), DIMENSION(H_mesh%gauss%n_w, 6) :: Hs, Bs

      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, 2, SIZE(list_mode)) :: scalar_gauss, scalar_gauss2, scalar_gauss3
      REAL(KIND = 8), DIMENSION(3, vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)) :: Gu_reg_gauss

      REAL(KIND = 8), DIMENSION(H_mesh%gauss%l_G * H_mesh%me, 6, SIZE(list_mode)) :: grad_one_over_mu_gauss, j_gauss
      REAL(KIND = 8), DIMENSION(H_mesh%gauss%l_G * H_mesh%me, 6, SIZE(list_mode)) :: rot_B_gauss

      REAL(KIND = 8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode)) :: scalar_nodes_vv, scalar_nodes_vv2
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode)) :: scalar_nodes_vv3
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)) :: vector_nodes_vv, vector_nodes_vv2
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)) :: vector_nodes_vv3
      REAL(KIND = 8), DIMENSION(H_mesh%np, 2, SIZE(list_mode)) :: scalar_nodes_H, lambda_H, scalar_nodes_H2
      REAL(KIND = 8), DIMENSION(H_mesh%np, 6, SIZE(list_mode)) :: vector_nodes_H, vector_nodes_H2, vector_nodes_H3

      INTEGER :: n, i, k, kp, mode, code, index, m, l, nu_mat
      INTEGER, DIMENSION(2) :: ns
      INTEGER :: n1, n2, n3, n123, cas

      !---------petsc & mpi
      PetscErrorCode                   :: ierr
      MPI_Comm, DIMENSION(:), POINTER :: communicator_mhd, communicator_ns
      TYPE(petsc_csr_LA), INTENT(IN) :: vv_1_LA, vv_3_LA, LA_H_1, LA_H_3
      INTEGER, POINTER, DIMENSION(:) :: vv_1_ifrom, vv_3_ifrom, H_1_ifrom, H_3_ifrom
      Vec,                                    SAVE :: cb_DR_1, cb_DR_2, vx_1, vx_1_ghost
      Vec,                                    SAVE :: cb_DR_1_H, cb_DR_2_H, Hx_1, Hx_1_ghost
      Vec,                                    SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
      Vec,                                    SAVE :: Hb_3_145, Hb_3_236, Hx_3, Hx_3_ghost

      INTEGER :: rank, nb

      TYPE(solver_param), SAVE :: my_par_DR
      INTEGER, SAVE :: m_max_c
      INTEGER, SAVE :: vv_m_max_pad
      INTEGER, SAVE :: H_m_max_pad
      INTEGER, SAVE :: vv_nb_procs_F, vv_nb_procs_S, vv_rang_S, vv_rang_F
      INTEGER, SAVE :: H_nb_procs_F, H_nb_procs_S, H_rang_S, H_rang_F
      INTEGER, SAVE :: vv_bloc_size
      INTEGER, SAVE :: vv_bloc_size_gauss
      INTEGER, SAVE :: H_bloc_size
      INTEGER, SAVE :: H_bloc_size_gauss

      !-------for regularization
      Mat, DIMENSION(:), POINTER, SAVE :: vel_mat
      KSP, DIMENSION(:), POINTER, SAVE :: vel_ksp
      Mat, DIMENSION(:), POINTER, SAVE :: H_mat
      KSP, DIMENSION(:), POINTER, SAVE :: H_ksp

      !------for gauss to nodes
      Mat,                        SAVE :: vv_cc_DR_mat
      KSP,                        SAVE :: vv_cc_DR_ksp
      Mat,                        SAVE :: H_cc_DR_mat
      KSP,                        SAVE :: H_cc_DR_ksp

      REAL(KIND = 8), SAVE :: Hhh
      !--------global general variables
      LOGICAL, SAVE :: once = .TRUE.
      LOGICAL, SAVE :: once_lc = .TRUE.
      LOGICAL, SAVE :: once_lc_H = .TRUE.

      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: once_coef
      LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: once_coef_H

      REAL(KIND = 8), DIMENSION(:), POINTER, SAVE :: vv_one
      REAL(KIND = 8), DIMENSION(:), POINTER, SAVE :: vv_zero
      REAL(KIND = 8), DIMENSION(:), POINTER, SAVE :: H_one
      REAL(KIND = 8), DIMENSION(:), POINTER, SAVE :: H_zero

      TYPE(dyn_real_ksp), dimension(:), ALLOCATABLE, SAVE :: vel_ksp_many
      TYPE(dyn_real_mat), dimension(:), ALLOCATABLE, SAVE :: vel_mat_many
      TYPE(dyn_real_ksp), dimension(:), ALLOCATABLE, SAVE :: H_ksp_many
      TYPE(dyn_real_mat), dimension(:), ALLOCATABLE, SAVE :: H_mat_many
      TYPE(dyn_int_line), DIMENSION(3), SAVE :: vv_js_D_DR
      TYPE(dyn_int_line), DIMENSION(3), SAVE :: H_js_D_DR
      TYPE(dyn_int_line), DIMENSION(:), POINTER, SAVE :: vv_mode_global_js_D_DR
      TYPE(dyn_int_line), DIMENSION(:), POINTER, SAVE :: H_mode_global_js_D_DR
      TYPE(dyn_real_line), DIMENSION(:), ALLOCATABLE, SAVE :: vel_global_D_DR
      TYPE(dyn_real_line), DIMENSION(:), ALLOCATABLE, SAVE :: mag_global_D_DR



      !-------writing in file
      REAL(KIND = 8) :: scal, Em_mean
      INTEGER :: lblank
      CHARACTER(len = 1) :: tit1, tit2
      CHARACTER(len = 3) :: tit_DR
      CHARACTER(len = 4) ::  tit_S, tit_F
      CHARACTER(len = 250) :: out_name
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE, SAVE :: means


      !----------USEFUL VARIABLES-----------------------------------------------------------------------------
      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)) :: field_gauss_temp
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode))  :: field_node_temp
      REAL(KIND = 8), DIMENSION( vv_mesh%np, 6, SIZE(list_mode)) :: v_grad_v_nodes
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode))  :: res_scalar_nodes_vv

      REAL(KIND = 8) :: ray, res_integral, norm_loc
      INTEGER, DIMENSION(1) :: zero_mode
      INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
      INTEGER :: m_i, m_j, m_k
      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me) :: vv_r_dr_dth_dz
      ! REAL(KIND = 8), TARGET, ALLOCATABLE, DIMENSION( :, :, :) :: field_vv_node_i,field_vv_node_j,field_vv_node_k
      REAL(KIND = 8), TARGET, ALLOCATABLE, DIMENSION( :, :, :) :: field_H_node_i, field_H_node_j,field_H_node_k 
      REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)) :: field_vv_gauss_i,field_vv_gauss_j,field_vv_gauss_k
      REAL(KIND = 8), DIMENSION(3, vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)) :: Gu_gauss_i,Gu_gauss_j,Gu_gauss_k
! TEST VB 
      REAL(KIND = 8), DIMENSION(user%ending_it+1,user%ending_it+1) :: L_ij
      REAL(KIND = 8), DIMENSION(user%ending_it+1,user%ending_it+1,user%ending_it+1) :: Q_ijk
! TEST VB
      CHARACTER(len = 50) :: opt_I,opt_J,opt_K, opt_dir
      CHARACTER(len = 4) :: tit_I,tit_J,tit_K


      REAL(KIND = 8), DIMENSION(user%ending_it+1) :: u_symmetries
      REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)) :: vv_dummy
      REAL(KIND = 8), DIMENSION(6) :: type_sym



      !------MPI calls
      CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, code)
      CALL MPI_COMM_SIZE(communicator_ns(1), vv_nb_procs_S, code)
      CALL MPI_COMM_SIZE(communicator_ns(2), vv_nb_procs_F, code)
      CALL MPI_COMM_RANK(communicator_ns(1), vv_rang_S, code)
      CALL MPI_COMM_RANK(communicator_ns(2), vv_rang_F, code)
      CALL MPI_COMM_SIZE(communicator_mhd(1), H_nb_procs_S, code)
      CALL MPI_COMM_SIZE(communicator_mhd(2), H_nb_procs_F, code)
      CALL MPI_COMM_RANK(communicator_mhd(1), H_rang_S, code)
      CALL MPI_COMM_RANK(communicator_mhd(2), H_rang_F, code)



      IF(rank==0) THEN
      WRITE(*, *) '=============================================================='
      WRITE(*, *) 'BEGIN POD GALERKIN'
      WRITE(*, *) '=============================================================='
      END IF

      111 FORMAT(1500(e22.9, 2x))


      !------------------------------------------------------------------------------------------------
      !-------------INIT : THINGS DONE ONCE------------------------------------------------------------
      !------------------------------------------------------------------------------------------------

      IF (once) THEN
      once = .FALSE.
      ! IF(number_coeff_DR.NE.1) CALL error_petsc('Helmholtz: please start with first coeff')

      !Allocating means ofr saved quantities
      ALLOCATE(means(20 + user%nb_coeff_DR * 20))

      !-------------
      ALLOCATE(vv_one(vv_mesh%me), vv_zero(vv_mesh%np))

      vv_one(:) = 1.d0
      vv_zero(:) = 0.d0

      IF (user%if_recompute_DR) THEN
         ALLOCATE(once_coef(user%nb_coeff_DR))
         DO i = 1, user%nb_coeff_DR
            once_coef(i) = .TRUE.
         END DO
      END IF

      !-------------FFT VARIABLES----------------------------------------------------
      m_max_c = SIZE(list_mode)

      vv_bloc_size = vv_mesh%np / vv_nb_procs_F + 1
      vv_bloc_size_gauss = vv_mesh%gauss%l_G * vv_mesh%me / vv_nb_procs_F + 1

      vv_m_max_pad = 3 * m_max_c * vv_nb_procs_F / 2

      my_par_DR%verbose = .FALSE.
      my_par_DR%it_max = 100
      my_par_DR%rel_tol = 1.d-10
      my_par_DR%abs_tol = 1.d-10
      my_par_DR%solver = 'GMRES'
      my_par_DR%precond = 'MUMPS'

      !-----CREATE PETSC VECTORS AND GHOSTS-----------------------------------------
      CALL create_my_ghost(vv_mesh, vv_1_LA, vv_1_ifrom)
      n = vv_mesh%dom_np
      CALL VecCreateGhost(communicator_ns(1), n, &
            PETSC_DETERMINE, SIZE(vv_1_ifrom), vv_1_ifrom, vx_1, ierr)
      CALL VecGhostGetLocalForm(vx_1, vx_1_ghost, ierr)
      CALL VecDuplicate(vx_1, cb_DR_1, ierr)
      CALL VecDuplicate(vx_1, cb_DR_2, ierr)

      CALL create_my_ghost(vv_mesh, vv_3_LA, vv_3_ifrom)
      n = 3 * vv_mesh%dom_np
      CALL VecCreateGhost(communicator_ns(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
      CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
      CALL VecDuplicate(vx_3, vb_3_145, ierr)
      CALL VecDuplicate(vx_3, vb_3_236, ierr)


      !-------------ASSEMBLE VELOCITY AND MAG MATRICES FOR REG---------------------------------------
      ALLOCATE(vel_mat(2 * m_max_c), vel_ksp(2 * m_max_c))

      IF (user%if_recompute_DR) THEN
         ALLOCATE(vel_ksp_many(user%nb_coeff_DR), vel_mat_many(user%nb_coeff_DR))
         vel_ksp_many(1)%DKSP => vel_ksp
         vel_mat_many(1)%DMat => vel_mat

         DO i = 2, user%nb_coeff_DR
            ALLOCATE(vel_ksp_many(i)%DKSP(2 * m_max_c))
            ALLOCATE(vel_mat_many(i)%DMat(2 * m_max_c))
         END DO
      END IF

      !---------PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR !VELOCITY-----
      CALL vector_glob_js_D(vv_mesh, list_mode, vv_3_LA, user%vv_list_dirichlet_sides_DR, & !HF032019+CN LC
            vv_js_D_DR, vv_mode_global_js_D_DR)

      ALLOCATE(vel_global_D_DR(m_max_c))
      DO i = 1, m_max_c
         ALLOCATE(vel_global_D_DR(i)%DRL(SIZE(vv_mode_global_js_D_DR(i)%DIL)))
      END DO
      !------------------------------------------------------------------------------
      DO i = 1, m_max_c
         DO k = 1, 2
            nu_mat = 2 * (i - 1) + k
            CALL create_local_petsc_matrix(communicator_ns(1), vv_3_LA, vel_mat(nu_mat), clean = .FALSE.)
         END DO
      END DO

      !---------ASSEMBLING MASS MATRIX FOR gauss to nodes : Mass matrix identical for all i
      CALL create_local_petsc_matrix(communicator_ns(1), vv_1_LA, vv_cc_DR_mat, clean = .FALSE.)
      CALL qs_00_M (vv_mesh, 1.d0, vv_1_LA, vv_cc_DR_mat)
      IF (vvz_per%n_bord/=0) THEN
         CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vv_cc_DR_mat, vv_1_LA)
      END IF
      CALL init_solver(my_par_DR, vv_cc_DR_ksp, vv_cc_DR_mat, communicator_ns(1), &
            solver = my_par_DR%solver, precond = my_par_DR%precond)

      !---------Compute rr_gauss
      ALLOCATE(vv_rr_gauss(2, vv_mesh%gauss%l_G * vv_mesh%dom_me))
      index = 0
      DO m = 1, vv_mesh%dom_me
         j_loc = vv_mesh%jj(:, m)
         DO l = 1, vv_mesh%gauss%l_G
            index = index + 1
            vv_rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
            vv_rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))
         END DO
      END DO

      !-------Allocation for saved variables
      !===vv.mesh
      !gauss
      ALLOCATE(u_gauss(vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)))
      ALLOCATE(Gu_gauss(3, vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)))

      !nodes
      ALLOCATE(div_u_cube(vv_mesh%np, 6, SIZE(list_mode)))
      ALLOCATE(grad_u_square(vv_mesh%np, 6, SIZE(list_mode)))
      !HF20072020

      IF (user%if_B) THEN

         !-------------
         ALLOCATE(H_one(H_mesh%me), H_zero(H_mesh%np))

         H_one(:) = 1.d0
         H_zero(:) = 0.d0


         !-------------FFT VARIABLES----------------------------------------------------

         H_bloc_size = H_mesh%np / H_nb_procs_F + 1
         H_bloc_size_gauss = H_mesh%gauss%l_G * H_mesh%me / H_nb_procs_F + 1

         H_m_max_pad = 3 * m_max_c * H_nb_procs_F / 2



         !-----CREATE PETSC VECTORS AND GHOSTS-----------------------------------------
         CALL create_my_ghost(H_mesh, LA_H_1, H_1_ifrom)
         n = H_mesh%dom_np
         CALL VecCreateGhost(communicator_mhd(1), n, &
               PETSC_DETERMINE, SIZE(H_1_ifrom), H_1_ifrom, Hx_1, ierr)
         CALL VecGhostGetLocalForm(Hx_1, Hx_1_ghost, ierr)
         CALL VecDuplicate(Hx_1, cb_DR_1_H, ierr)
         CALL VecDuplicate(Hx_1, cb_DR_2_H, ierr)

         CALL create_my_ghost(H_mesh, LA_H_3, H_3_ifrom)
         n = 3 * H_mesh%dom_np
         CALL VecCreateGhost(communicator_mhd(1), n, &
               PETSC_DETERMINE, SIZE(H_3_ifrom), H_3_ifrom, Hx_3, ierr)
         CALL VecGhostGetLocalForm(Hx_3, Hx_3_ghost, ierr)
         CALL VecDuplicate(Hx_3, Hb_3_145, ierr)
         CALL VecDuplicate(Hx_3, Hb_3_236, ierr)


         !-------------ASSEMBLE VELOCITY AND MAG MATRICES FOR REG---------------------------------------

         ALLOCATE(H_mat(2 * m_max_c), H_ksp(2 * m_max_c))

         IF (user%if_recompute_DR) THEN
            ALLOCATE(H_ksp_many(user%nb_coeff_DR), H_mat_many(user%nb_coeff_DR))
            H_ksp_many(1)%DKSP => H_ksp
            H_mat_many(1)%DMat => H_mat

            DO i = 2, user%nb_coeff_DR
               ALLOCATE(H_ksp_many(i)%DKSP(2 * m_max_c))
               ALLOCATE(H_mat_many(i)%DMat(2 * m_max_c))
            END DO
         END IF

         !---------PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR !VELOCITY-----
         CALL vector_glob_js_D(H_mesh, list_mode, LA_H_3, user%H_list_dirichlet_sides_DR, & !HF032019+CN LC
               H_js_D_DR, H_mode_global_js_D_DR)

         ALLOCATE(mag_global_D_DR(m_max_c))
         DO i = 1, m_max_c
            ALLOCATE(mag_global_D_DR(i)%DRL(SIZE(H_mode_global_js_D_DR(i)%DIL)))
         END DO

         !------------------------------------------------------------------------------
         DO i = 1, m_max_c
            DO k = 1, 2
               nu_mat = 2 * (i - 1) + k
               CALL create_local_petsc_matrix(communicator_mhd(1), LA_H_3, H_mat(nu_mat), clean = .FALSE.)
            END DO
         END DO

         !---------ASSEMBLING MASS MATRIX FOR gauss to nodes : Mass matrix identical for all i
         CALL create_local_petsc_matrix(communicator_mhd(1), LA_H_1, H_cc_DR_mat, clean = .FALSE.)

         CALL qs_00_M (H_mesh, 1.d0, LA_H_1, H_cc_DR_mat)
         !IF (vvz_per%n_bord/=0) THEN
         !    CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vv_cc_DR_mat, vv_1_LA)
         !END IF

         CALL init_solver(my_par_DR, H_cc_DR_ksp, H_cc_DR_mat, communicator_mhd(1), &
               solver = my_par_DR%solver, precond = my_par_DR%precond)

         !---------Compute rr_gauss
         ALLOCATE(H_rr_gauss(2, H_mesh%gauss%l_G * H_mesh%dom_me))
         index = 0
         DO m = 1, H_mesh%dom_me
            j_loc = H_mesh%jj(:, m)
            DO l = 1, H_mesh%gauss%l_G
               index = index + 1
               H_rr_gauss(1, index) = SUM(H_mesh%rr(1, j_loc) * H_mesh%gauss%ww(:, l))
               H_rr_gauss(2, index) = SUM(H_mesh%rr(2, j_loc) * H_mesh%gauss%ww(:, l))
            END DO
         END DO

         !-------Allocation for saved variables

         !===vv.mesh
         !gauss

         !nodes
         ALLOCATE(Bn_vv(vv_mesh%np, 6, SIZE(list_mode)))
         ALLOCATE(u_cross_b(vv_mesh%np, 6, SIZE(list_mode)))



         !===H.mesh
         !gauss

         !nodes
         ALLOCATE(one_over_mu_H(H_mesh%np, 2, SIZE(list_mode)))
         ALLOCATE(rot_B_H(H_mesh%np, 6, SIZE(list_mode)))
         ALLOCATE(j_cross_b_H(H_mesh%np, 6, SIZE(list_mode)))
         ALLOCATE(grad_one_over_mu(H_mesh%np, 6, SIZE(list_mode)))
         ALLOCATE(j_H(H_mesh%np, 6, SIZE(list_mode)))
      END IF

      END IF ! ONCE
      

      !=============================================================================================
      !                       START COMPUTE L_ij
      !=============================================================================================
      !================================================================================
      !                             INIT OTHER VARIABLES
      !================================================================================
      
      
      vv_bloc_size_gauss = vv_mesh%gauss%l_G * vv_mesh%me / vv_nb_procs_F + 1
      m_max_c = SIZE(list_mode)
      vv_m_max_pad = 3 * m_max_c * vv_nb_procs_F / 2
      opt_dir = user%folder_for_binaries

      !---------Compute rr_gauss
      index = 0
      DO m = 1, vv_mesh%dom_me
         j_loc = vv_mesh%jj(:, m)
         DO l = 1, vv_mesh%gauss%l_G
            index = index + 1
            vv_rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
            vv_rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))
            vv_r_dr_dth_dz(index)= vv_rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
         END DO
      END DO

      !=============================================================================================
      !                       COMPUTING SYMMETRIES WHEN ASKED (TO INCREASE EFFICIENCY)
      !=============================================================================================

      IF (user%if_Rpisym) THEN
         ! SYMETRIE Rpi
         type_sym(1) = 1.d0
         type_sym(2) = -1.d0
         type_sym(3) = -1.d0
         type_sym(4) = 1.d0
         type_sym(5) = -1.d0
         type_sym(6) = 1.d0
         ! SYMETRIE Rpi
            DO m_i=user%starting_it, user%ending_it
            
         !=============================================================================================
         !                       FOR VELOCITY POD MODES
         !=============================================================================================
               IF (user%restart_from_bin) THEN
                  ! init OPT_I := {it:04d}
                  WRITE(tit_I, '(i4)') m_i
                  lblank = eval_blank(4, tit_I)
                  DO l = 1, lblank - 1
                     tit_I(l:l) = '0'
                  END DO
                  opt_I = '_I' // tit_I
   
                  ! get field on gauss point, in fourier representation
                  CALL FFT_INV_REAL_GAUSS_PDF(communicator_ns, vv_mesh, field_vv_gauss_k, vv_nb_procs_F, vv_bloc_size_gauss, user%field_name, opt_I)
                  CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_vv_gauss_k, field_node_temp, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)
                  IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                     WRITE(*,*) "In pod_galerkin.F90 -> finding symmetry of velocity " , m_i
                  END IF
               ELSE
                  WRITE(*,*) "It has to be done from gauss points, save gauss points in file, then retry !"
               END IF
   
   
               CALL symm_champ(communicator_ns(1), field_node_temp, vv_mesh, vv_dummy, 'u') ! needs to be performed on nodes
               
               !===Compute Rpi-symmetric field
               DO k = 1, 6
                  vv_dummy(:, k, :) = type_sym(k) * vv_dummy(:, k, :)
               END DO
               !===Compute Rpi-symmetric field
               
               CALL FFT_PAR_DOT_PROD_DCL(communicator_ns(2), field_node_temp, vv_dummy, res_scalar_nodes_vv, &
               vv_nb_procs_F, vv_bloc_size, vv_m_max_pad)
   
               CALL nodes_to_gauss(vv_mesh, res_scalar_nodes_vv, scalar_gauss) 
   
               ! volume integral == over mode zero !
               res_integral=0
               norm_loc=0
               IF (MINVAL(list_mode)==0) THEN
                  zero_mode = MINLOC(list_mode)
                  norm_loc = SUM( scalar_gauss(:,1,zero_mode(1)) * vv_r_dr_dth_dz )
               ELSE
                  norm_loc = 0
               END IF         
               CALL MPI_ALLREDUCE(norm_loc, res_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator_ns(1), code)

               u_symmetries(m_i + 1) = res_integral ! Close to +1 if symmetric, -1 if anti-symmetric due to normalization

            END DO
   
            ! EASIEST WAY TO MAKE SURE ALL PROCS KNOW THE SYMMETRY OF EACH MODE
            IF (vv_rang_S==0 .AND. vv_rang_F==0 ) THEN
               out_name = "./sym_u.fbin_stream"
               OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
               WRITE(10) u_symmetries(:)
               CLOSE(10)
            END IF
   
            OPEN(UNIT=10, FILE="./sym_u.fbin_stream", FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown') 
            DO m_i = user%starting_it, user%ending_it
               READ(10) u_symmetries(m_i + 1)
            END DO
            CLOSE(10)
   
         END IF
      !================================================================================
      !                               LOADING FIELDS
      !================================================================================
      IF (vv_rang_F==0 .AND. vv_rang_S==0) THEN
      WRITE(*,*) "In initialization.F90 -> run_pod : reading in", user%folder_for_binaries
      END IF

      IF (user%if_u) THEN
         DO m_i=user%starting_it, user%ending_it
            !======================================================
            !         START GET DATA for m_i
            !======================================================
            IF (user%restart_from_bin) THEN
               ! init OPT_I := {it:04d}
               WRITE(tit_I, '(i4)') m_i
               lblank = eval_blank(4, tit_I)
               DO l = 1, lblank - 1
                  tit_I(l:l) = '0'
               END DO
               opt_I = '_I' // tit_I
               
               IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                  WRITE(*,*) "In pod_galerkin.F90 -> compute_L_ij_Q_ijk : reading for it m_i" , m_i  
               END IF
               ! get field on gauss point, in fourier representation
               CALL FFT_INV_REAL_GAUSS_PDF(communicator_ns, vv_mesh, field_vv_gauss_i, vv_nb_procs_F, vv_bloc_size_gauss, user%field_name, opt_I)
               CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_vv_gauss_i, field_node_temp, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)

            ELSE
               WRITE(*,*) "In pod_galerkin, NOT IMPLEMENTED ERROR : Importing fields on node points!"
            END IF


            !======================================================
            !         END GET DATA for m_i
            !======================================================

            !======================================================
            !         START COMPUTE GRAD_U GAUSS for m_i
            !======================================================
            ! !-----computation on gauss points on vv_mesh
            DO i = 1, m_max_c
               mode = list_mode(i)
               index = 0
               DO m = 1, vv_mesh%dom_me
                  j_loc = vv_mesh%jj(:, m)
                  DO k = 1, 6
                     Vs(:, k) = field_node_temp(j_loc, k, i)
                  END DO

                  DO l = 1, vv_mesh%gauss%l_G
                     index = index + 1
                     dw_loc = vv_mesh%gauss%dw(:, :, l, m)
                     ray = vv_rr_gauss(1, index)
                     
                     !----------------- grad_u on gauss points-------------------------
                     !===d_i V_r
                     Gu_gauss_i(1, index, 1, i) = SUM(Vs(:, 1) * dw_loc(1, :))     ! d_r V_r^c
                     Gu_gauss_i(1, index, 2, i) = SUM(Vs(:, 2) * dw_loc(1, :))    ! d_r V_r^s
                     Gu_gauss_i(1, index, 3, i) = mode / ray * field_vv_gauss_i(index, 2, i)&!1/r d_th V_r^s
                           - 1.d0 / ray * field_vv_gauss_i(index, 3, i) ! -1/r V_th^c
                     Gu_gauss_i(1, index, 4, i) = -mode / ray * field_vv_gauss_i(index, 1, i)&!1/r d_th V_r^c
                           - 1.d0 / ray * field_vv_gauss_i(index, 4, i) ! -1/r V_th^s
                     Gu_gauss_i(1, index, 5, i) = SUM(Vs(:, 1) * dw_loc(2, :))     ! d_z V_r^c
                     Gu_gauss_i(1, index, 6, i) = SUM(Vs(:, 2) * dw_loc(2, :))     ! d_z V_r^s
                     !===d_i V_th
                     Gu_gauss_i(2, index, 1, i) = SUM(Vs(:, 3) * dw_loc(1, :))     ! d_r V_th^c
                     Gu_gauss_i(2, index, 2, i) = SUM(Vs(:, 4) * dw_loc(1, :))     ! d_r V_th^s
                     Gu_gauss_i(2, index, 3, i) = mode / ray * field_vv_gauss_i(index, 4, i)&!1/r d_th V_th^s
                           + 1.d0 / ray * field_vv_gauss_i(index, 1, i) ! +1/r V_r^c
                     Gu_gauss_i(2, index, 4, i) = -mode / ray * field_vv_gauss_i(index, 3, i)&!1/r d_th V_th^c
                           + 1.d0 / ray * field_vv_gauss_i(index, 2, i) ! +1/r V_r^s
                     Gu_gauss_i(2, index, 5, i) = SUM(Vs(:, 3) * dw_loc(2, :))     ! d_z V_th^c
                     Gu_gauss_i(2, index, 6, i) = SUM(Vs(:, 4) * dw_loc(2, :))     ! d_z V_th^s
                     !===d_i V_z
                     Gu_gauss_i(3, index, 1, i) = SUM(Vs(:, 5) * dw_loc(1, :))     ! d_r V_z^c
                     Gu_gauss_i(3, index, 2, i) = SUM(Vs(:, 6) * dw_loc(1, :))     ! d_r V_z^s
                     Gu_gauss_i(3, index, 3, i) = mode / ray * field_vv_gauss_i(index, 6, i) !1/r d_th V_z^s
                     Gu_gauss_i(3, index, 4, i) = -mode / ray * field_vv_gauss_i(index, 5, i) !1/r d_th V_z^c
                     Gu_gauss_i(3, index, 5, i) = SUM(Vs(:, 5) * dw_loc(2, :))     ! d_z V_z^c
                     Gu_gauss_i(3, index, 6, i) = SUM(Vs(:, 6) * dw_loc(2, :))     ! d_z V_z^s

                  END DO
               END DO
            END DO

            IF (vv_rang_F==0 .AND. vv_rang_S==0) THEN
            WRITE(*,*) "In pod_galerkin.F90 -> compute_L_ij_Q_ijk : first gradient passed"
            END IF
            !======================================================
            !         END COMPUTE GRAD_U GAUSS for m_i
            !======================================================
      
            !======================================================
            !                 AND VOILA for m_i
            !                    GET DATA m_j
            !======================================================  
            DO m_j=user%starting_it, user%ending_it

               !======================================================
               !         START GET DATA for m_j
               !======================================================
               IF (user%restart_from_bin) THEN
                  ! init OPT_I := {it:04d}
                  WRITE(tit_J, '(i4)') m_j
                  lblank = eval_blank(4, tit_J)
                  DO l = 1, lblank - 1
                     tit_J(l:l) = '0'
                  END DO
                  opt_J = '_I' // tit_J

                  IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                     WRITE(*,*) "In pod_galerkin.F90 -> compute_L_ij_Q_ijk : reading for it m_j" , m_j   
                  END IF 
                  ! get field on gauss point, in fourier representation
                  CALL FFT_INV_REAL_GAUSS_PDF(communicator_ns, vv_mesh, field_vv_gauss_j, vv_nb_procs_F, vv_bloc_size_gauss, user%field_name, opt_J)
                  CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_vv_gauss_j, field_node_temp, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)
               ELSE
                  WRITE(*,*) "It has to be done from gauss points, save gauss points in file, then retry !"
               END IF
               !======================================================
               !         END GET DATA for m_j
               !======================================================

               IF ((user%if_Rpisym == .false.) .OR. (user%if_Rpisym .AND. u_symmetries(m_i + 1) * u_symmetries(m_j + 1) > 0)) THEN
                  IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                     WRITE(*,*) "In pod_galerkin.F90 -> compute Lij : m_j " , m_j ,' compatible with ' , m_i  
                  END IF
               !======================================================
               !         START COMPUTE GRAD_U GAUSS for m_j
               !======================================================
               ! !-----computation on gauss points on vv_mesh

                  DO i = 1, m_max_c
                     mode = list_mode(i)
                     index = 0
                     DO m = 1, vv_mesh%dom_me
                        j_loc = vv_mesh%jj(:, m)
                        DO k = 1, 6
                           Vs(:, k) = field_node_temp(j_loc, k, i)
                        END DO
                        DO l = 1, vv_mesh%gauss%l_G
                           index = index + 1
                           dw_loc = vv_mesh%gauss%dw(:, :, l, m)
                           ray = vv_rr_gauss(1, index)
                           
                           !----------------- grad_u on gauss points-------------------------
                           !===d_i V_r
                           Gu_gauss_j(1, index, 1, i) = SUM(Vs(:, 1) * dw_loc(1, :))     ! d_r V_r^c
                           Gu_gauss_j(1, index, 2, i) = SUM(Vs(:, 2) * dw_loc(1, :))    ! d_r V_r^s
                           Gu_gauss_j(1, index, 3, i) = mode / ray * field_vv_gauss_j(index, 2, i)&!1/r d_th V_r^s
                                 - 1.d0 / ray * field_vv_gauss_j(index, 3, i) ! -1/r V_th^c
                           Gu_gauss_j(1, index, 4, i) = -mode / ray * field_vv_gauss_j(index, 1, i)&!1/r d_th V_r^c
                                 - 1.d0 / ray * field_vv_gauss_j(index, 4, i) ! -1/r V_th^s
                           Gu_gauss_j(1, index, 5, i) = SUM(Vs(:, 1) * dw_loc(2, :))     ! d_z V_r^c
                           Gu_gauss_j(1, index, 6, i) = SUM(Vs(:, 2) * dw_loc(2, :))     ! d_z V_r^s
                           !===d_i V_th
                           Gu_gauss_j(2, index, 1, i) = SUM(Vs(:, 3) * dw_loc(1, :))     ! d_r V_th^c
                           Gu_gauss_j(2, index, 2, i) = SUM(Vs(:, 4) * dw_loc(1, :))     ! d_r V_th^s
                           Gu_gauss_j(2, index, 3, i) = mode / ray * field_vv_gauss_j(index, 4, i)&!1/r d_th V_th^s
                                 + 1.d0 / ray * field_vv_gauss_j(index, 1, i) ! +1/r V_r^c
                           Gu_gauss_j(2, index, 4, i) = -mode / ray * field_vv_gauss_j(index, 3, i)&!1/r d_th V_th^c
                                 + 1.d0 / ray * field_vv_gauss_j(index, 2, i) ! +1/r V_r^s
                           Gu_gauss_j(2, index, 5, i) = SUM(Vs(:, 3) * dw_loc(2, :))     ! d_z V_th^c
                           Gu_gauss_j(2, index, 6, i) = SUM(Vs(:, 4) * dw_loc(2, :))     ! d_z V_th^s
                           !===d_i V_z
                           Gu_gauss_j(3, index, 1, i) = SUM(Vs(:, 5) * dw_loc(1, :))     ! d_r V_z^c
                           Gu_gauss_j(3, index, 2, i) = SUM(Vs(:, 6) * dw_loc(1, :))     ! d_r V_z^s
                           Gu_gauss_j(3, index, 3, i) = mode / ray * field_vv_gauss_j(index, 6, i) !1/r d_th V_z^s
                           Gu_gauss_j(3, index, 4, i) = -mode / ray * field_vv_gauss_j(index, 5, i) !1/r d_th V_z^c
                           Gu_gauss_j(3, index, 5, i) = SUM(Vs(:, 5) * dw_loc(2, :))     ! d_z V_z^c
                           Gu_gauss_j(3, index, 6, i) = SUM(Vs(:, 6) * dw_loc(2, :))     ! d_z V_z^s
                        END DO
                     END DO
                  END DO
                  !======================================================
                  !         END COMPUTE GRAD_U GAUSS for m_j
                  !======================================================
                  !======================================================
                  !    COMPUTE L_ij=GRAD_U(m_i):GRAD_U(m_j) terms
                  !======================================================
                  !tensor contraction
                  CALL FFT_TENSOR_CONTRACTION_gauss(communicator_ns(2), &
                  Gu_gauss_i(1, :, :, :), Gu_gauss_i(2, :, :, :), Gu_gauss_i(3, :, :, :), &
                  Gu_gauss_j(1, :, :, :), Gu_gauss_j(2, :, :, :), Gu_gauss_j(3, :, :, :), &
                     scalar_gauss, vv_nb_procs_F, vv_bloc_size_gauss, vv_m_max_pad)


                     

                  ! volume integral == over mode zero !
                  res_integral=0
                  norm_loc=0
                  IF (MINVAL(list_mode)==0) THEN
                     zero_mode = MINLOC(list_mode)
                     norm_loc = SUM( scalar_gauss(:,1,zero_mode(1)) * vv_r_dr_dth_dz )
                  ELSE
                     norm_loc = 0;
                  END IF         
                  CALL MPI_ALLREDUCE(norm_loc, res_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator_ns(1), code)

               ELSE
                  IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                     WRITE(*,*) "In pod_galerkin.F90 -> compute Lij : m_j " , m_j ,' not compatible with ' , m_i  
                  END IF
                  res_integral = 0
               END IF ! if user%if_Rpisym

               IF (vv_rang_S==0 .AND. vv_rang_F==0 ) THEN
                  L_ij(m_i + 1,m_j + 1) = res_integral
               END IF

               !======================================================
               !                 AND VOILA for m_j
               !======================================================
               DO m_k=user%starting_it, user%ending_it
                  IF ((user%if_Rpisym == .false.) .OR. (user%if_Rpisym .AND. u_symmetries(m_i + 1) * u_symmetries(m_j + 1) *u_symmetries(m_k + 1) > 0)) THEN
                     IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                        WRITE(*,*) "In pod_galerkin.F90 -> compute  Qijk : m_k " , m_k ,'compatible with ' , m_i , m_j  
                     END IF
                     !======================================================
                     !         START GET DATA for m_k
                     !======================================================
                     IF (user%restart_from_bin) THEN
                        ! init OPT_I := {it:04d}
                        WRITE(tit_K, '(i4)') m_k
                        lblank = eval_blank(4, tit_K)
                        DO l = 1, lblank - 1
                           tit_K(l:l) = '0'
                        END DO
                        opt_K = '_I' // tit_K

                        IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                           WRITE(*,*) "In pod_galerkin.F90 -> compute_L_ij_Q_ijk : reading for it m_k" , m_k   
                        END IF


                        ! get field on gauss point, in fourier representation
                        CALL FFT_INV_REAL_GAUSS_PDF(communicator_ns, vv_mesh, field_vv_gauss_k, vv_nb_procs_F, vv_bloc_size_gauss, user%field_name, opt_K)
                        CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_vv_gauss_k, field_node_temp, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)
                     ELSE
                        WRITE(*,*) "It has to be done from gauss points, save gauss points in file, then retry !"
                     END IF
                     !======================================================
                     !         END GET DATA for m_k
                     !======================================================

                     !======================================================
                     !         START COMPUTE GRAD_U GAUSS for m_k
                     !======================================================
                     ! !-----computation on gauss points on vv_mesh
                     DO i = 1, m_max_c
                        mode = list_mode(i)
                        index = 0
                        DO m = 1, vv_mesh%dom_me
                           j_loc = vv_mesh%jj(:, m)
                           DO k = 1, 6
                              Vs(:, k) = field_node_temp(j_loc, k, i)
                           END DO
                           DO l = 1, vv_mesh%gauss%l_G
                              index = index + 1
                              dw_loc = vv_mesh%gauss%dw(:, :, l, m)
                              ray = vv_rr_gauss(1, index)
                              
                              !----------------- grad_u on gauss points-------------------------
                              !===d_i V_r
                              Gu_gauss_k(1, index, 1, i) = SUM(Vs(:, 1) * dw_loc(1, :))     ! d_r V_r^c
                              Gu_gauss_k(1, index, 2, i) = SUM(Vs(:, 2) * dw_loc(1, :))    ! d_r V_r^s
                              Gu_gauss_k(1, index, 3, i) = mode / ray * field_vv_gauss_k(index, 2, i)&!1/r d_th V_r^s
                                    - 1.d0 / ray * field_vv_gauss_k(index, 3, i) ! -1/r V_th^c
                              Gu_gauss_k(1, index, 4, i) = -mode / ray * field_vv_gauss_k(index, 1, i)&!1/r d_th V_r^c
                                    - 1.d0 / ray * field_vv_gauss_k(index, 4, i) ! -1/r V_th^s
                              Gu_gauss_k(1, index, 5, i) = SUM(Vs(:, 1) * dw_loc(2, :))     ! d_z V_r^c
                              Gu_gauss_k(1, index, 6, i) = SUM(Vs(:, 2) * dw_loc(2, :))     ! d_z V_r^s
                              !===d_i V_th
                              Gu_gauss_k(2, index, 1, i) = SUM(Vs(:, 3) * dw_loc(1, :))     ! d_r V_th^c
                              Gu_gauss_k(2, index, 2, i) = SUM(Vs(:, 4) * dw_loc(1, :))     ! d_r V_th^s
                              Gu_gauss_k(2, index, 3, i) = mode / ray * field_vv_gauss_k(index, 4, i)&!1/r d_th V_th^s
                                    + 1.d0 / ray * field_vv_gauss_k(index, 1, i) ! +1/r V_r^c
                              Gu_gauss_k(2, index, 4, i) = -mode / ray * field_vv_gauss_k(index, 3, i)&!1/r d_th V_th^c
                                    + 1.d0 / ray * field_vv_gauss_k(index, 2, i) ! +1/r V_r^s
                              Gu_gauss_k(2, index, 5, i) = SUM(Vs(:, 3) * dw_loc(2, :))     ! d_z V_th^c
                              Gu_gauss_k(2, index, 6, i) = SUM(Vs(:, 4) * dw_loc(2, :))     ! d_z V_th^s
                              !===d_i V_z
                              Gu_gauss_k(3, index, 1, i) = SUM(Vs(:, 5) * dw_loc(1, :))     ! d_r V_z^c
                              Gu_gauss_k(3, index, 2, i) = SUM(Vs(:, 6) * dw_loc(1, :))     ! d_r V_z^s
                              Gu_gauss_k(3, index, 3, i) = mode / ray * field_vv_gauss_k(index, 6, i) !1/r d_th V_z^s
                              Gu_gauss_k(3, index, 4, i) = -mode / ray * field_vv_gauss_k(index, 5, i) !1/r d_th V_z^c
                              Gu_gauss_k(3, index, 5, i) = SUM(Vs(:, 5) * dw_loc(2, :))     ! d_z V_z^c
                              Gu_gauss_k(3, index, 6, i) = SUM(Vs(:, 6) * dw_loc(2, :))     ! d_z V_z^s

                           END DO
                        END DO
                     END DO
                     !======================================================
                     !         END COMPUTE GRAD_U GAUSS for m_k
                     !======================================================
                     
                     !======================================================
                     !              START COMPUTE Q_ijk
                     !====================================================== 
                  
                     !--------------Computation of (V.GRAD)V on gauss points as vector_gauss2----------------------
                     CALL FFT_B_DOT_GRAD_A_gauss(communicator_ns(2), field_vv_gauss_j, Gu_gauss_k(1, :, :, :), &
                     Gu_gauss_k(2, :, :, :), Gu_gauss_k(3, :, :, :), &
                     field_gauss_temp, vv_nb_procs_F, vv_bloc_size_gauss, vv_m_max_pad)

                     CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_gauss_temp, v_grad_v_nodes, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)
                     CALL gauss_to_nodes(vv_mesh, vv_1_LA, vv_one, vv_zero, field_vv_gauss_i, field_node_temp, vv_cc_DR_ksp, cb_DR_1, cb_DR_2, vx_1, vx_1_ghost, my_par_DR)

                     !---------Computation of scalar product of v.v as scalar_nodes_vv--------------------------------
                     CALL FFT_PAR_DOT_PROD_DCL(communicator_ns(2), field_node_temp, v_grad_v_nodes, res_scalar_nodes_vv, &
                     vv_nb_procs_F, vv_bloc_size, vv_m_max_pad)

                     ! node to gauss
                     CALL nodes_to_gauss(vv_mesh, res_scalar_nodes_vv, scalar_gauss)

                     ! volume integral == over mode zero !
                     res_integral=0
                     norm_loc=0
                     IF (MINVAL(list_mode)==0) THEN
                        zero_mode = MINLOC(list_mode)
                        norm_loc = SUM( scalar_gauss(:,1,zero_mode(1)) * vv_r_dr_dth_dz )
                     ELSE
                        norm_loc = 0;
                     END IF         
                     CALL MPI_ALLREDUCE(norm_loc, res_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator_ns(1), code)

                  ELSE
                     IF ( (vv_rang_S == 0 ) .AND. (vv_rang_F==0) ) THEN
                        WRITE(*,*) "In pod_galerkin.F90 -> compute Qijk : m_k " , m_k ,'not compatible with ' , m_i , m_j  
                     END IF
                     res_integral = 0

                  END IF ! user%if_Rpisym

                  IF (vv_rang_S==0 .AND. vv_rang_F==0 ) THEN
                     Q_ijk(m_i + 1,m_j + 1,m_k + 1) = res_integral
                  END IF
                  
                  !======================================================
                  !              END COMPUTE Q_ijk
                  !======================================================

                  !======================================================
                  !              UPDATE OUTPUTS
                  !======================================================          

                  IF (vv_rang_S==0 .AND. vv_rang_F==0 ) THEN
                     out_name = "./L_ij.fbin_stream"
                     OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
                     WRITE(10) L_ij(:,:)
                     CLOSE(10)
                     out_name = "./Q_ijk.fbin_stream"
                     OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
                     WRITE(10) Q_ijk(:,:,:)
                     CLOSE(10)
                  END IF

                  !======================================================
                  !              END UPDATE OUTPUTS
                  !====================================================== 

               END DO
            END DO 
         END DO
      ELSE 
         WRITE(*,*) "========================== NOT IMPLEMENTED ERROR =========================="
         WRITE(*,*) "In initialization.F90 -> run_pod :: Implementation restricted to 'u' field "
         WRITE(*,*) "==========================================================================="      
      END IF
      
      !=============================================================================================
      !                       END COMPUTE L_ij
      !=============================================================================================

      !================================================================================
      !                      WRITING L_ij and Q_ijk         
      !================================================================================
      IF (vv_rang_S==0 .AND. vv_rang_F==0 ) THEN
         out_name = "./L_ij.fbin_stream"
         OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
         WRITE(10) L_ij(:,:)
         CLOSE(10)
         out_name = "./Q_ijk.fbin_stream"
         OPEN(UNIT = 10, FILE = out_name, FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
         WRITE(10) Q_ijk(:,:,:)
         CLOSE(10)
      END IF


   END SUBROUTINE compute_L_ij_and_Q_ijk


END MODULE pod_galerkin