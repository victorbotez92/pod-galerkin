MODULE user_data_module
   USE def_type_mesh
   TYPE personalized_data
      !===I declare my own data here==================================================
      !===suites_to_bins=============================================================
      LOGICAL :: if_suites_to_bins_fourier, if_suites_to_bins_phys, if_suites_to_paraview, read_all_suites_from_file, bins_to_paraview
      LOGICAL :: if_suites_to_bins_fourier_per_mode
      LOGICAL :: if_u, if_B, if_p
      LOGICAL :: restart_from_bin,restart_galerkin
      INTEGER :: starting_it, ending_it, suites_to_bins_freq
      LOGICAL :: read_old_nsample_suite_fourier, read_old_nsample_suite_phys, suites_to_mean_suite, mean_suite_renorm
      INTEGER :: nsample_read, nsample_name
      CHARACTER(len = 250) :: file_for_suite_to_bins
      CHARACTER(len = 250) :: folder_for_binaries, field_name, field_name_B
      LOGICAL :: mag_galerkin,vel_galerkin
      INTEGER :: nb_coeff_DR
      LOGICAL :: if_Rpisym
      LOGICAL :: if_DR_postpro
      LOGICAL :: if_recompute_DR

      INTEGER :: vv_nb_dirichlet_sides_DR
      INTEGER :: H_nb_dirichlet_sides_DR
      TYPE(dyn_int_line), DIMENSION(3) :: vv_list_dirichlet_sides_DR
      TYPE(dyn_int_line), DIMENSION(3) :: H_list_dirichlet_sides_DR

      
      !.......Continue here ................................
   END TYPE personalized_data
END MODULE user_data_module

MODULE user_data
   USE user_data_module
   IMPLICIT NONE
   PUBLIC :: read_user_data
   TYPE(personalized_data), PUBLIC :: user
   PRIVATE

CONTAINS

   SUBROUTINE read_user_data(data_file)
      USE my_util
      USE chaine_caractere
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: data_file
      INTEGER :: unit_file = 22
      LOGICAL :: test
      INTEGER :: k
      
    
      OPEN(UNIT = unit_file, FILE = data_file, FORM = 'formatted', STATUS = 'unknown')

      !===pod_galerking===============================================================
      CALL find_string(unit_file, '===Should we do it for u (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%if_u
      ELSE
         user%if_u = .false.
      END IF
      CALL find_string(unit_file, '===Should we do it for p (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%if_p
      ELSE
         user%if_p = .false.
      END IF
      CALL find_string(unit_file, '===Should we do it for B (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%if_B
      ELSE
         user%if_B = .false.
      END IF
      CALL find_string(unit_file, '===Starting index for I', test)
      IF (test) THEN
         READ (unit_file, *) user%starting_it
      ELSE
         user%starting_it = -1
      ENDIF
      CALL find_string(unit_file, '===Ending index for I', test)
      IF (test) THEN
         READ (unit_file, *) user%ending_it
      ELSE
         user%ending_it = -2
      ENDIF
      
      

      CALL find_string(unit_file, '===Should we restart from "binaries" file (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%restart_from_bin
      ELSE
         user%restart_from_bin = .false.
      END IF

      CALL find_string(unit_file, '===Is it a restart ? (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%restart_galerkin
      ELSE
         user%restart_galerkin = .false.
      END IF

      CALL find_string(unit_file, "===Should we do Galerkin with magnetic field ?", test)
      IF (test) THEN
         READ (unit_file, *) user%mag_galerkin
      ELSE
         user%mag_galerkin = .false.
      END IF
      CALL find_string(unit_file, "===How are the magnetic POD called ?", test)
                  ! IF (vv_rang_F==0 .AND. vv_rang_S==0) THEN
      IF (test .AND. user%mag_galerkin==.true.) THEN
         READ (unit_file, *) user%field_name_B
      END IF



      CALL find_string(unit_file, "===Should we do Galerkin with velocity field ?", test)
      IF (test) THEN
         READ (unit_file, *) user%vel_galerkin
      ELSE
         user%vel_galerkin = .true.
      END IF

      
      CALL find_string(unit_file, "===Specific field name ", test)
      IF (test) THEN
         READ (unit_file, *) user%field_name
      ELSE
         IF (user%if_p) THEN
            user%field_name = 'p'
         END IF 
         IF (user%if_B) THEN
            user%field_name = 'B'
         END IF 
         IF (user%if_u) THEN
            user%field_name = 'u'
         END IF 
      END IF

      CALL find_string(unit_file, '===Folder where to put/get binaries', test)
      IF (test) THEN
         READ (unit_file, *) user%folder_for_binaries
      ELSE
         user%folder_for_binaries = '.'
      END IF
      CALL find_string(unit_file, '===Should we do it from the data_for_suites file (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%read_all_suites_from_file
      ELSE
         user%read_all_suites_from_file = .false.
      END IF
      IF (user%read_all_suites_from_file .OR. user%bins_to_paraview) THEN
         CALL find_string(unit_file, '===Name of file containing all info of fields', test)
         IF (test) THEN
            READ (unit_file, *) user%file_for_suite_to_bins
         ELSE
            WRITE(*, *) 'Error: No file provided for suites_to_bins_fourier from a file containing all info of fields'
         END IF
      END IF

      CALL find_string(unit_file, '===Should we Rpi-symmetrize (true/false)', test)
      IF (test) THEN
         READ (unit_file, *) user%if_Rpisym
      ELSE
         user%if_Rpisym = .false.
      END IF

      CALL find_string(unit_file, '===Should we do DR post processing computation (true/false)', test)
      IF (test) THEN
         !===HF may 2020
         READ (unit_file, *) user%if_DR_postpro
      ELSE
         user%if_DR_postpro = .false.
      END IF
      CALL find_string(unit_file, '===Should we save helmolhtz matrices to save time (true/false)', test)
      IF(test) THEN
            READ (unit_file, *) user%if_recompute_DR
      ELSE
            user%if_recompute_DR = .false.
      END IF

      CALL find_string(unit_file, '===How many boundary pieces for Dirichlet BCs for DR_MHD?', test)
        IF (test) THEN
            READ(unit_file, *) user%H_nb_dirichlet_sides_DR
        ELSE
            user%H_nb_dirichlet_sides_DR = 0
        END IF
        IF (user%H_nb_dirichlet_sides_DR >0) THEN
            DO k = 1, 3
                ALLOCATE(user%H_list_dirichlet_sides_DR(k)%DIL(user%H_nb_dirichlet_sides_DR))
                CALL read_until(unit_file, '===List of boundary pieces for Dirichlet BCs for DR_MHD')
                READ(unit_file, *) user%H_list_dirichlet_sides_DR(k)%DIL
            END DO
        ELSE
            DO k = 1, 3
                ALLOCATE(user%H_list_dirichlet_sides_DR(k)%DIL(0))
            END DO
        END IF


      
      !===End template=================================================================
      
      !!! TODO : DATA EXAMPLE suite to bin + :
      !Sets the boundaries for the filtering (should be the same as the ones for the simulation)
      ! ===How many boundary pieces for Dirichlet BCs for DR
      ! 2
      ! ===List of boundary pieces for Dirichlet BCs for DR
      ! 10 4
      ! ===How many boundary pieces for Dirichlet BCs for DR_MHD?
      ! 2
      ! ===List of boundary pieces for Dirichlet BCs for DR_MHD
      ! 12 4

      CLOSE(unit_file)
   END SUBROUTINE read_user_data

END MODULE user_data
