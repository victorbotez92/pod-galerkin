MODULE projections

  USE def_type_mesh

  IMPLICIT NONE
  CONTAINS
  !---------------------------------------------------------------------------
  SUBROUTINE projection_mag_field(mesh_in, mesh_out, vn, connectivity_structure, if_restriction, coupling_variable) !projection of mag field on another mesh subroutine added
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                               :: mesh_in, mesh_out
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
    INTEGER, DIMENSION(:), INTENT(IN)                         :: connectivity_structure
    LOGICAL, INTENT(IN)                                       :: if_restriction
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
    INTEGER                                                   :: j, m

    DO j = 1, SIZE(connectivity_structure)
       IF (connectivity_structure(j) == -1) CYCLE
       IF (if_restriction) THEN
          coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
       ELSE
          coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
       END IF
    END DO
    IF (if_restriction) THEN
       IF (mesh_in%gauss%n_w/=mesh_out%gauss%n_w) THEN
          DO m = 1, mesh_out%me
             coupling_variable(mesh_out%jj(4,m),:,:) = (coupling_variable(mesh_out%jj(2,m),:,:) + coupling_variable(mesh_out%jj(3,m),:,:))/2
             coupling_variable(mesh_out%jj(5,m),:,:) = (coupling_variable(mesh_out%jj(3,m),:,:) + coupling_variable(mesh_out%jj(1,m),:,:))/2
             coupling_variable(mesh_out%jj(6,m),:,:) = (coupling_variable(mesh_out%jj(1,m),:,:) + coupling_variable(mesh_out%jj(2,m),:,:))/2 
          END DO
       END IF
    END IF

  END SUBROUTINE projection_mag_field

END MODULE projections
