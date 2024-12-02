MODULE pdf
    USE petsc
#include "petsc/finclude/petsc.h"
    USE input_data
    IMPLICIT NONE
    !  PUBLIC:: FFT_PAR_REAL_GAUSS_PDF
CONTAINS

    SUBROUTINE FFT_PAR_REAL_GAUSS_PDF(communicator, mesh, V1_in, nb_procs, bloc_size, &
            type_field, opt_I, opt_dir, opt_nb_plane, opt_norm)

        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data
        IMPLICIT NONE
        INCLUDE 'fftw3.f'
        ! Format: V_1in(1:np,1:6,1:m_max_c)
        ! INPUT ARE COSINE AND SINE COEFFICIENTS
        ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2

        CHARACTER(*), INTENT(IN) :: type_field
        CHARACTER(*), INTENT(IN) :: opt_I, opt_dir
        TYPE(mesh_type), INTENT(IN) :: mesh
        LOGICAL, OPTIONAL :: opt_norm

        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V1_in
        !    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE    :: V_out
        !TEST LC
        !TEST LC
        INTEGER, INTENT(IN) :: nb_procs, bloc_size
        !    INTEGER, INTENT(IN)                            :: rank_S, rank_F
        INTEGER, OPTIONAL :: opt_nb_plane
        INTEGER :: np, nb_field, &
                m_max, m_max_c, MPID, m_max_pad, N_r_pad, np_tot
        INTEGER(KIND = 8) :: fftw_plan_multi_c2r
        COMPLEX(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: cu
        INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
        REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: dist_field, combined_field
        REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: V_out

        INTEGER :: fft_dim, howmany, istride, ostride, idist, odist, l, lblank
        INTEGER, DIMENSION(1) :: dim, inembed, onembed
        LOGICAL, SAVE :: once_fft = .TRUE.
        CHARACTER(len = 1) :: tit_field
        !LOGICAL,  DIMENSION(11),                 SAVE   ::  &
        !once=[.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.]
        !TEST LC
        INTEGER :: rang_S, rang_F
        !   REAL(KIND=8)                      :: r, z
        !   INTEGER                           :: l, lblank, ne, type_field
        !   CHARACTER(len=3)                  :: tit_S, tit_F, tit_DR
        !   CHARACTER(len=250)                :: out_name
        !   CHARACTER(len=200)                :: filename
        !TEST LC
        ! Recall complexes must be rescaled
        ! End FFTW parameters
        !v5.0!#include "petsc/finclude/petsc.h"
        PetscErrorCode :: ierr
        !TEST LC
        MPI_Comm, DIMENSION(:), POINTER :: communicator


        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)

        np = SIZE(V1_in, 1)
        nb_field = SIZE(V1_in, 2)    ! Number of fields
        m_max_c = SIZE(V1_in, 3)    ! Number of complex (cosines + sines) coefficients per point
        m_max = m_max_c * nb_procs ! Number of complex coefficients per point per processor
        np_tot = nb_procs * bloc_size

        !TEST PDF LC CN
        IF (once_fft) THEN
            once_fft = .FALSE.
            !TEST LC DEBUG to uncomment
            !     WRITE(*,*) 'np,nb_field,m_max_c,m_max,np_tot ', np,nb_field,m_max_c,m_max,np_tot
            !TEST LC DEBUG to uncomment
        END IF
        !TEST PDF LC CN

        IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
            CALL error_petsc('Bug in FFT_PAR_REAL_GAUSS_PDF: MOD(nb_field,2)/=0 .OR. m_max_c==0')
        END IF

        !===Bloc_size is the number of points that are handled by one processor
        !===once the Fourier modes are all collected on the processor
        !   IF (MODULO(np,nb_procs)==0) THEN
        !      bloc_size = np/nb_procs
        !   ELSE
        !      CALL error_petsc('Bug in FFT_PAR_REAL: np is not a multiple of nb_procs')
        !   END IF

        IF (PRESENT(opt_nb_plane)) THEN
            IF (opt_nb_plane> 2 * m_max - 1) THEN
                m_max_pad = (opt_nb_plane + 1) / 2
            ELSE
                m_max_pad = m_max
            END IF
        ELSE
            m_max_pad = m_max
        END IF
        N_r_pad = 2 * m_max_pad - 1

        ALLOCATE(cu(m_max_pad, nb_field / 2, bloc_size))
        ALLOCATE(dist_field(m_max_c, nb_field, np_tot))
        ALLOCATE(combined_field(m_max_c, nb_field, np_tot))
        ALLOCATE(V_out(N_r_pad, nb_field / 2, bloc_size))

        DO i = 1, m_max_c
            dist_field(i, :, 1:np) = TRANSPOSE(V1_in(:, :, i))
        END DO

        IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 1.d100

        longueur_tranche = bloc_size * m_max_c * nb_field

        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
                MPID, communicator(2), code)

        cu = 0.d0
        DO n = 1, bloc_size
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                jindex = n + shiftc
                DO nf = 1, nb_field / 2
                    !===Put real and imaginary parts in a complex
                    !===nf=1,2,3 => V1_in
                    !===INPUT ARE COSINE AND SINE COEFFICIENTS
                    !===THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                    cu(shiftl + 1:shiftl + m_max_c, nf, n) = 0.5d0 * CMPLX(combined_field(:, 2 * nf - 1, jindex), &
                            -combined_field(:, 2 * nf, jindex), KIND = 8)
                END DO
            END DO
        END DO
        cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)
        !===Padding is done by initialization of cu: cu = 0
        !===This is equivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0

        !===Set the parameters for dfftw
        fft_dim = 1
        istride = 1
        ostride = 1
        idist = N_r_pad
        inembed(1) = N_r_pad
        DIM(1) = N_r_pad
        odist = m_max_pad
        onembed(1) = m_max_pad
        howmany = bloc_size * nb_field / 2

        !TEST LC
        !    IF (ALLOCATED(V_out)) DEALLOCATE(V_out)
        !    ALLOCATE(V_out(N_r_pad,nb_field/2,bloc_size))
        !TEST LC
        CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
                onembed, ostride, odist, V_out, inembed, istride, idist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_c2r)

        IF (present(opt_norm)) THEN
            IF (opt_norm) THEN
                V_out(:, 1, :) = SQRT(V_out(:, 1, :) * V_out(:, 1, :) + V_out(:, 2, :) * V_out(:, 2, :) + V_out(:, 3, :) * V_out(:, 3, :))
                CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, V_out(:, 1, :), type_field, opt_I, TRIM(ADJUSTL(opt_dir)))
            ELSE
                IF (nb_field / 2 > 1) THEN
                    DO n = 1, nb_field / 2
                        WRITE(tit_field, '(i1)') n
                        CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, V_out(:, n, :),  trim(type_field) // tit_field, opt_I,  TRIM(ADJUSTL(opt_dir)))
                    END DO
                ELSE
                    CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, V_out(:, 1, :), trim(type_field), opt_I, TRIM(ADJUSTL(opt_dir)))
                END IF
            END IF
        ELSE
            IF (nb_field / 2 > 1) THEN
                DO n = 1, nb_field / 2
                    WRITE(tit_field, '(i1)') n
                    CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, V_out(:, n, :),  trim(type_field) // tit_field, opt_I,  TRIM(ADJUSTL(opt_dir)))
                END DO
            ELSE
                CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, V_out(:, 1, :), trim(type_field), opt_I, TRIM(ADJUSTL(opt_dir)))
            END IF
        END IF


        ! clean up

        CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

        DEALLOCATE(cu, dist_field, combined_field)

    END SUBROUTINE FFT_PAR_REAL_GAUSS_PDF

    SUBROUTINE FFT_INV_REAL_GAUSS_PDF(communicator, mesh, V_out, nb_procs, bloc_size, type_field, opt_I, opt_nb_plane)
        !This a de-aliased version of the code, FEB 4, 2011, JLG
        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data
        IMPLICIT NONE
        INCLUDE 'fftw3.f'

        ! Format: V_1in(1:np,1:6,1:m_max_c)
        ! INPUT ARE COSINE AND SINE COEFFICIENTS
        ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: V_out
        CHARACTER(*), INTENT(IN) :: type_field
        CHARACTER(*), INTENT(IN) :: opt_I
        TYPE(mesh_type), INTENT(IN) :: mesh

        INTEGER, OPTIONAL :: opt_nb_plane

        INTEGER, INTENT(IN) :: bloc_size, nb_procs
        INTEGER :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad, m_max_pad
        INTEGER(KIND = 8) :: fftw_plan_multi_r2c
        COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size) :: intermediate
        REAL(KIND = 8), DIMENSION(SIZE(V_out, 3), 2 * SIZE(V_out, 2), bloc_size * nb_procs) :: dist_field
        COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size, SIZE(V_out, 3) * nb_procs) :: combined_prod_cu
        COMPLEX(KIND = 8), DIMENSION(SIZE(V_out, 2) / 2, bloc_size, SIZE(V_out, 3) * nb_procs) :: dist_prod_cu
        REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: V_read
        COMPLEX(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: prod_cu
        INTEGER :: i_field
        INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code
        REAL(KIND = 8) :: t
        ! FFTW parameters
        INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
        INTEGER, DIMENSION(1) :: dim, inembed, onembed
        CHARACTER(len = 1) :: tit_field

        ! Recall complexes must be rescaled
        ! End FFTW parameters
        !#include "petsc/finclude/petsc.h"
        MPI_Comm, DIMENSION(:), POINTER :: communicator



        N_r_pad = 2 * m_max_pad - 1
        np = SIZE(V_out, 1)
        nb_field = SIZE(V_out, 2) ! Number of fields
        m_max_c = SIZE(V_out, 3) ! Number of complex (cosines + sines) coefficients per point
        m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
        IF (PRESENT(opt_nb_plane)) THEN
            IF (opt_nb_plane> 2 * m_max - 1) THEN
                m_max_pad = (opt_nb_plane + 1) / 2
            ELSE
                m_max_pad = m_max
            END IF
        ELSE
            m_max_pad = m_max
        END IF
        N_r_pad = 2 * m_max_pad - 1
        np_tot = nb_procs * bloc_size
        ALLOCATE(V_read(N_r_pad, nb_field / 2, bloc_size))
        ALLOCATE(prod_cu(m_max_pad, SIZE(V_out, 2) / 2, bloc_size))
        IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
            WRITE(*, *) ' BUG '
            STOP
        END IF
        V_read = 0.d0
        IF (nb_field / 2 > 1) THEN
            DO n = 1, nb_field / 2
                WRITE(tit_field, '(i1)') n
                CALL READ_REAL_GAUSS_PDF(communicator, mesh, V_read(:, n, :), trim(type_field) // tit_field, opt_I)
            END DO
        ELSE
            CALL READ_REAL_GAUSS_PDF(communicator, mesh, V_read(:, 1, :), trim(type_field), opt_I)
        END IF
        fft_dim = 1; istride = 1; ostride = 1;
        idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
        odist = m_max_pad; onembed(1) = m_max_pad

        howmany = bloc_size * nb_field
        howmany = howmany / 2
        CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, V_read, &
                inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)

        CALL dfftw_execute(fftw_plan_multi_r2c)

        ! clean up
        CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

        !!$       prod_cu = prod_cu/N_r !Scaling
        prod_cu = (1.d0 / N_r_pad) * prod_cu !Scaling
        combined_prod_cu(:, :, 1) = prod_cu(1, :, :)
        DO n = 2, m_max
            !combined_prod_cu(:,:,n)=prod_cu(n,:,:)
            combined_prod_cu(:, :, n) = 2 * CONJG(prod_cu(n, :, :))
        END DO
        t = MPI_WTIME()
        longueur_tranche = bloc_size * m_max_c * nb_field
        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (combined_prod_cu, longueur_tranche, MPID, dist_prod_cu, longueur_tranche, &
                MPID, communicator(2), code)
        DO i = 1, m_max_c
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                intermediate = dist_prod_cu(:, :, shiftl + i)
                DO n = 1, bloc_size
                    IF (n + shiftc > np) CYCLE
                    DO i_field = 1, nb_field / 2
                        V_out(n + shiftc, i_field * 2 - 1, i) = REAL (intermediate(i_field, n), KIND = 8)
                        V_out(n + shiftc, i_field * 2, i) = AIMAG(intermediate(i_field, n))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE FFT_INV_REAL_GAUSS_PDF


    SUBROUTINE nodes_to_gauss(vv_mesh, field_in, field_out)
        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: vv_mesh
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field_in
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, SIZE(field_in, 2), SIZE(field_in, 3)), INTENT(OUT) :: field_out
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, SIZE(field_in, 2)) :: Vs
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        INTEGER :: m_max_c, i, index, m, k, l, nb_field
        m_max_c = SIZE(field_in, 3)
        nb_field = SIZE(field_in, 2)
        IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
            CALL error_petsc('nodes_to_gauss PROBLEM')
        END IF
        !=== Computation of dield_in on gauss points
        DO i = 1, m_max_c
            index = 0
            DO m = 1, vv_mesh%dom_me
                j_loc = vv_mesh%jj(:, m)

                DO k = 1, nb_field
                    Vs(:, k) = field_in(j_loc, k, i)
                END DO

                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1
                    DO k = 1, nb_field
                        field_out(index, k, i) = SUM(Vs(:, k) * vv_mesh%gauss%ww(:, l))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE nodes_to_gauss


    SUBROUTINE WRITE_REAL_GAUSS_PDF(communicator, mesh, field, type_field, opt_I, opt_dir)
        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data
        REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: field
        CHARACTER(*), INTENT(IN) :: type_field
        CHARACTER(200) :: filename
        CHARACTER(*), INTENT(IN) :: opt_I, opt_dir
        TYPE(mesh_type), INTENT(IN) :: mesh
        MPI_Comm, DIMENSION(:), POINTER :: communicator
        INTEGER :: code, nb_angles, ntot_gauss, nb_procs_S, nb_procs_F, lblank, k, l, rang_S, rang_F, bloc_size_gauss, n

        CHARACTER(len = 4) :: tit_S, tit_F
        CHARACTER(len = 250) :: out_name

        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

        bloc_size_gauss = mesh%gauss%l_G * mesh%me / nb_procs_F + 1

        WRITE(tit_S, '(i4)') rang_S
        lblank = eval_blank(4, tit_S)
        DO l = 1, lblank - 1
            tit_S(l:l) = '0'
        END DO
        CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
        filename = inputs%file_name
        out_name = TRIM(ADJUSTL(opt_dir)) // '/phys_' // trim(type_field) // '_S' // tit_S // trim(opt_I) // '.' // filename

        IF (rang_F==nb_procs_F - 1) THEN
            ntot_gauss = mesh%gauss%l_G * mesh%me
        ELSE
            ntot_gauss = rang_F * bloc_size_gauss + bloc_size_gauss
        END IF

        ntot_gauss = ntot_gauss - rang_F * bloc_size_gauss
        nb_angles = SIZE(field(:, 1))

        DO k = 1, nb_angles
            DO n = 1, nb_procs_F
                IF (rang_F == n - 1) THEN
                    IF (rang_F == 0 .AND. k == 1) THEN
                        OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                                FORM = 'unformatted', ACCESS = 'stream', STATUS = 'replace')
                    ELSE
                        OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                                FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
                    END IF
                    IF (rang_F == nb_procs_F - 1) THEN
                        WRITE(10) field(k, 1:ntot_gauss)
                    ELSE
                        WRITE(10) field(k, :)
                    END IF
                    CLOSE(10)

                END IF
                CALL MPI_BARRIER(communicator(2), code)
            END DO
        END DO

    END SUBROUTINE WRITE_REAL_GAUSS_PDF

    SUBROUTINE READ_REAL_GAUSS_PDF(communicator, mesh, field, type_field, opt_I)
        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data
        REAL(KIND = 8), DIMENSION(:, :) :: field
        REAL(KIND = 8), DIMENSION(SIZE(field, 1), SIZE(field, 2)) :: field_trash
        CHARACTER(*), INTENT(IN) :: type_field
        CHARACTER(200) :: filename
        CHARACTER(*), INTENT(IN) :: opt_I
        TYPE(mesh_type), INTENT(IN) :: mesh
        MPI_Comm, DIMENSION(:), POINTER :: communicator
        INTEGER :: code, nb_angles, ntot_gauss, ntot_gauss_last, nb_procs_S, nb_procs_F, lblank, k, l, rang_S, rang_F, bloc_size_gauss, n, tot

        CHARACTER(len = 4) :: tit_S, tit_F
        CHARACTER(len = 250) :: out_name

        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

        bloc_size_gauss = mesh%gauss%l_G * mesh%me / nb_procs_F + 1

        WRITE(tit_S, '(i4)') rang_S
        lblank = eval_blank(4, tit_S)
        DO l = 1, lblank - 1
            tit_S(l:l) = '0'
        END DO

        filename = inputs%file_name
        out_name = trim(user%folder_for_binaries) // '/phys_' // trim(type_field) // '_S' // tit_S // trim(opt_I) // '.' // filename

        ! added n_tot_last<ntot_gauss for last procF : RB MC 11/06/24 
        ntot_gauss_last = mesh%gauss%l_G * mesh%me - (nb_procs_F-1)*bloc_size_gauss
        ntot_gauss = bloc_size_gauss
        nb_angles = SIZE(field, 1)
        OPEN(UNIT = 10, FILE = out_name, &
                FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
        tot=0
        DO k = 1, nb_angles
            DO n = 1, nb_procs_F
                ! IF (n == nb_procs_F) THEN
                !     tot= tot + ntot_gauss_last
                !     WRITE(*,*) "RANK S", rang_S, "LAST RANK F ", rang_F, 'ntot_gauss =',ntot_gauss_last
                ! ELSE
                !     tot= tot + ntot_gauss
                !     WRITE(*,*) "RANK S", rang_S, "RANK F ", rang_F, 'ntot_gauss =',ntot_gauss
                ! END IF
                IF (rang_F == n - 1) THEN
                    IF (n == nb_procs_F) THEN
                    ! IF (rang_F == nb_procs_F - 1) THEN deleted , RB MC 12/06/24
                        READ(10) field(k, 1:ntot_gauss_last)
                    ELSE
                        READ(10) field(k, 1:ntot_gauss)
                    END IF
                ELSE
                    IF (n == nb_procs_F) THEN
                    ! IF (rang_F == nb_procs_F - 1) THEN deleted , RB MC 12/06/24
                        READ(10) field_trash(k, 1:ntot_gauss_last)
                    ELSE
                        READ(10) field_trash(k, 1:ntot_gauss)
                    END IF
                END IF
            END DO
        END DO
        CLOSE(10)
        ! WRITE(*,*) 'tot =',tot
    END SUBROUTINE READ_REAL_GAUSS_PDF

    SUBROUTINE WRITE_PLAN_GAUSS(communicator, mesh, field, type_field, opt_I, opt_dir)
        USE my_util
        USE def_type_mesh
        USE chaine_caractere
        USE user_data

        TYPE(mesh_type), INTENT(IN) :: mesh
        REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: field
        CHARACTER(*) :: type_field
        CHARACTER(*) :: opt_I, opt_dir
        MPI_Comm, DIMENSION(:), POINTER :: communicator
        INTEGER :: code, nb_procs_S, rang_S, lblank, l, ntot

        CHARACTER(len = 4) :: tit_S, tit_F
        CHARACTER(len = 250) :: out_name
        CHARACTER(len = 200) :: filename

        !===MPI comms
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_RANK(communicator(1), rang_S, code)

        WRITE(tit_S, '(i4)') rang_S
        lblank = eval_blank(4, tit_S)
        DO l = 1, lblank - 1
            tit_S(l:l) = '0'
        END DO

        CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
        filename = inputs%file_name
        out_name = TRIM(ADJUSTL(opt_dir)) // '/' // trim(type_field) // '_S' // tit_S // trim(opt_I) // '.' // filename
        OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', FORM = 'unformatted', ACCESS = 'stream', STATUS = 'replace')

        WRITE(10) field(:)

    END SUBROUTINE WRITE_PLAN_GAUSS


    SUBROUTINE WRITE_ALL_GAUSS_MODES(communicator, list_mode, mesh, field, name_field, opt_I, opt_dir)
        USE def_type_mesh
        USE chaine_caractere
        USE user_data

        TYPE(mesh_type), INTENT(IN) :: mesh
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: field
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        CHARACTER(*) :: name_field
        CHARACTER(*) :: opt_I, opt_dir

        INTEGER :: ntot, nb_field, nb_modes, rang_S, rang_F, nb_procs_F, nb_procs_S, code, lblank, l, n, i, k
        CHARACTER(len = 1) :: tit_n
        CHARACTER(len = 4) :: tit_S, tit_F
        CHARACTER(len = 250) :: out_name
        CHARACTER(len = 200) :: filename
        CHARACTER(len = 50) :: type_field

        MPI_Comm, DIMENSION(:), POINTER :: communicator

        ntot = SIZE(field, 1)
        nb_field = SIZE(field, 2)
        nb_modes = SIZE(field, 3)

        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)

        WRITE(tit_S, '(i4)') rang_S
        lblank = eval_blank(4, tit_S)
        DO l = 1, lblank - 1
            tit_S(l:l) = '0'
        END DO

        filename = inputs%file_name
        DO n = 1, nb_field
            IF (nb_field > 2) THEN
                IF (MOD(n, 2) == 0) THEN
                    WRITE(tit_n, '(i1)') n / 2
                    type_field = trim(name_field) // tit_n // 's'
                ELSE
                    WRITE(tit_n, '(i1)') (n + 1) / 2
                    type_field = trim(name_field) // tit_n // 'c'
                END IF
            ELSE
                IF (n == 1) THEN
                    type_field = trim(name_field) // 'c'
                ELSE
                    type_field = trim(name_field) // 's'
                END IF
            END IF

            WRITE(tit_F, '(i4)') list_mode(i)
            lblank = eval_blank(4, tit_F)
            DO l = 1, lblank - 1
                tit_F(l:l) = '0'
            END DO

            CALL system('mkdir -p ' // TRIM(ADJUSTL(opt_dir)))
            out_name = TRIM(ADJUSTL(opt_dir)) // '/fourier_' // trim(type_field) // '_S' // tit_S // trim(opt_I) // '.' // filename
            DO i = 1, nb_procs_F * nb_modes
                DO  k = 1, SIZE(list_mode)
                    IF (list_mode(k) == i - 1) THEN
                        IF (list_mode(k) == 0) THEN
                            OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                                    FORM = 'unformatted', ACCESS = 'stream', STATUS = 'replace')
                        ELSE
                            OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                                    FORM = 'unformatted', ACCESS = 'stream', STATUS = 'unknown')
                        END IF

                        WRITE(10) field(:, n, k)
                        CLOSE(10)
                    END IF
                    CALL MPI_BARRIER(communicator(2), code)
                END DO
            END DO
        END DO
    END SUBROUTINE WRITE_ALL_GAUSS_MODES

    Subroutine R_mrgrnk (XDONT, IRNGT)
        ! __________________________________________________________
        !   MRGRNK = Merge-sort ranking of an array
        !   For performance reasons, the first 2 passes are taken
        !   out of the standard loop, and use dedicated coding.
        ! __________________________________________________________
        ! _________________________________________________________
        Real(KIND = 8), Dimension (:), Intent (In) :: XDONT
        Integer, Dimension (:), Intent (Out) :: IRNGT
        ! __________________________________________________________
        Real(KIND = 8) :: XVALA, XVALB
        !
        Integer, Dimension (SIZE(IRNGT)) :: JWRKT
        Integer :: LMTNA, LMTNC, IRNG1, IRNG2
        Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
        !
        NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
        Select Case (NVAL)
        Case (:0)
            Return
        Case (1)
            IRNGT (1) = 1
            Return
        Case Default
            Continue
        End Select
        !
        !  Fill-in the index array, creating ordered couples
        !
        Do IIND = 2, NVAL, 2
            If (XDONT(IIND - 1) <= XDONT(IIND)) Then
                IRNGT (IIND - 1) = IIND - 1
                IRNGT (IIND) = IIND
            Else
                IRNGT (IIND - 1) = IIND
                IRNGT (IIND) = IIND - 1
            End If
        End Do
        If (Modulo(NVAL, 2) /= 0) Then
            IRNGT (NVAL) = NVAL
        End If
        !
        !  We will now have ordered subsets A - B - A - B - ...
        !  and merge A and B couples into     C   -   C   - ...
        !
        LMTNA = 2
        LMTNC = 4
        !
        !  First iteration. The length of the ordered subsets goes from 2 to 4
        !
        Do
            If (NVAL <= 2) Exit
            !
            !   Loop on merges of A and B into C
            !
            Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD + 4) > NVAL) Then
                    If ((IWRKD + 2) >= NVAL) Exit
                    !
                    !   1 2 3
                    !
                    If (XDONT(IRNGT(IWRKD + 2)) <= XDONT(IRNGT(IWRKD + 3))) Exit
                    !
                    !   1 3 2
                    !
                    If (XDONT(IRNGT(IWRKD + 1)) <= XDONT(IRNGT(IWRKD + 3))) Then
                        IRNG2 = IRNGT (IWRKD + 2)
                        IRNGT (IWRKD + 2) = IRNGT (IWRKD + 3)
                        IRNGT (IWRKD + 3) = IRNG2
                        !
                        !   3 1 2
                        !
                    Else
                        IRNG1 = IRNGT (IWRKD + 1)
                        IRNGT (IWRKD + 1) = IRNGT (IWRKD + 3)
                        IRNGT (IWRKD + 3) = IRNGT (IWRKD + 2)
                        IRNGT (IWRKD + 2) = IRNG1
                    End If
                    Exit
                End If
                !
                !   1 2 3 4
                !
                If (XDONT(IRNGT(IWRKD + 2)) <= XDONT(IRNGT(IWRKD + 3))) Cycle
                !
                !   1 3 x x
                !
                If (XDONT(IRNGT(IWRKD + 1)) <= XDONT(IRNGT(IWRKD + 3))) Then
                    IRNG2 = IRNGT (IWRKD + 2)
                    IRNGT (IWRKD + 2) = IRNGT (IWRKD + 3)
                    If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD + 4))) Then
                        !   1 3 2 4
                        IRNGT (IWRKD + 3) = IRNG2
                    Else
                        !   1 3 4 2
                        IRNGT (IWRKD + 3) = IRNGT (IWRKD + 4)
                        IRNGT (IWRKD + 4) = IRNG2
                    End If
                    !
                    !   3 x x x
                    !
                Else
                    IRNG1 = IRNGT (IWRKD + 1)
                    IRNG2 = IRNGT (IWRKD + 2)
                    IRNGT (IWRKD + 1) = IRNGT (IWRKD + 3)
                    If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD + 4))) Then
                        IRNGT (IWRKD + 2) = IRNG1
                        If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD + 4))) Then
                            !   3 1 2 4
                            IRNGT (IWRKD + 3) = IRNG2
                        Else
                            !   3 1 4 2
                            IRNGT (IWRKD + 3) = IRNGT (IWRKD + 4)
                            IRNGT (IWRKD + 4) = IRNG2
                        End If
                    Else
                        !   3 4 1 2
                        IRNGT (IWRKD + 2) = IRNGT (IWRKD + 4)
                        IRNGT (IWRKD + 3) = IRNG1
                        IRNGT (IWRKD + 4) = IRNG2
                    End If
                End If
            End Do
            !
            !  The Cs become As and Bs
            !
            LMTNA = 4
            Exit
        End Do
        !
        !  Iteration loop. Each time, the length of the ordered subsets
        !  is doubled.
        !
        Do
            If (LMTNA >= NVAL) Exit
            IWRKF = 0
            LMTNC = 2 * LMTNC
            !
            !   Loop on merges of A and B into C
            !
            Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                    If (JINDA >= NVAL) Exit
                    IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
                !
                !   Shortcut for the case when the max of A is smaller
                !   than the min of B. This line may be activated when the
                !   initial set is already close to sorted.
                !
                !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
                !
                !  One steps in the C subset, that we build in the final rank array
                !
                !  Make a copy of the rank array for the merge iteration
                !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
                !
                Do
                    IWRK = IWRK + 1
                    !
                    !  We still have unprocessed values in both A and B
                    !
                    If (XVALA > XVALB) Then
                        IRNGT (IWRK) = IRNGT (IINDB)
                        IINDB = IINDB + 1
                        If (IINDB > IWRKF) Then
                            !  Only A still with unprocessed values
                            IRNGT (IWRK + 1:IWRKF) = JWRKT (IINDA:LMTNA)
                            Exit
                        End If
                        XVALB = XDONT (IRNGT(IINDB))
                    Else
                        IRNGT (IWRK) = JWRKT (IINDA)
                        IINDA = IINDA + 1
                        If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                        XVALA = XDONT (JWRKT(IINDA))
                    End If
                    !
                End Do
            End Do
            !
            !  The Cs become As and Bs
            !
            LMTNA = 2 * LMTNA
        End Do
        !
        Return
        !
    End Subroutine R_mrgrnk

END MODULE pdf
