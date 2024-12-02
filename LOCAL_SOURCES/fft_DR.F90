MODULE fft_DR
    !===Hugues Faller
    USE petsc
    USE my_util
    USE input_data
#include "petsc/finclude/petsc.h"
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
CONTAINS
    SUBROUTINE FFT_B_DOT_GRAD_A_gauss(communicator, V1_in, V2_in, V3_in, V4_in, &
            V1_out, nb_procs, bloc_size, m_max_pad)
        !FFT(B.GRAD(A)) ONLY ON GAUSS POINTS!
        !This a de-aliased version of the code, FEB 4, 2011, JLG
        !===V1_in is B
        !===V2_in is first line of gradA (V2_in is a vector) (d_j u_r, j={r,theta,z})
        !===V2=(d_r A_r, 1/r d_th A_r - A_th/r, d_z A_r )
        !===V3_in is second line of gradA (d_j u_theta, j={r,theta,z})
        !===V3=(d_r A_th, 1/r d_th A_th + A_r/r, d_z A_th )
        !===V4_in is third line of gradA (d_j u_z, j={r,theta,z})
        !===V4=(d_r A_z, 1/r d_th A_z, d_z A_z )
        !=== CAREFULL : you need "full" grad a : with crossed terms
        IMPLICIT NONE
        ! Format: V_1in(1:np,1:6,1:m_max_c)
        ! INPUT ARE COSINE AND SINE COEFFICIENTS
        ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V1_in, V2_in, V3_in, V4_in
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: V1_out
        INTEGER, INTENT(IN) :: nb_procs, bloc_size, m_max_pad
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, SIZE(V1_in, 2) * 4 / 2, bloc_size) :: cu
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, SIZE(V1_in, 2) * 4 / 2, bloc_size) :: ru
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, SIZE(V1_in, 2) / 2, bloc_size) :: prod_ru_1
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, SIZE(V1_in, 2) / 2, bloc_size) :: prod_cu_1
        COMPLEX(KIND = 8), DIMENSION(SIZE(V1_in, 2) / 2, bloc_size) :: intermediate_1
        REAL(KIND = 8), DIMENSION(SIZE(V1_in, 3), 4 * SIZE(V1_in, 2), bloc_size * nb_procs) :: dist_field
        REAL(KIND = 8), DIMENSION(SIZE(V1_in, 3), 4 * SIZE(V1_in, 2), bloc_size * nb_procs) :: combined_field
        COMPLEX(KIND = 8), DIMENSION(SIZE(V1_in, 2) / 2, bloc_size, SIZE(V1_in, 3) * nb_procs) :: combined_prod_cu_1
        COMPLEX(KIND = 8), DIMENSION(SIZE(V1_in, 2) / 2, bloc_size, SIZE(V1_in, 3) * nb_procs) :: dist_prod_cu_1
        INTEGER :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
        INTEGER(KIND = 8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
        INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank, l, i_field
        ! FFTW parameters
        INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
        INTEGER, DIMENSION(1) :: dim, inembed, onembed
        ! Recall complexes must be rescaled
        ! End FFTW parameters
        !v5.0!#include "petsc/finclude/petsc.h"
        MPI_Comm :: communicator
        CALL MPI_COMM_RANK(communicator, rank, code)

        np = SIZE(V1_in, 1)
        nb_field = SIZE(V1_in, 2) ! Number of fields
        m_max_c = SIZE(V1_in, 3) ! Number of complex (cosines + sines) coefficients per point
        m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
        np_tot = nb_procs * bloc_size
        N_r_pad = 2 * m_max_pad - 1

        IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
            WRITE(*, *) ' BUG in FFT_B_DOT_GRAD_A_gauss '
            STOP
        END IF

        ! Packing all 3 complex components of input fields
        ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
        ! so that after distributing the data to the processes, each one will obtain a part
        ! on nodal points
        ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT

        DO i = 1, m_max_c
            dist_field(i, 1:nb_field, 1:np) = TRANSPOSE(V1_in(:, :, i))
            dist_field(i, nb_field + 1:2 * nb_field, 1:np) = TRANSPOSE(V2_in(:, :, i))
            dist_field(i, 2 * nb_field + 1:3 * nb_field, 1:np) = TRANSPOSE(V3_in(:, :, i))
            dist_field(i, 3 * nb_field + 1:4 * nb_field, 1:np) = TRANSPOSE(V4_in(:, :, i))
            !      dist_field(i,           1:nb_field,  1:np) = TRANSPOSE(V1_in(1:np,1:nb_field,i))
            !      dist_field(i,  nb_field+1:2*nb_field,1:np) = TRANSPOSE(V2_in(1:np,1:nb_field,i))
            !      dist_field(i,2*nb_field+1:3*nb_field,1:np) = TRANSPOSE(V3_in(1:np,1:nb_field,i))
            !      dist_field(i,3*nb_field+1:4*nb_field,1:np) = TRANSPOSE(V4_in(1:np,1:nb_field,i))
        END DO

        IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 0.d0

        longueur_tranche = bloc_size * m_max_c * nb_field * 4

        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
                MPID, communicator, code)

        !JLG, FEB 4, 2011
        cu = 0.d0
        !JLG, FEB 4, 2011
        DO n = 1, bloc_size
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                jindex = n + shiftc
                DO nf = 1, nb_field * 4 / 2
                    ! Put real and imaginary parts in a complex
                    ! nf=1,2,3 => V1_in
                    ! nf=4,5,6 => V2_in
                    ! nf=7,8,9 => V3_in
                    ! nf=10,11,12 => V4_in
                    ! INPUT ARE COSINE AND SINE COEFFICIENTS
                    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                    cu(shiftl + 1:shiftl + m_max_c, nf, n) = 0.5d0 * CMPLX(combined_field(:, 2 * nf - 1, jindex), &
                            -combined_field(:, 2 * nf, jindex), KIND = 8)
                END DO
            END DO
        END DO
        cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)


        ! Set the parameters for dfftw
        fft_dim = 1; istride = 1; ostride = 1;
        idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
        odist = m_max_pad; onembed(1) = m_max_pad

        howmany = bloc_size * nb_field * 4 / 2 !4 vectors in

        CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
                onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_c2r)
        CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

        !===ru(:,1:3,:)   is B in physical space
        !===ru(:,4:6,:)   is first line of gradA in physical space
        ! ie d_i(A_r)
        !===ru(:,7:9,:) is second line of gradA in physical space
        ! ie d_i(A_theta)
        !===ru(:,10:12,:) is third line of gradA in physical space
        ! ie d_i(A_z)

        !===Compute (B.GRAD)A
        IF (nb_field==6) THEN
            prod_ru_1(:, 1, :) = ru(:, 1, :) * ru(:, 4, :) & != B_r  d_r(A_r)
                    + ru(:, 2, :) * ru(:, 5, :) & !+ B_th (d_th(A_r)/r - A_th/r)
                    + ru(:, 3, :) * ru(:, 6, :)   !+ B_z  d_z(A_r)
            prod_ru_1(:, 2, :) = ru(:, 1, :) * ru(:, 7, :) & != B_r  d_r(A_th)
                    + ru(:, 2, :) * ru(:, 8, :) & !+ B_th (d_th(A_th)/r + A_r/r)
                    + ru(:, 3, :) * ru(:, 9, :)   !+ B_z  d_z(A_th)
            prod_ru_1(:, 3, :) = ru(:, 1, :) * ru(:, 10, :) &!= B_r  d_r(A_z)
                    + ru(:, 2, :) * ru(:, 11, :) &!+ B_th d_th(A_z)/r
                    + ru(:, 3, :) * ru(:, 12, :)  !+ B_z  d_z(A_z)
        END IF
        !===End (B.GRAD)A

        howmany = bloc_size * nb_field / 2 !1 vector out
        CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru_1, &
                inembed, istride, idist, prod_cu_1, onembed, ostride, odist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_r2c)
        CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

        prod_cu_1 = prod_cu_1 * (1.d0 / N_r_pad) !Scaling

        !Now we need to redistribute the Fourier coefficients on each processor
        combined_prod_cu_1(:, :, 1) = prod_cu_1(1, :, :)

        DO n = 2, m_max
            combined_prod_cu_1(:, :, n) = 2 * CONJG(prod_cu_1(n, :, :))
        END DO

        longueur_tranche = bloc_size * m_max_c * nb_field
        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (combined_prod_cu_1, longueur_tranche, MPID, dist_prod_cu_1, longueur_tranche, &
                MPID, communicator, code)

        DO i = 1, m_max_c
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                intermediate_1 = dist_prod_cu_1(:, :, shiftl + i)
                DO n = 1, bloc_size
                    IF (n + shiftc > np) CYCLE
                    DO i_field = 1, nb_field / 2
                        V1_out(n + shiftc, i_field * 2 - 1, i) = REAL (intermediate_1(i_field, n), KIND = 8)
                        V1_out(n + shiftc, i_field * 2, i) = AIMAG(intermediate_1(i_field, n))
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE FFT_B_DOT_GRAD_A_gauss

    SUBROUTINE FFT_PAR_POWER(communicator, c_in, p, if_abs, c_out, nb_procs, bloc_size, m_max_pad, temps)
        !FFT ((FFT(-1) c_in)**p)  = c_out
        !when if_abs is true, compute p-power of absolute value of FFT(-1) c_in [for
        !non integer p]
        !This a de-aliased version of the code, FEB 4, 2011, JLG
        IMPLICIT NONE
        ! Format: c_in(1:np,1:2,1:m_max_c)
        ! INPUT ARE COSINE AND SINE COEFFICIENTS
        ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: c_in
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: c_out
        REAL(KIND = 8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: temps
        INTEGER, INTENT(IN) :: nb_procs, bloc_size, m_max_pad
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, 1, bloc_size) :: cu
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, 1, bloc_size) :: ru
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, bloc_size) :: prod_cu
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, bloc_size) :: prod_ru
        REAL(KIND = 8), DIMENSION(SIZE(c_in, 3), SIZE(c_in, 2), bloc_size * nb_procs) :: dist_field, combined_field !LC comment: 4 becomes SIZE(c_in,2)
        COMPLEX(KIND = 8), DIMENSION(bloc_size, SIZE(c_in, 3) * nb_procs) :: dist_prod_cu, combined_prod_cu
        COMPLEX(KIND = 8), DIMENSION(bloc_size) :: intermediate
        INTEGER :: np, np_tot, m_max, m_max_c, MPID, N_r_pad, nb_field
        INTEGER(KIND = 8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
        INTEGER :: nb, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank
        REAL(KIND = 8) :: t
        ! FFTW parameters
        INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
        INTEGER, DIMENSION(1) :: dim, inembed, onembed
        ! Recall complexes must be rescaled
        !HFCN
        REAL(KIND = 8) :: p
        LOGICAL :: if_abs
        ! End FFTW parameters
        !v5.0!#include "petsc/finclude/petsc.h"
        MPI_Comm :: communicator

        IF (PRESENT(temps)) temps = 0.d0

        np = SIZE(c_in, 1)
        nb_field = SIZE(c_in, 2) !LC comment: definiton nb_field = nb_input components = 2*nb_scalar_in + 6*nb_vector_in
        m_max_c = SIZE(c_in, 3) ! Number of complex (cosines + sines) coefficients per point
        m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
        np_tot = nb_procs * bloc_size
        N_r_pad = 2 * m_max_pad - 1

        IF (m_max_c==0) THEN
            WRITE(*, *) ' BUG in FFT_PAR_POWER '
            STOP
        END IF

        ! Packing all 3 complex components of both v1 and v2 input fields
        ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
        ! so that after distributing the data to the processes, each one will obtain a part
        ! on nodal points
        ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT

        t = MPI_WTIME()
        DO i = 1, m_max_c
            dist_field(i, 1:nb_field, 1:np) = TRANSPOSE(c_in(:, :, i))  !LC comment: 1:nb_field
        END DO
        IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t

        IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 1.d100

        longueur_tranche = bloc_size * m_max_c * nb_field !LC comment: use nb_field=2 (1 scalar as input so 2 Fourier components)

        t = MPI_WTIME()
        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
                MPID, communicator, code)
        IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() - t

        t = MPI_WTIME()
        !JLG, FEB 4, 2011
        cu = 0.d0
        !JLG, FEB 4, 2011
        DO n = 1, bloc_size
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                jindex = n + shiftc
                ! Put real and imaginary parts in a complex
                ! INPUT ARE COSINE AND SINE COEFFICIENTS
                ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                cu(shiftl + 1:shiftl + m_max_c, 1, n) = 0.5d0 * CMPLX(combined_field(:, 1, jindex), &
                        -combined_field(:, 2, jindex), KIND = 8)
            END DO
        END DO
        cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)
        !JLG, FEB 4, 2011
        !Padding is done by initialization of cu: cu = 0
        !This is eequivalent to cu(m_max+1:m_max_pad,:,:) = 0.d0
        !JLG, FEB 4, 2011

        IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t

        ! Set the parameters for dfftw
        fft_dim = 1; istride = 1; ostride = 1;
        !JLG, FEB 4, 2011
        idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
        odist = m_max_pad; onembed(1) = m_max_pad
        !JLG, FEB 4, 2011

        howmany = bloc_size * nb_field / 2 !LC comment: input = one scalar so bloc_size*nb_field/2=bloc_size elements (vector would be times 3).
        !nb_field is divided by two because we talking in physical space (1 scalar = 1 component)

        t = MPI_WTIME()
        CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
                onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_c2r)

        ! clean up
        CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

        !POWER p
        IF(if_abs) THEN
            prod_ru(:, :) = ABS(ru(:, 1, :))**p
        ELSE
            !CALL MPI_COMM_RANK(PETSC_COMM_WORLD,rank,code)
            !if(rank==0)  write(*,*) 'p, intp, intp+0.1', p, int(p), int(p+0.1d0)
            prod_ru(:, :) = ru(:, 1, :)**(int(p + 0.1d0))
        END IF
        !POWER p

        howmany = howmany !LC comment:input = output = 1 scalar

        CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
                inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_r2c)

        ! clean up
        CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

        !JLG, FEB 4, 2011
        prod_cu = prod_cu * (1.d0 / N_r_pad) !Scaling
        !JLG, FEB 4, 2011
        IF (PRESENT(temps)) temps(2) = temps(2) + MPI_WTIME() - t

        !Now we need to redistribute the Fourier coefficients on each processor
        t = MPI_WTIME()
        combined_prod_cu(:, 1) = prod_cu(1, :)
        DO n = 2, m_max
            combined_prod_cu(:, n) = 2 * CONJG(prod_cu(n, :))
        END DO

        IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t

        t = MPI_WTIME()
        longueur_tranche = bloc_size * m_max_c * 2
        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (combined_prod_cu, longueur_tranche, MPID, dist_prod_cu, longueur_tranche, &
                MPID, communicator, code)
        IF (PRESENT(temps)) temps(1) = temps(1) + MPI_WTIME() - t

        t = MPI_WTIME()
        DO i = 1, m_max_c
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                intermediate = dist_prod_cu(:, shiftl + i)
                DO n = 1, bloc_size
                    IF (n + shiftc > np) CYCLE
                    c_out(n + shiftc, 1, i) = REAL (intermediate(n), KIND = 8)
                    c_out(n + shiftc, 2, i) = AIMAG(intermediate(n))
                END DO
            END DO
        END DO
        IF (PRESENT(temps)) temps(3) = temps(3) + MPI_WTIME() - t
    END SUBROUTINE FFT_PAR_POWER

    SUBROUTINE FFT_TENSOR_CONTRACTION_gauss(communicator, V1_in, V2_in, V3_in, &
            V4_in, V5_in, V6_in, S1_out, nb_procs, bloc_size, m_max_pad, &
            opt_rr, opt_rank_S, opt_coef)
        USE pdf
        !S1_out=FFT(FFT-1(A):FFT-1(B)) A B being tensors
        !===V1_in is first  ligne of tensor A(1,1:np,1:6,list_modes)
        !===   if gradient : d_i(A_r) i={r, theta, z}
        !===V2_in is second ligne of tensor A(2,1:np,1:6,list_modes)
        !===   if gradient : d_i(A_th) i={r,theta, z}
        !===V3_in is third  ligne of tensor A(3,1:np,1:6,list_modes)
        !===   if gradient : d_i(A_z) i={r, theta, z}
        !===V4_in is first  ligne of tensor B(1,1:np,1:6,list_modes)
        !===   if gradient : d_i(B_r) i={r, theta, z}
        !===V5_in is second ligne of tensor B(2,1:np,1:6,list_modes)
        !===   if gradient : d_i(B_th) i={r,theta, z}
        !===V6_in is third  ligne of tensor B(3,1:np,1:6,list_modes)
        !===   if gradient : d_i(B_z) i={r, theta, z}
        !===S1_out is a scaler : (1:np,1:2,list_mode)

        IMPLICIT NONE
        ! Format: V_1in(1:np,1:6,1:m_max_c)
        ! INPUT ARE COSINE AND SINE COEFFICIENTS
        ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V1_in, V2_in, V3_in, V4_in, V5_in, V6_in
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: S1_out
        INTEGER, INTENT(IN) :: nb_procs, bloc_size, m_max_pad
        REAL(Kind = 8), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: opt_rr
        INTEGER, OPTIONAL, INTENT(IN) :: opt_rank_S, opt_coef
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, SIZE(V1_in, 2) * 6 / 2, bloc_size) :: cu
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, SIZE(V1_in, 2) * 6 / 2, bloc_size) :: ru
        REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, 1, bloc_size) :: prod_ru
        COMPLEX(KIND = 8), DIMENSION(m_max_pad, bloc_size) :: prod_cu
        COMPLEX(KIND = 8), DIMENSION(bloc_size) :: intermediate
        REAL(KIND = 8), DIMENSION(SIZE(V1_in, 3), 6 * SIZE(V1_in, 2), bloc_size * nb_procs) :: dist_field
        REAL(KIND = 8), DIMENSION(SIZE(V1_in, 3), 6 * SIZE(V1_in, 2), bloc_size * nb_procs) :: combined_field
        COMPLEX(KIND = 8), DIMENSION(bloc_size, SIZE(V1_in, 3) * nb_procs) :: combined_prod_cu
        COMPLEX(KIND = 8), DIMENSION(bloc_size, SIZE(V1_in, 3) * nb_procs) :: dist_prod_cu
        INTEGER :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
        INTEGER(KIND = 8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
        INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank_F, l, i_field
        ! FFTW parameters
        INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
        INTEGER, DIMENSION(1) :: dim, inembed, onembed
        ! Recall complexes must be rescaled
        ! End FFTW parameters
        !v5.0!#include "petsc/finclude/petsc.h"
        MPI_Comm :: communicator
        CALL MPI_COMM_RANK(communicator, rank_F, code)

        np = SIZE(V1_in, 1)
        nb_field = SIZE(V1_in, 2) ! Number of fields
        m_max_c = SIZE(V1_in, 3) ! Number of complex (cosines + sines) coefficients per point
        m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
        np_tot = nb_procs * bloc_size
        N_r_pad = 2 * m_max_pad - 1

        IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
            WRITE(*, *) ' BUG in FFT_TENSOR_CONTRACTION'
            STOP
        END IF

        ! Packing all complex components of input fields
        ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
        ! so that after distributing the data to the processes, each one will obtain a part
        ! on nodal points
        ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
        DO i = 1, m_max_c
            dist_field(i, 1:nb_field, 1:np) = TRANSPOSE(V1_in(:, :, i))
            dist_field(i, nb_field + 1:2 * nb_field, 1:np) = TRANSPOSE(V2_in(:, :, i))
            dist_field(i, 2 * nb_field + 1:3 * nb_field, 1:np) = TRANSPOSE(V3_in(:, :, i))
            dist_field(i, 3 * nb_field + 1:4 * nb_field, 1:np) = TRANSPOSE(V4_in(:, :, i))
            dist_field(i, 4 * nb_field + 1:5 * nb_field, 1:np) = TRANSPOSE(V5_in(:, :, i))
            dist_field(i, 5 * nb_field + 1:6 * nb_field, 1:np) = TRANSPOSE(V6_in(:, :, i))
        END DO

        !IF (np/=np_tot) dist_field(:,:,np+1:np_tot) = 1.d100
        !LC 2016/02/16 (we do a max on vel later)
        IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 0.d0
        !LC 2016/02/16

        longueur_tranche = bloc_size * m_max_c * nb_field * 6

        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
                MPID, communicator, code)

        !JLG, FEB 4, 2011
        cu = 0.d0
        !JLG, FEB 4, 2011
        DO n = 1, bloc_size
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                jindex = n + shiftc
                DO nf = 1, nb_field * 6 / 2
                    ! Put real and imaginary parts in a complex
                    ! nf= 1, 2, 3 => V1_in
                    ! nf= 4, 5, 6 => V2_in
                    ! nf= 7, 8, 9 => V3_in
                    ! nf=10,11,12 => V4_in
                    ! nf=13,14,15 => V5_in
                    ! nf=16,17,18 => V6_in
                    ! INPUT ARE COSINE AND SINE COEFFICIENTS
                    ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
                    cu(shiftl + 1:shiftl + m_max_c, nf, n) = 0.5d0 * CMPLX(combined_field(:, 2 * nf - 1, jindex), &
                            -combined_field(:, 2 * nf, jindex), KIND = 8)
                END DO
            END DO
        END DO
        cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)

        ! Set the parameters for dfftw
        fft_dim = 1; istride = 1; ostride = 1;
        idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
        odist = m_max_pad; onembed(1) = m_max_pad

        howmany = bloc_size * nb_field * 6 / 2

        CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
                onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_c2r)

        ! clean up
        CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

        !===ru(:, 1:3 ,:)   is A(1,1), A(1,2) and A(1,3) in physical space.
        !===ru(:, 4:6 ,:)   is A(2,1), A(2,2) and A(2,3) in physical space.
        !===ru(:, 7:9 ,:)   is A(3,1), A(3,2) and A(3,3) in physical space.
        !===ru(:,10:12,:)   is B(1,1), B(1,2) and B(1,3) in physical space.
        !===ru(:,13:15,:)   is B(2,1), B(2,2) and B(2,3) in physical space.
        !===ru(:,16:18,:)   is B(3,1), B(3,2) and B(3,3) in physical space.
        !===Compute A : B
        IF (nb_field==6) THEN
            prod_ru(:, 1, :) = ru(:, 1, :) * ru(:, 10, :) + ru(:, 2, :) * ru(:, 11, :) + ru(:, 3, :) * ru(:, 12, :) &
                    + ru(:, 4, :) * ru(:, 13, :) + ru(:, 5, :) * ru(:, 14, :) + ru(:, 6, :) * ru(:, 15, :) &
                    + ru(:, 7, :) * ru(:, 16, :) + ru(:, 8, :) * ru(:, 17, :) + ru(:, 9, :) * ru(:, 18, :)
        END IF
        !===End Compute A : B
        IF(PRESENT(opt_rr)) THEN
            IF(.NOT.(PRESENT(opt_rank_S))) CALL error_petsc('Wrong arguments FFT_TENSOR_CONTRACTION_gauss')
            IF(PRESENT(opt_coef)) THEN
                !CALL WRITE_REAL_GAUSS_PDF(opt_rank_S, rank_F, prod_ru,  &
                !     np, opt_rr, 'Dnu', nb_procs, bloc_size, opt_number_co=opt_coef)
            ELSE
                !CALL WRITE_REAL_GAUSS_PDF(opt_rank_S, rank_F, prod_ru,  &
                !         np, opt_rr, 'eps', nb_procs, bloc_size)
            END IF
        END IF

        howmany = bloc_size * 2 / 2
        CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
                inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
        CALL dfftw_execute(fftw_plan_multi_r2c)
        CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

        prod_cu = prod_cu * (1.d0 / N_r_pad) !Scaling

        !Now we need to redistribute the Fourier coefficients on each processor
        combined_prod_cu(:, 1) = prod_cu(1, :)

        DO n = 2, m_max
            combined_prod_cu(:, n) = 2 * CONJG(prod_cu(n, :))
        END DO

        longueur_tranche = bloc_size * m_max_c * 2
        MPID = MPI_DOUBLE_PRECISION
        CALL MPI_ALLTOALL (combined_prod_cu, longueur_tranche, MPID, dist_prod_cu, longueur_tranche, &
                MPID, communicator, code)

        DO i = 1, m_max_c
            DO nb = 1, nb_procs
                shiftc = (nb - 1) * bloc_size
                shiftl = (nb - 1) * m_max_c
                intermediate = dist_prod_cu(:, shiftl + i)
                DO n = 1, bloc_size
                    IF (n + shiftc > np) CYCLE
                    S1_out(n + shiftc, 1, i) = REAL (intermediate(n), KIND = 8)
                    S1_out(n + shiftc, 2, i) = AIMAG(intermediate(n))
                END DO
            END DO
        END DO
    END SUBROUTINE FFT_TENSOR_CONTRACTION_gauss

    ! SUBROUTINE FFT_INCR_VEL(communicator, V1_in, V2_in, V3_in, V_out, nb_procs, bloc_size, m_max_pad, &
    !         opt_I, type_field, mesh, opt_write)
    !     USE pdf
    !     USE def_type_mesh
    !     !S1_out=FFT(MAX(FFT-1(grad(A)ij))) grad(A) being tensor
    !     !===V1_in is first  ligne of tensor A(1,1:np,1:6,list_modes)
    !     !===    d_i(v_r) i={r, theta, z}
    !     !===V2_in is second ligne of tensor A(2,1:np,1:6,list_modes)
    !     !===    d_i(v_th) i={r,theta, z}
    !     !===V3_in is third  ligne of tensor A(3,1:np,1:6,list_modes)
    !     !===    d_i(v_z) i={r, theta, z}
    !     !===V_out corresonds to three scalars : (1:np,1:6,list_mode)
    !     !=== 1:2 is max(d_i u_j) 3:4 is max(d_i u_j + d_j u_i)/2
    !     !=== 5:6 is max(d_i u_j- d_j u_i)/2
    !     !This a de-aliased version of the code, FEB 4, 2011, JLG
    !     IMPLICIT NONE
    !     !
    !     ! Format: V_1in(1:np,1:6,1:m_max_c)
    !     ! INPUT ARE COSINE AND SINE COEFFICIENTS
    !     ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    !     REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: V1_in, V2_in, V3_in ! VECTOR
    !     REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: V_out
    !     INTEGER, INTENT(IN) :: bloc_size, m_max_pad, nb_procs
    !     LOGICAL, OPTIONAL :: opt_write
    !     CHARACTER(*), INTENT(IN), OPTIONAL :: type_field
    !     CHARACTER(*), INTENT(IN), OPTIONAL :: opt_I
    !     TYPE(mesh_type), INTENT(IN), OPTIONAL :: mesh

    !     COMPLEX(KIND = 8), DIMENSION(m_max_pad, SIZE(V1_in, 2) * 3 / 2, bloc_size) :: cu
    !     REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, SIZE(V1_in, 2) * 3 / 2, bloc_size) :: ru
    !     REAL(KIND = 8), DIMENSION(SIZE(V1_in, 3), 3 * SIZE(V1_in, 2), bloc_size * nb_procs) :: dist_field, combined_field
    !     !    REAL(KIND=8), DIMENSION(2*m_max_pad-1,bloc_size)                     :: prod_ru
    !     REAL(KIND = 8), DIMENSION(2 * m_max_pad - 1, 3, bloc_size) :: prod_ru
    !     COMPLEX(KIND = 8), DIMENSION(m_max_pad, 3, bloc_size) :: prod_cu
    !     COMPLEX(KIND = 8), DIMENSION(3, bloc_size) :: intermediate
    !     COMPLEX(KIND = 8), DIMENSION(3, bloc_size, SIZE(V1_in, 3) * nb_procs) :: combined_prod_cu
    !     COMPLEX(KIND = 8), DIMENSION(3, bloc_size, SIZE(V1_in, 3) * nb_procs) :: dist_prod_cu

    !     REAL(KIND = 8), DIMENSION(3 * SIZE(V1_in, 2) / 2) :: temp_grad
    !     REAL(KIND = 8), DIMENSION(3 * SIZE(V1_in, 2) / 2) :: temp_sym_grad, temp_antisym_grad
    !     INTEGER :: np, np_tot, nb_field, m_max, m_max_c, MPID, N_r_pad
    !     INTEGER(KIND = 8) :: fftw_plan_multi_c2r, fftw_plan_multi_r2c
    !     INTEGER :: nb, nf, shiftc, shiftl, jindex, longueur_tranche, i, n, code, rank, l, i_field
    !     !REAL(KIND=8), DIMENSION(2*m_max_pad-1, bloc_size) ::  norm_vel_loc
    !     !REAL(KIND=8)    :: max_velocity
    !     ! FFTW parameters
    !     INTEGER :: fft_dim, howmany, istride, ostride, idist, odist
    !     INTEGER, DIMENSION(1) :: dim, inembed, onembed
    !     ! Recall complexes must be rescaled
    !     ! End FFTW parameters
    !     !#include "petsc/finclude/petsc.h"
    !     MPI_Comm, DIMENSION(:), POINTER :: communicator
    !     CALL MPI_COMM_RANK(communicator(2), rank, code)

    !     np = SIZE(V1_in, 1)
    !     nb_field = SIZE(V1_in, 2) ! Number of fields
    !     m_max_c = SIZE(V1_in, 3) ! Number of complex (cosines + sines) coefficients per point
    !     m_max = m_max_c * nb_procs! Number of comlex coefficients per point per processor
    !     N_r_pad = 2 * m_max_pad - 1
    !     np_tot = nb_procs * bloc_size

    !     IF (MOD(nb_field, 2)/=0 .OR. m_max_c==0) THEN
    !         WRITE(*, *) ' BUG in FFT_INCR_VEL '
    !         STOP
    !     END IF

    !     ! Bloc_size is the number of points that are handled by one processor
    !     ! once the Fourier modes are all collected
    !     ! Computation of bloc_size and np_tot
    !     ! fin de la repartition des points

    !     ! Packing all 3 complex components of both v1 and v2 input fields
    !     ! into dist_field, where the dimension indexing the nodal points varies the least rapidly,
    !     ! so that after distributing the data to the processes, each one will obtain a part
    !     ! on nodal points
    !     ! TRANSPOSE pr que la variable i associee aux modes soit la 1ere sur laquelle on va faire la FFT
    !     DO i = 1, m_max_c
    !         !       dist_field(i,1:nb_field,1:np) = TRANSPOSE(V1_in(:,:,i))
    !         dist_field(i, 1:nb_field, 1:np) = TRANSPOSE(V1_in(:, :, i))
    !         dist_field(i, nb_field + 1:2 * nb_field, 1:np) = TRANSPOSE(V2_in(:, :, i))
    !         dist_field(i, 2 * nb_field + 1:3 * nb_field, 1:np) = TRANSPOSE(V3_in(:, :, i))
    !     END DO

    !     IF (np/=np_tot) dist_field(:, :, np + 1:np_tot) = 0.d0

    !     longueur_tranche = bloc_size * m_max_c * nb_field * 3

    !     MPID = MPI_DOUBLE_PRECISION
    !     CALL MPI_ALLTOALL (dist_field, longueur_tranche, MPID, combined_field, longueur_tranche, &
    !             MPID, communicator(2), code)

    !     !JLG, FEB 4, 2011
    !     cu = 0.d0
    !     !JLG, FEB 4, 2011
    !     DO n = 1, bloc_size
    !         DO nb = 1, nb_procs
    !             shiftc = (nb - 1) * bloc_size
    !             shiftl = (nb - 1) * m_max_c
    !             jindex = n + shiftc
    !             DO nf = 1, nb_field * 3 / 2
    !                 ! Put real and imaginary parts in a complex
    !                 ! nf=1,2,3 => V1_in
    !                 ! nf=4,5,6 => V2_in
    !                 ! nf=7,8,9 => V3_in
    !                 ! INPUT ARE COSINE AND SINE COEFFICIENTS
    !                 ! THEY ARE PUT IN COMPLEX FORMAT: c_0 = a_0 + i*0 and c_n = (a_n-i*b_n)/2
    !                 cu(shiftl + 1:shiftl + m_max_c, nf, n) = 0.5d0 * CMPLX(combined_field(:, 2 * nf - 1, jindex), &
    !                         -combined_field(:, 2 * nf, jindex), KIND = 8)
    !             END DO
    !         END DO
    !     END DO
    !     cu(1, :, :) = 2 * CMPLX(REAL(cu(1, :, :), KIND = 8), 0.d0, KIND = 8)

    !     ! Set the parameters for dfftw
    !     fft_dim = 1; istride = 1; ostride = 1;
    !     idist = N_r_pad;   inembed(1) = N_r_pad; DIM(1) = N_r_pad
    !     odist = m_max_pad; onembed(1) = m_max_pad

    !     howmany = bloc_size * nb_field * 3 / 2

    !     CALL dfftw_plan_many_dft_c2r(fftw_plan_multi_c2r, fft_dim, dim, howmany, cu, &
    !             onembed, ostride, odist, ru, inembed, istride, idist, FFTW_ESTIMATE)
    !     CALL dfftw_execute(fftw_plan_multi_c2r)

    !     ! clean up
    !     CALL dfftw_destroy_plan(fftw_plan_multi_c2r)

    !     !===ru(:, 1:3 ,:)   is d_i(v_r)  i={r, theta, z} in physical space.
    !     !===ru(:, 4:6 ,:)   is d_i(v_th) i={r, theta, z} in physical space.
    !     !===ru(:, 7:9 ,:)   is d_i(v_z)  i={r, theta, z} in physical space.


    !     !HF here do max en fonction du cas symétrique ou non puis fft dans l'autre sens'
    !     !===Structure functions
    !     IF (nb_field==6) THEN
    !         DO i = 1, 2 * m_max_pad - 1
    !             DO n = 1, bloc_size
    !                 if (size(ru, 2) - size(temp_grad).ne. 0) call error_petsc('size pb in incr vel')
    !                 temp_grad = ru(i, :, n)
    !                 temp_sym_grad = temp_grad
    !                 temp_sym_grad(2) = ABS((temp_grad(2) + temp_grad(4)) / 2.d0)
    !                 temp_sym_grad(3) = ABS((temp_grad(3) + temp_grad(7)) / 2.d0)
    !                 temp_sym_grad(6) = ABS((temp_grad(6) + temp_grad(8)) / 2.d0)
    !                 temp_sym_grad(4) = ABS((temp_grad(2) + temp_grad(4)) / 2.d0)
    !                 temp_sym_grad(7) = ABS((temp_grad(3) + temp_grad(7)) / 2.d0)
    !                 temp_sym_grad(8) = ABS((temp_grad(6) + temp_grad(8)) / 2.d0)
    !                 temp_antisym_grad = 0.d0
    !                 temp_antisym_grad(2) = ABS((temp_grad(2) - temp_grad(4)) / 2.d0)
    !                 temp_antisym_grad(3) = ABS((temp_grad(3) - temp_grad(7)) / 2.d0)
    !                 temp_antisym_grad(6) = ABS((temp_grad(6) - temp_grad(8)) / 2.d0)
    !                 temp_antisym_grad(4) = ABS((temp_grad(2) - temp_grad(4)) / 2.d0)
    !                 temp_antisym_grad(7) = ABS((temp_grad(3) - temp_grad(7)) / 2.d0)
    !                 temp_antisym_grad(8) = ABS((temp_grad(6) - temp_grad(8)) / 2.d0)
    !                 !=== Don't forget absolute values on temp_grad
    !                 temp_grad = ABS(temp_grad(:))
    !                 !=== Don't forget absolute values on temp_grad
    !                 prod_ru(i, 1, n) = MAXVAL(temp_grad(:)) !max_{i,j}(d_i u_j)
    !                 prod_ru(i, 2, n) = MAXVAL(temp_sym_grad(:)) !max_{i,j}(d_i u_j + d_j u_i)/2
    !                 prod_ru(i, 3, n) = MAXVAL(temp_antisym_grad(:)) !max_{i,j}(d_i u_j - d_j u_i)/2
    !             END DO
    !         END DO
    !     ELSE
    !         CALL error_petsc('Wrong call of FFT_INCR_VEL')
    !     END IF

    !     IF (present(opt_write)) THEN
    !         IF (opt_write) THEN
    !             CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, prod_ru(:, 1, :), type_field, opt_I)
    !             CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, prod_ru(:, 2, :), type_field // 's', opt_I)
    !             CALL WRITE_REAL_GAUSS_PDF(communicator, mesh, prod_ru(:, 3, :), type_field // 'a', opt_I)
    !         END IF
    !     END IF

    !     howmany = bloc_size * nb_field / 2 !3 scalars out
    !     CALL dfftw_plan_many_dft_r2c(fftw_plan_multi_r2c, fft_dim, dim, howmany, prod_ru, &
    !             inembed, istride, idist, prod_cu, onembed, ostride, odist, FFTW_ESTIMATE)
    !     CALL dfftw_execute(fftw_plan_multi_r2c)
    !     CALL dfftw_destroy_plan(fftw_plan_multi_r2c)

    !     prod_cu = prod_cu * (1.d0 / N_r_pad) !Scaling

    !     !Now we need to redistribute the Fourier coefficients on each processor
    !     combined_prod_cu(:, :, 1) = prod_cu(1, :, :)

    !     DO n = 2, m_max
    !         combined_prod_cu(:, :, n) = 2 * CONJG(prod_cu(n, :, :))
    !     END DO

    !     longueur_tranche = bloc_size * m_max_c * nb_field
    !     MPID = MPI_DOUBLE_PRECISION
    !     CALL MPI_ALLTOALL (combined_prod_cu, longueur_tranche, MPID, dist_prod_cu, longueur_tranche, &
    !             MPID, communicator(2), code)

    !     DO i = 1, m_max_c
    !         DO nb = 1, nb_procs
    !             shiftc = (nb - 1) * bloc_size
    !             shiftl = (nb - 1) * m_max_c
    !             intermediate = dist_prod_cu(:, :, shiftl + i)
    !             DO n = 1, bloc_size
    !                 IF (n + shiftc > np) CYCLE
    !                 DO i_field = 1, nb_field / 2
    !                     V_out(n + shiftc, i_field * 2 - 1, i) = REAL (intermediate(i_field, n), KIND = 8)
    !                     V_out(n + shiftc, i_field * 2, i) = AIMAG(intermediate(i_field, n))
    !                 END DO
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE FFT_INCR_VEL

END MODULE fft_DR
