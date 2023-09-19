module cuda_rgf_mod 

use cublas
use cusolverdn

IMPLICIT NONE

    type type_matrix_complex
        complex(8), allocatable :: m(:, :)
    !! complex matrix
        integer :: size(2)
    !! matrix size
    end type type_matrix_complex

private
public :: cuda_rgf_variableblock_forward 

integer,parameter::dp=8
complex(dp),parameter:: czero=dcmplx(0.0d0,0.0d0)
REAL(dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1

contains

!!  Recursive Forward Green's solver
    subroutine cuda_rgf_variableblock_forward(nx, En, mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph, &
                                         sigma_r_ph, G_r, G_lesser, G_greater, Jdens, Gl, Gln, tr, tre)
        type(type_matrix_complex), intent(in) :: Hii(nx), H1i(nx + 1), Sii(nx), sigma_lesser_ph(nx), sigma_r_ph(nx)
        real(dp), intent(in)       :: En, mul(:, :), mur(:, :), TEMPr(:, :), TEMPl(:, :)
        integer, intent(in) :: nx !! lenght of the device
        type(type_matrix_complex), intent(inout):: G_greater(nx), G_lesser(nx), G_r(nx), Jdens(nx), Gl(nx), Gln(nx)
        real(dp), intent(out)      :: tr, tre
!---- local variables
        integer    :: M, ii, jj
        complex(dp) :: z
        real(dp)    :: tim
        complex(dp), allocatable :: sig(:, :), H00(:, :), H10(:, :)
        complex(dp), allocatable :: A(:, :), B(:, :), C(:, :), G00(:, :), GBB(:, :), sigmar(:, :), sigmal(:, :), GN0(:, :)
!!!$omp declare target device_type(any)
        z = dcmplx(En, 0.0d0)
!
! on the left contact
        ii = 1
        M = size(Hii(ii)%m, 1)
        allocate (H00(M, M))
        allocate (H10(M, M))
        allocate (G00(M, M))
        allocate (GBB(M, M))
        allocate (sigmal(M, M))
        allocate (sig(M, M))
!
!! $$H00 = H(i,i) + \Sigma_{ph}(i) * S(i,i)$$
        call MUL_c(sigma_r_ph(ii)%m, Sii(ii)%m, 'n', 'n', B)
!
        H00 = Hii(ii)%m + B
        H10 = H1i(ii)%m
        call sancho(M, En, Sii(ii)%m, H00, transpose(conjg(H10)), G00, GBB)

!!!$omp critical
!!        open (unit=10, file='sancho_g00.dat', position='append')
!!        write (10, *) En, 2, -aimag(trace(G00))
!!        close (10)
!!        open (unit=10, file='sancho_gbb.dat', position='append')
!!        write (10, *) En, 2, -aimag(trace(Gbb))
!!        close (10)
!!!$omp end critical

!! $$\Sigma^R = H(i,i+1) * G00 * H(i+1,i)$$
!! $$Gl(i) = [E*S(i,i) - H00 - \Sigma_R]^{-1}$$
        call triMUL_c(H10, G00, H10, sigmal, 'n', 'n', 'c')
        B = z*Sii(ii)%m - H00 - sigmal
        call invert(B, M)
        Gl(ii)%m = B
!
!! $$Gln(i) = Gl(i) * [\Sigma_{ph}^<(i)*S(i,i) + (-(\Sigma^R - \Sigma_R^\dagger)*ferm(..))] * Gl(i)^\dagger$$
        call MUL_c(sigma_lesser_ph(ii)%m, Sii(ii)%m, 'n', 'n', B)
        sig = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mur)/(BOLTZ*TEMPr))
!
        sig = sig + B
        call triMUL_c(Gl(ii)%m, sig, Gl(ii)%m, B, 'n', 'n', 'c')
        Gln(ii)%m = B
        deallocate (G00, GBB, sig, H10)
!
        allocate (A(M, M))
! inside device l -> r
        do ii = 2, nx - 1
            M = size(Hii(ii)%m, 1)
            if (size(H00, 1) .ne. M) then
                deallocate (H00, A)
                allocate (H00(M, M))
                allocate (A(M, M))
            end if
            call MUL_c(sigma_r_ph(ii)%m, Sii(ii)%m, 'n', 'n', B)
            H00 = Hii(ii)%m + B
!
!! $$H00 = H(i,i) + \Sigma_{ph}(i) * S(i,i)$$
!! $$Gl(i) = [E*S(i,i) - H00 - H(i,i-1) * Gl(i-1) * H(i-1,i)]^{-1}$$
            call triMUL_c(H1i(ii)%m, Gl(ii - 1)%m, H1i(ii)%m, B, 'n', 'n', 'c')
            A = z*Sii(ii)%m - H00 - B
            call invert(A, M)
            Gl(ii)%m = A
!
!! $$Gln(i) = Gl(i) * [\Sigma_{ph}^<(i)*S(i,i) + H(i,i+1)*Gln(i+1)*H(i+1,i)] * Gl(i)^\dagger$$
            call triMUL_c(H1i(ii)%m, Gln(ii - 1)%m, H1i(ii)%m, B, 'n', 'n', 'c')
            call MUL_c(sigma_lesser_ph(ii)%m, Sii(ii)%m, 'n', 'n', A)
            B = B + A
            call triMUL_c(Gl(ii)%m, B, Gl(ii)%m, A, 'n', 'n', 'c')
            Gln(ii)%m = A
        end do
!
! on the right contact
        ii = nx
        M = size(Hii(ii)%m, 1)
        allocate (H10(M, M))
        allocate (G00(M, M))
        allocate (GBB(M, M))
        allocate (sig(M, M))
        allocate (sigmar(M, M))
        if (size(H00, 1) .ne. M) then
            deallocate (H00)
            allocate (H00(M, M))
        end if
!
        call MUL_c(sigma_r_ph(ii)%m, Sii(ii)%m, 'n', 'n', B)
        H00 = Hii(ii)%m + B
        H10 = H1i(nx + 1)%m
!
        call sancho(M, En, Sii(ii)%m, H00, H10, G00, GBB)
!
        call triMUL_c(H10, G00, H10, sigmar, 'c', 'n', 'n')

!!!$omp critical
!!        open (unit=10, file='sancho_g00.dat', position='append')
!!        write (10, *) En, 1, -aimag(trace(G00))
!!        close (10)
!!        open (unit=10, file='sancho_gbb.dat', position='append')
!!        write (10, *) En, 1, -aimag(trace(Gbb))
!!        close (10)
!!!$omp end critical

        call triMUL_c(H1i(nx)%m, Gl(nx - 1)%m, H1i(nx)%m, B, 'n', 'n', 'c')
        A = z*Sii(ii)%m - H00 - B - sigmar
!
        call invert(A, M)
        G_r(ii)%m = A
        Gl(ii)%m = A
!
!! $$\Sigma^< = \Sigma_11^< + \Sigma_{ph}^< + \Sigma_s^<$$
        call triMUL_c(H1i(nx)%m, Gln(nx - 1)%m, H1i(nx)%m, B, 'n', 'n', 'c')
        call MUL_c(sigma_lesser_ph(nx)%m, Sii(nx)%m, 'n', 'n', A)
        sig = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mul)/(BOLTZ*TEMPl))
        sig = sig + A + B
!
!! $$G^< = G * \Sigma^< * G^\dagger$$
        call triMUL_c(G_r(ii)%m, sig, G_r(ii)%m, B, 'n', 'n', 'c')
!
        G_lesser(ii)%m = B
        G_greater(ii)%m = G_lesser(ii)%m + (G_r(ii)%m - transpose(conjg(G_r(ii)%m)))
!
        A = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mul)/(BOLTZ*TEMPl))
        call MUL_c(A, G_greater(ii)%m, 'n', 'n', B)
        A = -(sigmar - transpose(conjg(sigmar)))*(ferm((En - mul)/(BOLTZ*TEMPl)) - 1.0d0)
        call MUL_c(A, G_lesser(ii)%m, 'n', 'n', C)
!
        Jdens(ii)%m = B - C
!
        tim = 0.0d0
        do jj = 1, M
            tim = tim + dble(Jdens(ii)%m(jj, jj))
        end do
        tr = tim ! transmission
        deallocate (sigmar, sig, G00, GBB, H10)
        allocate (GN0(M, M))
!
! inside device r -> l
        do ii = nx - 1, 1, -1
            M = size(Hii(ii)%m, 1)
!! $$A = G^<(i+1) * H(i+1,i) * Gl(i)^\dagger + G(i+1) * H(i+1,i) * Gln(i)$$
            call triMUL_c(G_lesser(ii + 1)%m, H1i(ii)%m, Gl(ii)%m, A, 'n', 'n', 'c')
            call triMUL_c(G_r(ii + 1)%m, H1i(ii)%m, Gln(ii)%m, B, 'n', 'n', 'n')
            A = A + B
!! $$B = H(i,i+1) * A$$
!! $$Jdens(i) = -2 * B$$
            call MUL_c(H1i(ii)%m, A, 'c', 'n', B)
            Jdens(ii)%m = -2.0d0*B(:, :)
!
!! $$GN0 = Gl(i) * H(i,i+1) * G(i+1)$$
!! $$G(i) = Gl(i) + GN0 * H(i+1,i) * Gl(i)$$
            call MUL_c(Gl(ii)%m, H1i(ii)%m, 'n', 'c', B)
            call MUL_c(B, G_r(ii + 1)%m, 'n', 'n', GN0)
            call MUL_c(GN0, H1i(ii)%m, 'n', 'n', C)
            call MUL_c(C, Gl(ii)%m, 'n', 'n', A)
            G_r(ii)%m = Gl(ii)%m + A
!
!! $$G^<(i) = Gln(i) + Gl(i) * H(i,i+1) * G^<(i+1) * H(i+1,i) *Gl(i)^\dagger$$
            call MUL_c(Gl(ii)%m, H1i(ii)%m, 'n', 'c', B)
            call MUL_c(B, G_lesser(ii + 1)%m, 'n', 'n', C)
            call MUL_c(C, H1i(ii)%m, 'n', 'n', A)
            call MUL_c(A, Gl(ii)%m, 'n', 'c', C)
            G_lesser(ii)%m = Gln(ii)%m + C
!
!! $$G^<(i) = G^<(i) + GN0 * H(i+1,i) * Gln(i)$$
            call MUL_c(GN0, H1i(ii)%m, 'n', 'n', B)
            call MUL_c(B, Gln(ii)%m, 'n', 'n', C)
            G_lesser(ii)%m = G_lesser(ii)%m + C
!
!! $$G^<(i) = G^<(i) + Gln(i) * H(i,i+1) * GN0$$
            call MUL_c(Gln(ii)%m, H1i(ii)%m, 'n', 'c', B)
            call MUL_c(B, GN0, 'n', 'c', C)
            G_lesser(ii)%m = G_lesser(ii)%m + C
!
!! $$G^>(i) = G^<(i) + [G(i) - G(i)^\dagger]$$
            G_greater(ii)%m = G_lesser(ii)%m + (G_r(ii)%m - transpose(conjg(G_r(ii)%m)))
        end do
        ii = 1
! on the left contact
        A = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mur)/(BOLTZ*TEMPr))
        call MUL_c(A, G_greater(ii)%m, 'n', 'n', B)
        A = -(sigmal - transpose(conjg(sigmal)))*(ferm((En - mur)/(BOLTZ*TEMPr)) - 1.0d0)
        call MUL_c(A, G_lesser(ii)%m, 'n', 'n', C)
        tim = 0.0d0
        do jj = 1, M
            tim = tim + dble(B(jj, jj) - C(jj, jj))
        end do
        tre = tim
        deallocate (B, A, C, GN0, sigmal)
    
    end subroutine cuda_rgf_variableblock_forward


    subroutine invert(A, nn)
        integer :: info, nn
        integer, dimension(:), allocatable :: ipiv
        complex(8), dimension(nn, nn), intent(inout) :: A
        complex(8), dimension(:), allocatable :: work
        allocate (work(nn*nn))
        allocate (ipiv(nn))
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            A = czero
        else
            call zgetri(nn, A, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info
                A = czero
            end if
        end if
        deallocate (work)
        deallocate (ipiv)
    end subroutine invert

!!  Sancho-Rubio
    subroutine sancho(nm, E, S00, H00, H10, G00, GBB)
        integer i, j, k, nm, nmax, ii,jj
        COMPLEX(dp) :: z
        real(dp) :: E, error
        REAL(dp) :: TOL = 1.0D-100  ! [eV]
        COMPLEX(dp), INTENT(IN) ::  S00(nm, nm), H00(nm, nm), H10(nm, nm)
        COMPLEX(dp), INTENT(OUT) :: G00(nm, nm), GBB(nm, nm)
        COMPLEX(dp), ALLOCATABLE :: A(:, :), B(:, :), C(:, :), tmp(:, :)
        COMPLEX(dp), ALLOCATABLE :: H_BB(:, :), H_SS(:, :), H_01(:, :), H_10(:, :), Id(:, :)
        COMPLEX(dp), EXTERNAL :: ZLANGE
        complex(dp), parameter :: alpha = cmplx(1.0d0, 0.0d0)
        complex(dp), parameter :: beta = cmplx(0.0d0, 0.0d0)
        character*1 :: transa, transb
        complex(dp), dimension(:,:), allocatable :: work   
        integer, dimension(:), allocatable :: ipiv     
        integer :: info
        !
        Allocate (H_BB(nm, nm))
        Allocate (H_SS(nm, nm))
        Allocate (H_01(nm, nm))
        Allocate (H_10(nm, nm))
        Allocate (Id(nm, nm))
        Allocate (A(nm, nm))
        Allocate (B(nm, nm))
        Allocate (C(nm, nm))
        Allocate (tmp(nm, nm))
        Allocate (work(nm, nm))
        Allocate (ipiv(nm))
        nmax = 50
        z = dcmplx(E, 1.0d-5)
        Id = dcmplx(0.0d0,0.0d0)
        tmp = 0.0d0
        do i = 1, nm
            Id(i, i) = 1.0d0
            tmp(i, i) = dcmplx(0.0d0, 1.0d0)
        end do
        H_BB = H00
        H_10 = H10
        H_01 = TRANSPOSE(CONJG(H_10))
        H_SS = H00
        transa='N'
        transb='N'
        
        do i = 1, nmax            
            !$omp target enter data map(to:C,A,H_10,H_01,B,H_SS,H_BB,S00,ipiv,work,Id)  
            !$omp target data use_device_ptr(C,A,H_10,H_01,B,H_SS,H_BB,S00,ipiv,work,Id)
            ! A = z*S00 - H_BB            
            call zscal(nm*nm,beta,A,1)
            call zaxpy(nm*nm,-alpha,H_BB,1,A,1)
            call zaxpy(nm*nm,z,S00,1,A,1)

            ! inv(A)
            call zgetrf(nm, nm, A, nm, ipiv, info)
            call zcopy(nm*nm,Id,1,work,1)
            call zgetrs(transa,nm, nm, A, nm, ipiv, work, nm, info)
            call zcopy(nm*nm,work,1,A,1)
            
            call Zgemm(transa, transb, nm, nm, nm, alpha, A, nm, H_10, nm, beta, B, nm)                                    
            call Zgemm(transa, transb, nm, nm, nm, alpha, H_01, nm, B, nm, beta, C, nm)
            
            call zaxpy(nm*nm,alpha,C,1,H_SS,1)
            call zaxpy(nm*nm,alpha,C,1,H_BB,1)                                
            
            call Zgemm(transa, transb, nm, nm, nm, alpha, H_10, nm, B, nm, beta, C, nm)
            call Zgemm(transa, transb, nm, nm, nm, alpha, A, nm, H_01, nm, beta, B, nm)
            call Zgemm(transa, transb, nm, nm, nm, alpha, H_10, nm, B, nm, beta, A, nm)                        
            
            ! H_10 = C
            ! H_BB = H_BB + A
            call zcopy(nm*nm,C,1,H_10,1)
            call zaxpy(nm*nm,alpha,A,1,H_BB,1)                                
            call Zgemm(transa, transb, nm, nm, nm, alpha, H_01, nm, B, nm, beta, C, nm)            
            
            ! H_01 = C
            call zcopy(nm*nm,C,1,H_01,1)

            !$omp end target data 
            !$omp target update from(C,H_10,H_01,H_SS,H_BB)     
            !$omp target exit data map(delete:C,A,H_10,H_01,B,H_SS,H_BB,S00,work,ipiv,Id)

            ! NORM --> inspect the diagonal of A                                    
            error=0.0d0
            DO k = 1, nm
                DO j = 1, nm
                    error = error + sqrt(aimag(C(k, j))**2 + Dble(C(k, j))**2)
                END DO
            END DO
            !write(90,*)E,i,error
            tmp = H_SS
            ! call zcopy(nm*nm,H_SS,1,tmp,1)

            IF (abs(error) < TOL) THEN
                !write(90,*) 'SR: Exited, abs(error)=',i,abs(error)
                EXIT            
            END IF
            IF (i .EQ. nmax) THEN
                write (*, *) 'SEVERE warning: nmax reached in sancho!!!', error
            END IF
        end do

        G00 = z*S00 - H_SS
        
        call invert(G00, nm)
        !
        GBB = z*S00 - H_BB
        
        call invert(GBB, nm)
        !
        Deallocate (tmp)
        Deallocate (A)
        Deallocate (B)
        Deallocate (C)
        Deallocate (H_BB)
        Deallocate (H_SS)
        Deallocate (H_01)
        Deallocate (H_10)
        Deallocate (Id)
        Deallocate (work)
        Deallocate (ipiv)
    end subroutine sancho

    subroutine triMUL_C(A, B, C, R, trA, trB, trC)
        complex(8), intent(in), dimension(:, :) :: A, B, C
        complex(8), intent(inout), allocatable :: R(:, :)
        character, intent(in) :: trA, trB, trC
        complex(8), allocatable, dimension(:, :) :: tmp
        integer :: n, m, k, kb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trB is wrong: ", trB
            call abort()
        end if
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in triMUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        call MUL_C(A, B, trA, trB, tmp)
        call MUL_C(tmp, C, 'n', trC, R)
        deallocate (tmp)
    end subroutine triMUL_C

    subroutine MUL_C(A, B, trA, trB, R)
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8), intent(inout), allocatable :: R(:, :)
        CHARACTER, intent(in) :: trA, trB
        integer :: n, m, k, kb, lda, ldb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trB is wrong: ", trB
            call abort()
        end if
        lda = size(A, 1)
        ldb = size(B, 1)
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in MUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        if (allocated(R)) then
            if ((size(R, 1) .ne. m) .or. (size(R, 2) .ne. n)) then
                deallocate (R)
                Allocate (R(m, n))
            end if
        else
            Allocate (R(m, n))
        end if
        R = dcmplx(0.0d0, 0.0d0)
        call zgemm(trA, trB, m, n, k, dcmplx(1.0d0, 0.0d0), A, lda, B, ldb, dcmplx(0.0d0, 0.0d0), R, m)
    end subroutine MUL_C

!!  Fermi distribution function
    elemental Function ferm(a)
        Real(dp), intent(in) ::  a
        real(dp) :: ferm
        ferm = 1.0d0/(1.0d0 + Exp(a))
    End Function ferm



end module cuda_rgf_mod
