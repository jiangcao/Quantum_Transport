! Copyright (c) 2023 Jiang Cao, ETH Zurich
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors
!    may be used to endorse or promote products derived from this software without
!    specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
module linalg

    use static

    implicit none

    private

    public :: invert, invert_banded, cross, eig, eigv, eigv_feast

CONTAINS

    ! matrix inversion
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

    ! find the inverse of a banded matrix A by solving a system of linear equations
    !   on exit, A contains the banded matrix of inv(A)
    !   banded format see [https://netlib.org/lapack/lug/node124.html]
    subroutine invert_banded(A, nn, nb)
        integer, intent(in)::nn, nb
        complex(8), intent(inout)::A(3*nb + 1, nn)
        complex(8), allocatable::work(:), B(:, :), X(:, :)
        integer, allocatable::ipiv(:)
        integer::info, lda, lwork, ldb, i, nrhs
        allocate (ipiv(nn))
        allocate (work(nn*nn))
        lda = 3*nb + 1
        call zgbtrf(nn, nn, nb, nb, A, lda, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgbtrf failed, info=', info
            call abort()
        end if
        ldb = 1
        allocate (B(ldb, nn))
        allocate (X(lda, nn))
        nrhs = ldb
        do i = 1, nn
            B = 0.0d0
            B(1, i) = 1.0d0
            call zgbtrs('N', nn, nb, nb, nrhs, A, lda, ipiv, B, ldb, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgbtrs failed, info=', info
                call abort()
            end if
            X(1:nb*2 + 1, i) = B(1, i - nb:i + nb)
        end do
        A = X
        deallocate (B, work, ipiv, X)
    end subroutine invert_banded

    ! vector cross-product
    FUNCTION cross(a, b)
        REAL(8), DIMENSION(3) :: cross
        REAL(8), DIMENSION(3), INTENT(IN) :: a, b
        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)
    END FUNCTION cross

    ! calculate eigen-values of a Hermitian matrix A
    FUNCTION eig(NN, A)
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:, :) :: A
        ! -----
        REAL(8) :: eig(NN)
        real(8) :: W(1:NN)
        integer :: INFO, LWORK, liwork, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:)
        lwork = max(1, 2*NN - 1)
        lrwork = max(1, 3*NN - 2)
        allocate (work(lwork))
        allocate (rwork(lrwork))
        !
        CALL zheev('N', 'U', NN, A, NN, W, WORK, LWORK, RWORK, INFO)
        !
        deallocate (work, rwork)
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: ZHEEV HAS FAILED. INFO=', INFO
            call abort()
        end if
        eig(:) = W(:)
    END FUNCTION eig

    ! calculate all eigen-values and eigen-vectors of a Hermitian matrix A
    !   within a given search interval, a wrapper to the FEAST function in MKL https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/feast-syev-feast-heev.html
    !   upon return A(:,1:m) will be modified and contains the eigen-vectors
    FUNCTION eigv_feast(NN, A, emin, emax, m)
        include 'mkl.fi'
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:, :) :: A
        REAL(8), INTENT(IN) :: emin, emax ! lower and upper bounds of the interval to be searched for eigenvalues
        REAL(8) :: eigv_feast(NN)
        integer, intent(out) :: m ! total number of eigenvalues found
        ! -----
        real(8) :: epsout
        integer :: fpm(128), m0, loop, info
        complex(8), allocatable :: x(:, :)
        real(8), allocatable :: w(:), res(:)
        m0 = max(nn/10, 1000)
        allocate (x(nn, m0))
        allocate (w(m0))
        allocate (res(m0))
        !
        call feastinit(fpm)
        fpm(1) = 1 ! print runtime status to the screen
        !
        call zfeast_heev('U', nn, A, nn, fpm, epsout, loop, emin, emax, m0, W, x, m, res, info)
        !
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: zfeast_heev HAS FAILED. INFO=', INFO
            call abort()
        end if
        eigv_feast(1:m) = W(1:m)
        A(:, 1:m) = x(:, 1:m)
        deallocate (x, w, res)
    END FUNCTION eigv_feast

    ! calculate eigen-values and eigen-vectors of a Hermitian matrix A
    !   upon return A will be modified and contains the eigen-vectors
    FUNCTION eigv(NN, A)
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:, :) :: A
        ! -----
        REAL(8) :: eigv(NN)
        real(8) :: W(1:NN)
        integer :: INFO, LWORK, liwork, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:)
        lwork = max(1, 2*NN - 1)
        lrwork = max(1, 3*NN - 2)
        allocate (work(lwork))
        allocate (rwork(lrwork))
        !
        CALL zheev('V', 'U', NN, A, NN, W, WORK, LWORK, RWORK, INFO)
        !
        deallocate (work, rwork)
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: ZHEEV HAS FAILED. INFO=', INFO
            call abort
        end if
        eigv(:) = W(:)
    END FUNCTION eigv

end module linalg
