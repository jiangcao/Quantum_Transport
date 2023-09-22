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
! AUTHOR: Jiang Cao
!
module cuda_rgf_mod
    !! Recursive  Green's  function  solvers module

    implicit none

    private

    integer, parameter :: dp = 8
    REAL(dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1

    public :: sancho

contains



!!  Sancho-Rubio
subroutine sancho(nm, E, S00, H00, H10, G00, GBB)
    use cublas
    use linalg, only: invert
    integer i, j, k, nm, nmax
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
    do i = 1, nmax
        A = z*S00 - H_BB
        call invert(A, nm)
        call Zgemm('n', 'n', nm, nm, nm, alpha, A, nm, H_10, nm, beta, B, nm)
        call Zgemm('n', 'n', nm, nm, nm, alpha, H_01, nm, B, nm, beta, C, nm)
        H_SS = H_SS + C
        H_BB = H_BB + C
        call Zgemm('n', 'n', nm, nm, nm, alpha, H_10, nm, B, nm, beta, C, nm)
        call Zgemm('n', 'n', nm, nm, nm, alpha, A, nm, H_01, nm, beta, B, nm)
        call Zgemm('n', 'n', nm, nm, nm, alpha, H_10, nm, B, nm, beta, A, nm)
        H_10 = C
        H_BB = H_BB + A
        call Zgemm('n', 'n', nm, nm, nm, alpha, H_01, nm, B, nm, beta, C, nm)
        H_01 = C
        ! NORM --> inspect the diagonal of A
        error = 0.0d0
        DO k = 1, nm
            DO j = 1, nm
                error = error + sqrt(aimag(C(k, j))**2 + Dble(C(k, j))**2)
            END DO
        END DO
        !write(90,*)E,i,error
        tmp = H_SS
        IF (abs(error) < TOL) THEN
            !write(90,*) 'SR: Exited, abs(error)=',i,abs(error)
            EXIT
        ELSE
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
end subroutine sancho


end module cuda_rgf_mod