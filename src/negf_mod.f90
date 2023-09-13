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
module negf_mod
    !! Non-equilibrium Green's function (NEGF) module, upper-level driver for solving the NEGF equations


    implicit none

    private

    integer, parameter :: dp = 8

    public :: negf_solve

contains

    subroutine negf_solve(nx, Hii, H1i, Sii)
        use matrix_c, only: type_matrix_complex,malloc,free,sizeof
        use rgf_mod, only: rgf_variableblock_backward
        type(type_matrix_complex), intent(in), dimension(nx)::Hii,Sii
        type(type_matrix_complex), intent(in), dimension(nx+1)::H1i
        integer, intent(in)::nx  !! number of blocks        
        ! ----
        real(dp),allocatable,dimension(:)::temp !! temperatures
        real(dp),allocatable,dimension(:):: mu !! chemical potentials
        real(dp)::emin,emax !! min and max energy range        
        integer:: nen !! number of energy points
        integer::nm
        real(dp)::en
        real(dp),dimension(:,:),allocatable::mul,mur,templ,tempr
        type(type_matrix_complex),dimension(:),allocatable::sigma_lesser_ph,sigma_r_ph,G_r,G_lesser,G_greater,Jdens,Gl,Gln
        real(dp)::tr,tre
        en = 0.0d0 
        nm=size(Hii(1)%m,1)
        allocate(mul(nm,nm))
        allocate(templ(nm,nm))
        nm=size(Hii(nx)%m,1)
        allocate(mur(nm,nm))
        allocate(tempr(nm,nm))

        print *, sizeof(Hii)
        allocate(sigma_lesser_ph(nx))
        allocate(sigma_r_ph(nx))
        allocate(G_r(nx))
        allocate(G_lesser(nx))
        allocate(G_greater(nx))
        allocate(Jdens(nx))
        allocate(Gl(nx))
        allocate(Gln(nx))

        call malloc(sigma_lesser_ph,nx, sizeof(Hii))
        call malloc(sigma_r_ph,nx, sizeof(Hii))
        call malloc(G_r,nx, sizeof(Hii))
        call malloc(G_lesser,nx, sizeof(Hii))
        call malloc(G_greater,nx, sizeof(Hii))
        call malloc(Jdens,nx, sizeof(Hii))
        call malloc(Gl,nx, sizeof(Hii))
        call malloc(Gln,nx, sizeof(Hii))

        print *,'allocate memory done'
        call rgf_variableblock_backward(nx,En, mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph, &
            sigma_r_ph, G_r, G_lesser, G_greater, Jdens, Gl, Gln, tr, tre)

        print *,'free memory'
        deallocate(mul,mur,tempr,templ)
        call free(sigma_lesser_ph)
        call free(sigma_r_ph)
        call free(G_r)
        call free(G_lesser)
        call free(G_greater)
        call free(Jdens)
        call free(Gl)
        call free(Gln)

    end subroutine negf_solve

end module negf_mod
