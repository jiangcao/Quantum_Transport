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

    subroutine negf_solve(nx,nen,nk,emin,emax, Hii, H1i, Sii,temp,mu)
        use matrix_c, only: type_matrix_complex,malloc,free,sizeof
        use rgf_mod, only: rgf_variableblock_backward
        type(type_matrix_complex), intent(in), dimension(nx,nk)::Hii,Sii
        type(type_matrix_complex), intent(in), dimension(nx+1,nk)::H1i
        integer, intent(in)::nx  !! number of blocks        
        integer, intent(in)::nen  !! number of energy points        
        integer, intent(in)::nk  !! number of k points        
        real(dp),intent(in),dimension(2)::temp !! temperatures
        real(dp),intent(in),dimension(2):: mu !! chemical potentials
        ! ----
        real(dp)::emin,emax !! min and max energy range        
        integer::nm(2,nx), ik,ie
        real(dp)::en,dE
        real(dp),dimension(:,:),allocatable::mul,mur,templ,tempr
        type(type_matrix_complex),dimension(nx,nen,nk)::sigma_lesser_ph,sigma_r_ph,G_r,G_lesser,G_greater
        type(type_matrix_complex),dimension(nx)::Jdens,Gl,Gln
        real(dp)::tr,tre
        en = 0.0d0 
        nm=sizeof(Hii(:,1))
        allocate(mul(nm(1,1),nm(1,1)))
        allocate(templ(nm(1,1),nm(1,1)))
        allocate(mur(nm(1,nx),nm(1,nx)))
        allocate(tempr(nm(1,nx),nm(1,nx)))

        do ik=1,nk
            do ie=1,nen
                call malloc(sigma_lesser_ph(:,ie,ik),nx, nm)
                call malloc(sigma_r_ph(:,ie,ik),nx, nm)
                call malloc(G_r(:,ie,ik),nx, nm)
                call malloc(G_lesser(:,ie,ik),nx, nm)
                call malloc(G_greater(:,ie,ik),nx, nm)
            enddo
        enddo
        call malloc(Jdens,nx, nm)
        call malloc(Gl,nx, nm)
        call malloc(Gln,nx, nm)

        print *,'allocate memory done'

        dE = (emax-emin)/dble(nen-1)

        do ik=1,nk
            !$omp parallel default(shared) private(ie,en)
            !$omp do 
            do ie=1,nen
                En=emin+dble(ie-1)*dE
                call rgf_variableblock_backward(nx,En, mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph(:,ie,ik), &
                    sigma_r_ph(:,ie,ik), G_r(:,ie,ik), G_lesser(:,ie,ik), G_greater(:,ie,ik), Jdens, Gl, Gln, tr, tre)
            enddo
            !$omp end do 
            !$omp end parallel
        enddo

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
