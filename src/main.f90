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
PROGRAM main

    use negf_mod, only: negf_solve
    use matrix_c, only: type_matrix_complex, free
    use deviceHam_mod, only: devH_build_fromWannierFile

    implicit none

    type(type_matrix_complex), allocatable, dimension(:)::Hii, H1i, Sii
    integer::nx, ns

    ! MPI variables
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) comm_size
    integer ( kind = 4 ) comm_rank
    integer ( kind = 4 ) local_Nkz
    integer ( kind = 4 ) local_Nky
    integer ( kind = 4 ) local_NE
    integer ( kind = 4 ) first_local_energy


    include "mpif.h"
    call MPI_Init(ierr)

    nx = 5
    ns = 3

    call devH_build_fromWannierFile('ham_dat', Hii, H1i, Sii, nx, ns)
    print *,"solve"
    call negf_solve(nx, Hii, H1i, Sii)

    call free(Hii)
    call free(H1i)
    call free(Sii)

END PROGRAM main
