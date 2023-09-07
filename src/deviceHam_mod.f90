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
module deviceHam_mod

    implicit none

    private

    integer, parameter :: dp = 8

    public :: deviceHam_load_COOmatrix, deviceHam_build_blocks

contains

    subroutine deviceHam_load_COOmatrix(fname, H, nnz, nm, row, col, use0index, complex)
        character(len=*), intent(in)        :: fname !! input text file name
        complex(dp), allocatable, intent(out), dimension(:) :: H
        integer, allocatable, intent(out), dimension(:):: row, col
        integer, intent(out)::nnz
        integer, intent(out)::nm
        logical, intent(in), optional::use0index, complex
        logical :: l0index
        real(8) :: re, im
        integer::handle, io
        integer::M, NL, i, j, k
        handle = 101
        l0index = .false.
        if (present(use0index)) l0index = use0index
        open (unit=handle, file=fname)
        M = -HUGE(1)
        NL = 0
! Read through the file first to know the number of lines
        do
            read (handle, *, IOSTAT=IO) i, j, re
            if (IO < 0) exit
            if (max(i, j) > M) M = max(i, j)
            NL = NL + 1
        end do
        print '("Number of nonzeros = ",i18)', M
        allocate (H(M))
        allocate (row(M))
        allocate (col(M))
        H = 0.0d0
        nnz = NL
        nm = M
! Read again the file for the matrix
        rewind handle
        im = 0.0d0
        do k = 1, NL
            if ((present(complex)) .and. complex) then
                read (handle, *) i, j, re, im
            else
                read (handle, *) i, j, re
            end if
            if (l0index) then
                i = i + 1
                j = j + 1
            end if
            row(k) = i
            col(k) = j
            H(k) = dcmplx(re, im)
        end do
        close (handle)
    end subroutine deviceHam_load_COOmatrix

    subroutine deviceHam_build_blocks()
        use matrix_c, only: type_matrix_complex

    end subroutine deviceHam_build_blocks

end module deviceHam_mod
