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

    public :: devH_build_fromCOOfile, devH_build_fromWannierFile

contains

    subroutine load_COOmatrix(fname, H, nnz, nm, row, col, use0index, complex)
        character(len=*), intent(in)        :: fname !! input text file name
        complex(dp), allocatable, intent(out), dimension(:) :: H
        integer, allocatable, intent(out), dimension(:):: row, col
        integer, intent(out)::nnz
        integer, intent(out)::nm
        logical, intent(in), optional::use0index, complex
        logical :: l0index
        real(dp) :: re, im
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
        print '("Number of nonzeros = ",i18)', NL
        allocate (H(NL))
        allocate (row(NL))
        allocate (col(NL))
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
    end subroutine load_COOmatrix

    subroutine devH_build_fromCOOfile(fname, Hii, H1i, ext_left, ext_right, num_slice, use0index, complex, threshold)
        use matrix_c, only: type_matrix_complex
        character(len=*), intent(in)        :: fname !! input text file name
        integer, intent(in)::ext_left, ext_right!! extension on left/right side
        integer, intent(out):: num_slice !!number of slices in the central part
        type(type_matrix_complex), dimension(:), intent(inout), allocatable::Hii, H1i !! Hamiltonian blocks
        logical, intent(in), optional::use0index, complex
        real(dp), intent(in), optional::threshold
        ! ----
        integer::nx, nnz, nm
        complex(dp), allocatable, dimension(:)::H !! Hamiltonian matrix value in COO
        integer, allocatable, dimension(:)::row, col !! Hamiltonian matrix index in COO
        integer, allocatable, dimension(:, :)::Slices !! slicing information , refer to [[graph_partition]]
        call load_COOmatrix(fname, H, nnz, nm, row, col, use0index, complex)
        allocate (Hii(nx))
        allocate (H1i(nx - 1))
        num_slice = size(Slices, 2)
        deallocate (H, col, row)
    end subroutine devH_build_fromCOOfile

    subroutine devH_build_fromWannierFile(fname, Hii, H1i, Sii, nx, nslab,nk,k, lreorder_axis, axis)
        use matrix_c, only: type_matrix_complex, malloc
        use wannierHam3d, only: w90_load_from_file, w90_free_memory, w90_MAT_DEF, nb
        character(len=*), intent(in)        :: fname !! input text file name
        logical, intent(in), optional :: lreorder_axis !! whether to reorder axis
        integer, intent(in), optional :: axis(3) !! permutation order
        integer, intent(in)::nx, nslab !! extension on left/right side
        type(type_matrix_complex), dimension(:,:), intent(inout), allocatable::Hii, H1i !! Hamiltonian blocks
        type(type_matrix_complex), dimension(:,:), intent(inout), allocatable::Sii !! overlap matrix blocks
        real(dp),intent(in)::k(2,nk)
        integer,intent(in)::nk
        ! ----
        complex(dp), allocatable, dimension(:, :)::H00, H10
        real(dp)::kx, ky, kz
        integer::nm, i, im,ik
        integer, dimension(nx+1)::nmm
        open (unit=10, file=trim(fname), status='unknown')
        call w90_load_from_file(10, lreorder_axis, axis)
        close (10)
        nm = nb*nslab
        nmm(:) = nm
        allocate (H00(nm, nm))
        allocate (H10(nm, nm))
        kx = 0.0d0

        allocate(Hii(nx,nk))
        allocate(H1i(nx+1,nk))
        allocate(Sii(nx,nk))
        do ik=1,nk
            ky=k(1,ik)
            kz=k(2,ik)
            call w90_MAT_DEF(H00, H10, kx, ky, kz, nslab)
            call w90_free_memory()
    
            call malloc(Hii(:,ik), nx, nmm(1:nx))        
            call malloc(Sii(:,ik), nx, nmm(1:nx))
            call malloc(H1i(:,ik), nx+1, nmm)
    
            do i = 1, nx
                Hii(i,ik)%m = H00
                H1i(i,ik)%m = H10
                Sii(i,ik)%m = dcmplx(0.0d0, 0.0d0)
                do im = 1, nm
                    Sii(i,ik)%m(im, im) = 1.0d0
                end do
            end do
            H1i(nx+1,ik)%m = H10
        enddo
        deallocate (H00, H10)
    end subroutine devH_build_fromWannierFile

end module deviceHam_mod
