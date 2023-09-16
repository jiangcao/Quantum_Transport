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

    subroutine load_COOmatrix(fname, H, nnz, nm, row, col, use0index, iscomplex)
        character(len=*), intent(in) :: fname !! input text file name
        complex(dp), allocatable, intent(out), dimension(:) :: H
        integer, allocatable, intent(out), dimension(:):: row, col
        integer, intent(out)::nnz
        integer, intent(out)::nm
        logical, intent(in), optional::use0index, iscomplex
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
            if ((present(iscomplex)) .and. iscomplex) then
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

    subroutine devH_build_fromCOOfile(hfname, sfname, Hii, H1i, Sii, ext, contactBlockSize, nx, use0index, iscomplex, threshold)
        use matrix_c, only: type_matrix_complex, malloc
        use graph_partition, only: slice, convert_fromCOO
        character(len=*), intent(in) :: hfname !! input H file name
        character(len=*), intent(in), optional :: sfname !! input S file name
        integer, intent(in)::ext(2)!! number of extension blocks on left/right side
        integer, intent(in)::contactBlockSize(2)!! number of orbitals in the contact block
        integer, intent(out):: nx !! total number of slices
        type(type_matrix_complex), dimension(:), intent(inout), allocatable::Hii, H1i, Sii !! Hamiltonian blocks
        logical, intent(in), optional::use0index, iscomplex
        real(dp), intent(in), optional::threshold
        ! ----
        integer:: nnz, norb, i, newnnz, j, num_slices, nmax
        integer, allocatable, dimension(:, :)::nmii, nm1i
        complex(dp), allocatable, dimension(:)::H, newH !! Hamiltonian matrix value in COO
        complex(dp), allocatable, dimension(:)::S, newS !! overlap matrix value in COO
        integer, allocatable, dimension(:)::row, col !! Hamiltonian matrix index in COO
        integer, allocatable, dimension(:)::newrow, newcol !! Hamiltonian matrix index in COO
        integer, allocatable, dimension(:, :)::Slices !! slicing information , refer to [[graph_partition]]
        integer, allocatable, dimension(:, :)::g !! graph , refer to [[graph_partition]]
        integer, dimension(contactBlockSize(1))::E1 !! edge 1 , refer to [[graph_partition]]
        integer, dimension(contactBlockSize(2))::E2 !! edge 2 , refer to [[graph_partition]]
        !
        if (present(sfname)) call load_COOmatrix(trim(Sfname), S, nnz, norb, row, col, use0index, iscomplex)
        call load_COOmatrix(trim(hfname), H, nnz, norb, row, col, use0index, iscomplex)
        !
        ! optionally remove small matrix elements below threshold
        if (present(threshold)) then
            newnnz = count(abs(H) > threshold)
            allocate (newH(newnnz))
            if (present(sfname)) allocate (newS(newnnz))
            allocate (newrow(newnnz))
            allocate (newcol(newnnz))
            j = 0
            do i = 1, nnz
                if (abs(H(i)) > threshold) then
                    j = j + 1
                    newH(j) = H(i)
                    if (present(sfname)) newS(j) = S(i)
                    newrow(j) = row(i)
                    newcol(j) = col(i)
                end if
            end do
            deallocate (H, row, col)
            if (present(sfname)) deallocate (S)
            call move_alloc(newH, H)
            if (present(sfname)) call move_alloc(newS, S)
            call move_alloc(newrow, row)
            call move_alloc(newcol, col)
            nnz = newnnz
        end if
        ! convert sparse matrix to a graph
        call convert_fromCOO(nnz, row, col, g)
        forall (i=1:contactBlockSize(1)) E1(i) = i
        forall (i=1:contactBlockSize(2)) E2(i) = norb - contactBlockSize(2) + i
        ! slice the device
        nmax = 100
        ! try 1 contact slicing from left
        call slice(g, E1, Slices) ! slicing from left to right
        num_slices = size(Slices, 2)
        if (((Slices(1, num_slices) - 1) < contactBlockSize(2)) .or. &
            ((maxval(Slices(1, :)) - minval(Slices(1, :))) > (2*minval(Slices(1, :))))) then
            ! too small right contact , or very unbalanced slicing
            ! try 1 contact slicing from right
            call slice(g, E2, Slices) ! slicing from right to left
            num_slices = size(Slices, 2)
            if (((Slices(1, 1) - 1) < contactBlockSize(1)) .or. &
                ((maxval(Slices(1, :)) - minval(Slices(1, :))) > (2*minval(Slices(1, :))))) then
                ! too small left contact , or very unbalanced slicing
                call slice(g, E1, E2, NMAX, Slices) ! try 2 contact slicing
                num_slices = size(Slices, 2)
            end if
        end if
        print *, ' Slicing info: ', num_slices, ' slices, with block sizes = '
        print '(12 I8)', slices(1, :) - 1
        ! allocate the blocks : left extension|device|right extension
        nx = ext(1) + num_slices + ext(2)
        allocate (nmii(2, nx))
        allocate (nm1i(2, nx + 1))
        nmii(:, 1:ext(1)) = contactBlockSize(1)
        nmii(:, (nx - ext(2) + 1):nx) = contactBlockSize(2)
        nmii(1, ext(1):(nx - ext(2))) = Slices(1, :) - 1
        nmii(2, ext(1):(nx - ext(2))) = Slices(1, :) - 1
        nm1i(:, 1) = nmii(:, 1)
        nm1i(:, nx + 1) = nmii(:, nx)
        nm1i(1, 2:nx) = nmii(1, 2:nx)
        nm1i(2, 2:nx) = nmii(1, 1:(nx - 1))
        allocate (Hii(nx))
        allocate (Sii(nx))
        allocate (H1i(nx + 1))
        call malloc(Hii, nx, nmii)
        call malloc(Sii, nx, nmii)
        call malloc(H1i, nx + 1, nm1i)
        ! build the blocks

        deallocate (H, col, row, nmii, nm1i)
        if (present(sfname)) deallocate (S)
    end subroutine devH_build_fromCOOfile

    subroutine devH_build_fromWannierFile(fname, Hii, H1i, Sii, nx, nslab, nband, nk, k, length_x, lreorder_axis, axis)
        use matrix_c, only: type_matrix_complex, malloc
        use wannierHam3d, only: w90_load_from_file, w90_free_memory, w90_MAT_DEF, nb, Lx
        character(len=*), intent(in)        :: fname !! input text file name
        logical, intent(in), optional :: lreorder_axis !! whether to reorder axis
        integer, intent(in), optional :: axis(3) !! permutation order
        integer, intent(in)::nx !! number of slabs
        integer, intent(in):: nslab !! number of cells per slab
        integer, intent(out):: nband !! number of bands/orbitals per cell (from wannier90)
        type(type_matrix_complex), dimension(:, :), intent(inout), allocatable::Hii, H1i !! Hamiltonian blocks
        type(type_matrix_complex), dimension(:, :), intent(inout), allocatable::Sii !! overlap matrix blocks
        real(dp), intent(in)::k(2, nk)
        real(dp), intent(out)::length_x
        integer, intent(in)::nk
        ! ----
        complex(dp), allocatable, dimension(:, :)::H00, H10
        real(dp)::kx, ky, kz
        integer::nm, i, im, ik
        integer, dimension(2,nx + 1)::nmm
        open (unit=10, file=trim(fname), status='unknown')
        call w90_load_from_file(10, lreorder_axis, axis)
        close (10)
        nm = nb*nslab
        nband = nb
        length_x = Lx
        nmm = nm
        allocate (H00(nm, nm))
        allocate (H10(nm, nm))
        kx = 0.0d0
        !
        allocate (Hii(nx, nk))
        allocate (H1i(nx + 1, nk))
        allocate (Sii(nx, nk))
        do ik = 1, nk
            ky = k(1, ik)
            kz = k(2, ik)
            call w90_MAT_DEF(H00, H10, kx, ky, kz, nslab)
            !
            call malloc(Hii(:, ik), nx, nmm(:,1:nx))
            call malloc(Sii(:, ik), nx, nmm(:,1:nx))
            call malloc(H1i(:, ik), nx + 1, nmm)
            !
            do i = 1, nx
                Hii(i, ik)%m = H00
                H1i(i, ik)%m = H10
                Sii(i, ik)%m = dcmplx(0.0d0, 0.0d0)
                do im = 1, nm
                    Sii(i, ik)%m(im, im) = 1.0d0
                end do
            end do
            H1i(nx + 1, ik)%m = H10
        end do
        deallocate (H00, H10)
        call w90_free_memory()
    end subroutine devH_build_fromWannierFile

end module deviceHam_mod
