module Output

    implicit none

    private

    integer, parameter :: dp = 8

    public::write_spectrum_summed_over_k

contains

    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_summed_over_k(dataset, i, G, nen, en, nk, length, NB, NS, Lx, coeff, append)
        use matrix_c, only: type_matrix_complex, sizeof
        character(len=*), intent(in) :: dataset
        type(type_matrix_complex), intent(in), dimension(length, nen, nk)::G
        integer, intent(in)::i, nen, length, NB(length), NS(length), nk
        real(dp), intent(in)::Lx(sum(NS)), en(nen), coeff(2)
        logical, intent(in), optional::append
        ! ----
        integer:: ie, j, ib, k, ik, nm(2, length)
        real(dp)::xx
        complex(dp)::tr
        logical append_
        append_ = .false.
        if (present(append)) append_ = append
        if (append_) then
            open (unit=11, file=trim(dataset)//TRIM(STRING(i))//'.dat', status='unknown', position='append')
        else
            open (unit=11, file=trim(dataset)//TRIM(STRING(i))//'.dat', status='unknown')
        end if
        !
        nm = sizeof(G(:, 1, 1))
        do ie = 1, nen
            xx = 0.0d0
            ix = 0
            do j = 1, length
                do k = 1, NS(j)
                    ix = ix + 1
                    xx = xx + Lx(ix)
                    tr = 0.0d0
                    !$omp parallel default(shared) private(ik, ib) reduction(+:tr)
                    !$omp do
                    do ik = 1, nk
                        do ib = 1, nb(j)
                            tr = tr + G(j, ie, ik)%m(ib + (k - 1)*NB(j), ib + (k - 1)*NB(j))
                        end do
                    end do
                    !$omp end do
                    !$omp end parallel
                    tr = tr/dble(nk)
                    write (11, '(4E18.4)') xx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)
                end do
            end do
            write (11, *)
        end do
        close (11)
    end subroutine write_spectrum_summed_over_k

end module Output
