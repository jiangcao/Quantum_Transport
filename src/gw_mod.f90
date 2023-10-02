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
! AUTHOR: Jiang Cao, Alexandros Nikolaos Ziogas 
!
module gw_mod

    implicit none

    private

    public :: gw_ephoton_3d_ijs 

contains

! AUTHOR: Alexandros Nikolaos Ziogas 
subroutine energies_to_ijs_2(buf, tmp0, tmp1, NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  ! complex(8), pointer :: p0(:, :, :, :), p1(:, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  ! Assume (global) g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))

  ! 1a. Interpret as g(1:NM*NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1b. Interpret as g(1:local_Nij, 1:comm_size (Nij), 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1c. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (Nij), 1:comm_size (NE))
  p0(1:local_nij, 1:comm_size, 1:nx, 1:local_ne, 1:nk) => buf
  p1(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 5, 2, 3, 4])
  ! p0(1:nm*nm, 1:nx, 1:local_ne, 1:nk) => buf
  ! p1(1:nx, 1:local_ne, 1:nk, 1:nm*nm) => tmp0
  ! p1 = reshape(p0, shape(p1), order = [4, 1, 2, 3])

  ! 2. Redistribute to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE), 1:comm_size (Nij))
  count = local_nij * nx * local_ne * nk
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (Nij))
  p2(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp1
  p3(1:local_ne, 1:comm_size, 1:nk, 1:nx, 1:local_nij) => buf
  p3 = reshape(p2, shape(p3), order = [5, 4, 1, 3, 2])
  ! p2(1:nx, 1:local_ne, 1:nk, 1:local_nij, 1:comm_size) => tmp1
  ! p3(1:local_nij, 1:nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  ! p3 = reshape(p2, shape(p3), order = [2, 3, 5, 1, 4])

end subroutine energies_to_ijs_2



! AUTHOR: Jiang Cao, Alexandros Nikolaos Ziogas 
subroutine gw_ephoton_3d_ijs(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,&
        Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs, comm_size, comm_rank, local_NE, first_local_energy)
  use fft_mod, only : conv1d => conv1d2, corr1d => corr1d2
  integer ( kind = 4), intent(in) :: comm_size, comm_rank, local_NE, first_local_energy
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
  ! complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  ! complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  ! complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  ! complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  ! complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:)::Sii
  complex(8),allocatable,dimension(:,:,:,:)::Mii,M1i
  real(8)::tr(local_NE,nky*nkz),tre(local_NE,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io,j, global_ij, row, col, ij
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk,ne
  complex(8)::dE, B(nm,nm)
  character ( len = 20 ) :: filename
  character ( len = 8 ) :: fmt
  character ( len = 4 ) :: rank_str
  character ( len = 8 ) :: cmethod
  logical append

  real(8) :: local_energies(local_NE)

  integer ( kind = 4 ) :: ierr, local_Nij, first_local_ij, local_NX, first_local_block
  real(8) :: local_sum_tr, global_sum_tr, local_sum_tre, global_sum_tre
  integer ( kind = 8 ) :: num_g, num_p, num_err, num, count

  integer reqs(8), stats(8)

  complex(8), pointer :: g_r_buf(:), g_r_by_energies(:, :, :, :, :), g_r_by_blocks(:, :, :, :)
  complex(8), pointer :: g_greater_extended_buf(:), g_greater_buf(:), g_greater_by_energies(:, :, :, :, :), g_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: g_lesser_extended_buf(:), g_lesser_buf(:), g_lesser_by_energies(:, :, :, :, :), g_lesser_by_blocks(:, :, :, :)

  complex(8), pointer :: sigma_r_gw_buf(:), sigma_r_gw_by_energies(:, :, :, :, :), sigma_r_gw_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_greater_gw_buf(:), sigma_greater_gw_by_energies(:, :, :, :, :), sigma_greater_gw_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_lesser_gw_buf(:), sigma_lesser_gw_by_energies(:, :, :, :, :), sigma_lesser_gw_by_blocks(:, :, :, :)

  complex(8), pointer :: sigma_r_new_buf(:), sigma_r_new_by_energies(:, :, :, :, :), sigma_r_new_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_greater_new_buf(:), sigma_greater_new_by_energies(:, :, :, :, :), sigma_greater_new_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_lesser_new_buf(:), sigma_lesser_new_by_energies(:, :, :, :, :), sigma_lesser_new_by_blocks(:, :, :, :)

  complex(8), pointer :: P_retarded_buf(:), P_retarded_by_energies(:, :, :, :, :), P_retarded_by_blocks(:, :, :, :)
  complex(8), pointer :: P_greater_buf(:), P_greater_by_energies(:, :, :, :, :), P_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: P_lesser_buf(:), P_lesser_by_energies(:, :, :, :, :), P_lesser_by_blocks(:, :, :, :)

  complex(8), pointer :: W_retarded_buf(:), W_retarded_by_energies(:, :, :, :, :), W_retarded_by_blocks(:, :, :, :)
  complex(8), pointer :: W_greater_buf(:), W_greater_by_energies(:, :, :, :, :), W_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: W_lesser_buf(:), W_lesser_by_energies(:, :, :, :, :), W_lesser_by_blocks(:, :, :, :)

  real(8), allocatable, dimension(:, :, :, :, :) :: cur, jdens_local
  real(8), allocatable, dimension(:, :, :, :) :: tot_cur,tot_ecur,tot_cur_local,tot_ecur_local

  complex(8), pointer :: tmp0(:), tmp1(:), g_lesser_extended(:, :, :, :, :), g_greater_extended(:, :, :, :, :)

  complex(8), pointer :: g_greater_t_buf(:), g_greater_t_by_energies(:, :, :, :, :), g_greater_t_by_blocks(:, :, :, :)
  complex(8), pointer :: g_lesser_t_buf(:), g_lesser_t_by_energies(:, :, :, :, :), g_lesser_t_by_blocks(:, :, :, :)

  complex(8), pointer :: g_lesser_photon_left_send(:), g_lesser_photon_left_recv(:)
  complex(8), pointer :: g_lesser_photon_right_send(:), g_lesser_photon_right_recv(:)
  complex(8), pointer :: g_greater_photon_left_send(:), g_greater_photon_left_recv(:)
  complex(8), pointer :: g_greater_photon_right_send(:), g_greater_photon_right_recv(:)

  complex(8), pointer :: g_lesser_photon_buf(:), g_greater_photon_buf(:)
  complex(8), pointer :: g_lesser_photon(:, :, :, :), g_greater_photon(:, :, :, :)

  real(8) :: start, finish, it_start

  real(8), allocatable :: extended_local_energies(:)
  
  complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:,:),Ispec_ik(:,:,:,:),Itot_ik(:,:,:) ! collision integral variables
  
  complex(8),allocatable,dimension(:,:,:)::Pi_retarded_ik,Pi_lesser_ik,Pi_greater_ik,Pi_retarded ! photon Pi self-energies

  ! For debugging/validation
  ! complex(8), pointer :: g_r_buf2(:), g_r_by_energies2(:, :, :, :, :), g_r_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: g_greater_buf2(:), g_greater_by_energies2(:, :, :, :, :), g_greater_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: g_lesser_buf2(:), g_lesser_by_energies2(:, :, :, :, :), g_lesser_by_blocks2(:, :, :, :, :)

  ! complex(8), pointer :: P_retarded_buf2(:), P_retarded_by_energies2(:, :, :, :, :), P_retarded_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: P_greater_buf2(:), P_greater_by_energies2(:, :, :, :, :), P_greater_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: P_lesser_buf2(:), P_lesser_by_energies2(:, :, :, :, :), P_lesser_by_blocks2(:, :, :, :, :)  

  local_Nij = (nm * nm) / comm_size
  first_local_ij = local_Nij * comm_rank + 1

  local_NX = nx / comm_size
  first_local_block = local_NX * comm_rank + 1

  ! num_g = nm * nm * local_NE * nkz * nky * local_NX
  ! num_p = nm * nm * local_NE * nkz * nky * local_NX

  fmt = '(I4.4)'
  write ( rank_str, fmt ) comm_rank
  append = (comm_rank /= 0)

  if (comm_rank == 0) then
    print *,'======== green_rgf_solve_gw_ephoton_3D ========'
  endif

  do ie = 1, local_NE
    local_energies(ie) = en(ie + first_local_energy - 1)
  enddo

  nk = nky * nkz
  ne = nen
  nopphot=floor(hw / (En(2)-En(1)))

  ! real(8) :: extended_local_energies(local_NE + 2 * nopphot)
  allocate(extended_local_energies(local_NE + 2 * nopphot))

  do ie = 1, local_NE
    local_energies(ie) = en(ie + first_local_energy - 1)
  enddo
  extended_local_energies = 0.0d0
  extended_local_energies(nopphot + 1:nopphot + local_NE) = local_energies(:)
  if (comm_rank - 1 >= 0) then
    do ie = 1, nopphot
      extended_local_energies(ie) = en(ie + first_local_energy - nopphot - 1)
    enddo
  endif
  if (comm_rank +1 < comm_size) then
    ! do ie = nopphot + local_NE + 1, local_NE + 2 * nopphot
    do ie = 1, nopphot
      extended_local_energies(ie + nopphot + local_NE) = en(local_NE + first_local_energy + ie - 1)
    enddo
  endif

  allocate(g_r_buf(nm * nm * local_NE * nk * nx))
  allocate(g_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_buf(nm * nm * nx * local_NE * nk))
  ! allocate(g_greater_extended_buf(nm * nm * nx * (local_NE + 2 * nopphot) * nk))
  ! allocate(g_lesser_extended_buf(nm * nm * nx * (local_NE + 2 * nopphot) * nk))
  ! count = nm * nm * nx * nopphot * nk
  ! g_greater_buf(1:nm * nm * local_NE * nk * nx) => g_greater_extended_buf(count:count + nm * nm * local_NE * nk * nx - 1)
  ! g_lesser_buf(1:nm * nm * local_NE * nk * nx) => g_lesser_extended_buf(count:count + nm * nm * local_NE * nk * nx - 1)

  ! For photon computation/communication
  allocate(g_lesser_photon_buf(nm * nm * nx * (local_NE + 2 * nopphot)))
  allocate(g_greater_photon_buf(nm * nm * nx * (local_NE + 2 * nopphot)))
  g_lesser_photon_buf = dcmplx(0.0d0, 0.0d0)
  g_greater_photon_buf = dcmplx(0.0d0, 0.0d0)
  g_lesser_photon(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot) => g_lesser_photon_buf
  g_greater_photon(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot) => g_greater_photon_buf
  g_lesser_photon_left_send(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * nopphot + 1 : nm * nm * nx * 2 * nopphot)
  g_lesser_photon_left_recv(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(1 : nm * nm * nx * nopphot)
  g_lesser_photon_right_send(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * local_NE + 1 : nm * nm * nx * (local_NE + nopphot))
  g_lesser_photon_right_recv(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * (local_NE + nopphot) + 1 : nm * nm * nx * (local_NE + 2 * nopphot))
  g_greater_photon_left_send(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * nopphot + 1 : nm * nm * nx * 2 * nopphot)
  g_greater_photon_left_recv(1:nm * nm * nx* nopphot) => g_greater_photon_buf(1 : nm * nm * nx * nopphot)
  g_greater_photon_right_send(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * local_NE + 1 : nm * nm * nx * (local_NE + nopphot))
  g_greater_photon_right_recv(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * (local_NE + nopphot) + 1 : nm * nm * nx * (local_NE + 2 * nopphot))

  ! g_lesser_extended_buf = dcmplx(0.0d0,0.0d0)
  ! g_greater_extended_buf = dcmplx(0.0d0,0.0d0)
  ! g_lesser_extended(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot, 1:nk) => g_lesser_extended_buf
  ! g_greater_extended(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot, 1:nk) => g_greater_extended_buf

  g_r_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_r_buf
  g_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_buf
  g_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_buf

  ! g_r_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_r_buf
  ! g_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_greater_buf
  ! g_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_lesser_buf
  g_r_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_r_buf
  g_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_greater_buf
  g_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_lesser_buf

  allocate(sigma_r_gw_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_greater_gw_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_lesser_gw_buf(nm * nm * nx * local_NE * nk))

  sigma_r_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_r_gw_buf
  sigma_greater_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_greater_gw_buf
  sigma_lesser_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_lesser_gw_buf

  sigma_r_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_r_gw_buf
  sigma_greater_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_greater_gw_buf
  sigma_lesser_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_lesser_gw_buf

  allocate(sigma_r_new_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_greater_new_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_lesser_new_buf(nm * nm * nx * local_NE * nk))

  sigma_r_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_r_new_buf
  sigma_greater_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_greater_new_buf
  sigma_lesser_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_lesser_new_buf

  ! sigma_r_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_r_new_buf
  ! sigma_greater_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_greater_new_buf
  ! sigma_lesser_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_lesser_new_buf
  sigma_r_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_r_new_buf
  sigma_greater_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_greater_new_buf
  sigma_lesser_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_lesser_new_buf

  allocate(P_retarded_buf(nm * nm * nx * local_NE * nk))
  allocate(P_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(P_lesser_buf(nm * nm * nx * local_NE * nk))

  P_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_retarded_buf
  P_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_greater_buf
  P_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_lesser_buf

  ! P_retarded_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_retarded_buf
  ! P_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_greater_buf
  ! P_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_lesser_buf
  P_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_retarded_buf
  P_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_greater_buf
  P_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_lesser_buf

  allocate(W_retarded_buf(nm * nm * nx * local_NE * nk))  
  allocate(W_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(W_lesser_buf(nm * nm * nx * local_NE * nk))

  W_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_retarded_buf
  W_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_greater_buf
  W_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_lesser_buf

  ! W_retarded_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_retarded_buf
  ! W_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_greater_buf
  ! W_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_lesser_buf
  W_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_retarded_buf
  W_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_greater_buf
  W_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_lesser_buf

  allocate(cur(nm, nm, nx, local_NE, nk))

  allocate(tmp0(nm * nm * nx * local_NE * nk))
  allocate(tmp1(nm * nm * nx * local_NE * nk))

  allocate(g_greater_t_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_t_buf(nm * nm * nx * local_NE * nk))

  g_greater_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_t_buf
  g_lesser_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_t_buf

  ! g_greater_t_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_greater_t_buf
  ! g_lesser_t_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_lesser_t_buf
  g_greater_t_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_greater_t_buf
  g_lesser_t_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_lesser_t_buf

  sigma_greater_gw_buf = dcmplx(0.0d0,0.0d0)
  sigma_lesser_gw_buf = dcmplx(0.0d0,0.0d0)
  sigma_r_gw_buf = dcmplx(0.0d0,0.0d0)

  allocate(Mii(nm,nm,nx,nk))
  Mii(:,:,:,:)=czero
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo
  if (comm_rank == 0) then
    print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity
  endif
  ! nopphot=floor(hw / (En(2)-En(1)))
  if (comm_rank == 0) then
    print *,'nop photon=',nopphot
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  do iter=0,niter
    if (comm_rank == 0) then
      print *,'+ iter=',iter
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    start = MPI_Wtime()
    it_start = start

    if (comm_rank == 0) then
      print *, 'Computing G ...'
    endif      
        
    allocate(jdens_local(nb,nb,nx*ns,local_NE,nk))
    allocate(tot_cur_local(nb,nb,nx*ns,nk))
    allocate(tot_ecur_local(nb,nb,nx*ns,nk))
    allocate(tot_cur(nb,nb,nx*ns,nk))
    allocate(tot_ecur(nb,nb,nx*ns,nk))
        
    do ik=1,nk
      if (comm_rank == 0) then
        print *, ' ik=', ik,'/',nk
      endif
      !$omp parallel default(shared) private(ie)
      !$omp do 
      do ie=1, local_NE     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call green_RGF_RS( &
          TEMP, nm, nx, local_energies(ie), mu, Hii(:,:,:,ik), H1i(:,:,:,ik), &
          sigma_lesser_gw_by_energies(:,:,:,ie,ik),sigma_greater_gw_by_energies(:,:,:,ie,ik), sigma_r_gw_by_energies(:,:,:,ie,ik), &
          g_lesser_by_energies(:,:,:,ie,ik), g_greater_by_energies(:,:,:,ie,ik), g_r_by_energies(:,:,:,ie,ik), &
          tre(ie,ik), tr(ie,ik), cur(:,:,:,ie,ik) ) 
      enddo
      !$omp end do 
      !$omp end parallel
      call calc_block_current(Hii(:,:,:,ik),g_lesser_by_energies(:,:,:,:,ik),cur(:,:,:,:,ik),local_NE,local_energies,spindeg,nb,ns,nm,nx,tot_cur_local(:,:,:,ik),tot_ecur_local(:,:,:,ik),jdens_local(:,:,:,:,ik))    
    enddo
    

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("G computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing G ...'
    endif
    start = finish

    do i = 0, comm_size - 1
      if (i == comm_rank) then
        filename = 'gw_ldos'
        call write_spectrum_summed_over_k(filename,iter,g_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/), append)
!        filename = 'gw_ndos'
!        call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_pdos'   
!        call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        filename = 'gw_Jdens'
        call write_current_spectrum_summed_over_kz(filename,iter,jdens_local,local_NE,local_energies,nx*NS,NB,Lx,nk, append)        
        filename = 'gw_trL'
        call write_transmission_spectrum_k(filename,iter,tr,local_NE,local_energies,nk, append)
        filename = 'gw_trR'
        call write_transmission_spectrum_k(filename,iter,tre,local_NE,local_energies,nk, append)
!        call write_dos_summed_over_k('gw_dos',iter,G_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx, append)
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

    ! Reduce tr and tre
    local_sum_tr = sum(tr)
    local_sum_tre = sum(tre)
    call MPI_Reduce(local_sum_tr, global_sum_tr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(local_sum_tre, global_sum_tre, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Reduce current and energy-current
    call MPI_Reduce(tot_cur_local, tot_cur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(tot_ecur_local, tot_ecur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (comm_rank == 0) then
      open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
      write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
      close(101)
      call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
      call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)      
    endif
    
    deallocate(tot_cur,tot_ecur,jdens_local,tot_cur_local,tot_ecur_local)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("G storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing P ...'
    endif
    start = finish

    g_r_buf = dcmplx( 0.0d0*dble(g_r_buf), aimag(g_r_buf))
    g_lesser_buf = dcmplx( 0.0d0*dble(g_lesser_buf), aimag(g_lesser_buf))
    g_greater_buf = dcmplx( 0.0d0*dble(g_greater_buf), aimag(g_greater_buf))

    g_greater_t_by_energies = reshape(g_greater_by_energies, shape(g_greater_t_by_energies), order=[2, 1, 3, 4, 5])
    g_lesser_t_by_energies = reshape(g_lesser_by_energies, shape(g_lesser_t_by_energies), order=[2, 1, 3, 4, 5])

    ! Redistribute G
    call energies_to_ijs_2(g_r_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_lesser_t_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_greater_t_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    P_lesser_buf = dcmplx(0.0d0,0.0d0)
    P_greater_buf = dcmplx(0.0d0,0.0d0)    
    P_retarded_buf = dcmplx(0.0d0,0.0d0)    

    nopmax=nen/2-1
    
    cmethod='fft'
    call MKL_SET_NUM_THREADS(4)

    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)             
    do i = 1, local_Nij
      global_ij = i + first_local_ij - 1
      col = (global_ij - 1) / nm + 1  ! convert to 0-based indexing, divide, and convert back to 1-based indexing
      row = mod(global_ij - 1, nm) + 1  ! convert to 0-based indexing, mod, and convert back to 1-based indexing
      l = max(row - ndiag, 1)
      h = min(nm, row + ndiag)
      if (col >= l .and. col <= h) then
      !$omp parallel default(shared) private(ix,iop,nop,ie,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
      !$omp do
      do ix = 1, nx
        do iqy=1,nky        
          do iqz=1,nkz
            iq=iqz + (iqy-1)*nkz
            do iky=1,nky
              do ikz=1,nkz              
                ik=ikz + (iky-1)*nkz
                ikzd=ikz-iqz + nkz/2
                ikyd=iky-iqy + nky/2
                if (ikzd<1)   ikzd=ikzd+nkz
                if (ikzd>nkz) ikzd=ikzd-nkz
                if (ikyd<1)   ikyd=ikyd+nky
                if (ikyd>nky) ikyd=ikyd-nky                
                if (nky==1)   ikyd=1
                if (nkz==1)   ikzd=1
                ikd=ikzd + (ikyd-1)*nkz
                
                P_lesser_by_blocks(:,iq,ix,i) = P_lesser_by_blocks(:,iq,ix,i) + &
                                                corr1d(nen,G_lesser_by_blocks(:,ik,ix,i),G_greater_t_by_blocks(:,ikd,ix,i),method=cmethod)
                
                P_greater_by_blocks(:,iq,ix,i) = P_greater_by_blocks(:,iq,ix,i) + & 
                                                corr1d(nen,G_greater_by_blocks(:,ik,ix,i),G_lesser_t_by_blocks(:,ikd,ix,i),method=cmethod)
                
                P_retarded_by_blocks(:,iq,ix,i) = P_retarded_by_blocks(:,iq,ix,i) + & 
                                                  corr1d(nen,G_lesser_by_blocks(:,ik,ix,i),conjg(G_r_by_blocks(:,ikd,ix,i)),method=cmethod) + &
                                                  corr1d(nen,G_r_by_blocks(:,ik,ix,i),G_lesser_t_by_blocks(:,ikd,ix,i),method=cmethod)
                
                
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do 
      !$omp end parallel
    endif     
    enddo          

    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser_buf = P_lesser_buf * dE
    P_greater_buf = P_greater_buf * dE
    P_retarded_buf = P_retarded_buf * dE

    ! Redistribute P
    call ijs_to_energies_2(P_retarded_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(P_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(P_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    ! Zero-out off diagonals
    !$omp parallel default(shared) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
    !$omp do
    do ik = 1, nk
      do ie = 1, local_NE
        do ix = 1, nx
          do i=1,nm      
            l = max(i-ndiag,1)
            h = min(nm,i+ndiag)
            if (l > 1) then
              P_lesser_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_greater_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_retarded_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
            endif
            if (h < nm) then
              P_lesser_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_greater_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_retarded_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0) 
            endif    
          enddo      
        enddo
      enddo
    enddo    
    !$omp end do 
    !$omp end parallel

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("P computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing P ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_PR'
!        call write_spectrum_summed_over_k(filename,iter,P_retarded_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_PL'
!        call write_spectrum_summed_over_k(filename,iter,P_lesser_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_PG'
!        call write_spectrum_summed_over_k(filename,iter,P_greater_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("P storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing W ...'
    endif
    start = finish

    do iq=1,nky*nkz

      if (comm_rank == 0) then
        print *, ' iq=', iq,'/',nk
      endif

      !$omp parallel default(shared) private(nop, ie)
      !$omp do 
      do nop = max(-nopmax + nen/2, first_local_energy), min(nopmax + nen/2, first_local_energy + local_NE - 1)
        ie = nop - first_local_energy + 1
!        call green_calc_w_full( &
!          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
!          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
!          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))

        call green_rgf_calc_w( &
          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))  
      enddo
      !$omp end do 
      !$omp end parallel
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("W computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing W ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_WR'
!        call write_spectrum_summed_over_k(filename,iter,W_retarded_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_WL'
!        call write_spectrum_summed_over_k(filename,iter,W_lesser_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_WG'
!        call write_spectrum_summed_over_k(filename,iter,W_greater_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("W storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing SigGW ...'
    endif
    start = finish

    call energies_to_ijs_2(W_retarded_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(W_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(W_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    Sigma_greater_new_buf = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new_buf = dcmplx(0.0d0,0.0d0)
    Sigma_r_new_buf = dcmplx(0.0d0,0.0d0)
  
    cmethod='fft'
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)          
    do i = 1, local_Nij
      global_ij = i + first_local_ij - 1
      col = (global_ij - 1) / nm + 1  ! convert to 0-based indexing, divide, and convert back to 1-based indexing
      row = mod(global_ij - 1, nm) + 1  ! convert to 0-based indexing, mod, and convert back to 1-based indexing
      l = max(row - ndiag, 1)
      h = min(nm, row + ndiag)
      if (col >= l .and. col <= h) then
        !$omp parallel default(shared) private(ix,iop,nop,ie,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
        !$omp do  
        do ix=1,nx     
          do iky=1,nky
            do ikz=1,nkz   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  
                  sigma_lesser_new_by_blocks(:, ik, ix, i) = sigma_lesser_new_by_blocks(:, ik, ix, i) + &
                                                      conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i),W_lesser_by_blocks(:, iq, ix, i),method=cmethod)
                  
                  Sigma_greater_new_by_blocks(:, ik, ix, i) = Sigma_greater_new_by_blocks(:, ik, ix, i) + &
                                                      conv1d(nen,G_greater_by_blocks(:, ikd, ix, i),W_greater_by_blocks(:, iq, ix, i),method=cmethod)
                  
                  Sigma_r_new_by_blocks(:, ik, ix, i) = Sigma_r_new_by_blocks(:, ik, ix, i) + & 
                                                      conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                      conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_lesser_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                      conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i) ,method=cmethod)
                  
                  
                enddo
              enddo            
            enddo
          enddo
        enddo
        !$omp end do
        !$omp end parallel
      endif      
    enddo

    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new_buf = Sigma_lesser_new_buf  * dE
    Sigma_greater_new_buf = Sigma_greater_new_buf * dE
    Sigma_r_new_buf = Sigma_r_new_buf * dE
    Sigma_r_new_buf = dcmplx( dble(Sigma_r_new_buf), aimag(Sigma_greater_new_buf-Sigma_lesser_new_buf)/2.0d0 )

    call ijs_to_energies_2(sigma_r_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(sigma_lesser_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(sigma_greater_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    ! Zero-out off diagonals
    !$omp parallel default(shared) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
    !$omp do
    do ik = 1, nk
      do ie = 1, local_NE
        do ix = 1, nx
          do i=1,nm      
            l = max(i-ndiag,1)
            h = min(nm,i+ndiag)
            if (l > 1) then
              Sigma_lesser_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_greater_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_r_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
            endif
            if (h < nm) then
              Sigma_lesser_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_greater_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_r_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0) 
            endif    
          enddo      
        enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel   
 
    ! symmetrize the selfenergies
    ! do ie=1,nen
    do ie = 1, local_NE
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new_by_energies(:,:,ix,ie,ik))
          Sigma_r_new_by_energies(:,:,ix,ie,ik) = (Sigma_r_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new_by_energies(:,:,ix,ie,ik))
          Sigma_lesser_new_by_energies(:,:,ix,ie,ik) = (Sigma_lesser_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new_by_energies(:,:,ix,ie,ik))
          Sigma_greater_new_by_energies(:,:,ix,ie,ik) = (Sigma_greater_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo

    ! mixing with the previous one

    Sigma_r_gw_buf = Sigma_r_gw_buf + alpha_mix * (Sigma_r_new_buf - Sigma_r_gw_buf)
    Sigma_lesser_gw_buf  = Sigma_lesser_gw_buf + alpha_mix * (Sigma_lesser_new_buf - Sigma_lesser_gw_buf)
    Sigma_greater_gw_buf = Sigma_greater_gw_buf + alpha_mix * (Sigma_greater_new_buf -Sigma_greater_gw_buf)    
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("SigGW computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing SigGW ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_SigR'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_r_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_SigL'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_lesser_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_SigG'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_greater_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("SigGW storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing GW scattering rates ...'
    endif
    start = finish
      
    ! We need to convert G back to energies     
    call ijs_to_energies_2(g_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(g_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    !!!!! calculate collision integral
    
    allocate(Ispec(nm,nm,nx,local_NE))
    allocate(Itot(nm,nm,nx))
    allocate(Ispec_ik(nm,nm,nx,local_NE))
    allocate(Itot_ik(nm,nm,nx))
    !
    Ispec=czero
    Itot=czero
    do ik=1,nk
      call calc_block_collision(sigma_lesser_gw_by_energies(:,:,:,:,ik),sigma_greater_gw_by_energies(:,:,:,:,ik),G_lesser_by_energies(:,:,:,:,ik),G_greater_by_energies(:,:,:,:,ik),local_NE,local_energies,spindeg,nm,nx,Itot_ik,Ispec_ik)
      Ispec=Ispec+Ispec_ik
      Itot=Itot+Itot_ik
    enddo
    
    do i = 0, comm_size - 1
      if (i == comm_rank) then
        filename = 'gw_Scat'
        call write_spectrum(filename,iter,Ispec,local_NE,local_energies,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/),append)    
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo    
        
    finish = MPI_Wtime()
    if (comm_rank == 0) then
      print '("Scattering computation time = ", F0.3 ," seconds.")', finish-start      
    endif
    
    if (iter>=(niter-5)) then

      if (comm_rank == 0) then
        print *, 'Computing SigEPhoton ...'
      endif
      start = finish

      do ik=1, nk

        write (rank_str, fmt) ik

        if (comm_rank == 0) then
          print *, ' ik=', ik,'/',nk
        endif

        ! do i = 0, comm_size - 1
        !   if (i == comm_rank) then
        !     filename = 'gw_GL_ik' // rank_str // '_'
        !     call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies(:,:,:,:,ik:ik),local_NE,local_energies,1,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
        !     filename = 'gw_GG_ik' // rank_str // '_'  
        !     call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies(:,:,:,:,ik:ik),local_NE,local_energies,1,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        !   endif
        !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
        ! enddo

        ! Copy local_NE energies to the buffers
        g_lesser_photon_buf = dcmplx(0.0d0,0.0d0)
        g_greater_photon_buf = dcmplx(0.0d0,0.0d0)
        g_lesser_photon(:, :, :, nopphot + 1 : nopphot + local_NE) = g_lesser_by_energies(:, :, :, :, ik)
        g_greater_photon(:, :, :, nopphot + 1 : nopphot + local_NE) = g_greater_by_energies(:, :, :, :, ik)

        ! Communicate up to 2*nopphot energies to the left and right
        if (comm_rank - 1 >= 0) then
          call MPI_Isend(g_lesser_photon_left_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 0, MPI_COMM_WORLD, reqs(1), ierr)
          call MPI_Irecv(g_lesser_photon_left_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 1, MPI_COMM_WORLD, reqs(2), ierr)
          call MPI_Isend(g_greater_photon_left_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 2, MPI_COMM_WORLD, reqs(3), ierr)
          call MPI_Irecv(g_greater_photon_left_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 3, MPI_COMM_WORLD, reqs(4), ierr)
        endif
        if (comm_rank + 1 < comm_size) then
          call MPI_Isend(g_lesser_photon_right_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 1, MPI_COMM_WORLD, reqs(5), ierr)
          call MPI_Irecv(g_lesser_photon_right_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 0, MPI_COMM_WORLD, reqs(6), ierr)
          call MPI_Isend(g_greater_photon_right_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 3, MPI_COMM_WORLD, reqs(7), ierr)
          call MPI_Irecv(g_greater_photon_right_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 2, MPI_COMM_WORLD, reqs(8), ierr)
        endif
  
        if (comm_rank - 1 >= 0) then
          call MPI_Waitall(4, reqs(1:4), MPI_STATUSES_IGNORE, ierr)
        endif
        if (comm_rank +1 < comm_size) then
          call MPI_Waitall(4, reqs(5:8), MPI_STATUSES_IGNORE, ierr)
        endif

        ! do i = 0, comm_size - 1
        !   if (i == comm_rank) then
        !     filename = 'gw_GL_ext_ik' // rank_str // '_'
        !     call write_spectrum_summed_over_k(filename,iter,g_lesser_photon,local_NE + 2 * nopphot,extended_local_energies,1,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
        !     filename = 'gw_GG_ext_ik' // rank_str // '_'   
        !     call write_spectrum_summed_over_k(filename,iter,g_greater_photon,local_NE + 2 * nopphot,extended_local_energies,1,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        !   endif
        !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
        ! enddo

        call calc_sigma_ephoton_monochromatic_ext( &
          nm, NX, local_NE, En, nopphot, Mii(:,:,:,ik), &
          g_lesser_photon(:,:,:,:), g_greater_photon(:,:,:,:), &
          Sigma_r_new_by_energies(:,:,:,:,ik), Sigma_lesser_new_by_energies(:,:,:,:,ik), Sigma_greater_new_by_energies(:,:,:,:,ik), &
          pre_fact, intensity, hw)
      enddo

      ! call ijs_to_energies(sigma_r_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies(sigma_lesser_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies(sigma_greater_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

      Sigma_r_gw_buf       = Sigma_r_gw_buf + Sigma_r_new_buf 
      Sigma_lesser_gw_buf  = Sigma_lesser_gw_buf + Sigma_lesser_new_buf 
      Sigma_greater_gw_buf = Sigma_greater_gw_buf + Sigma_greater_new_buf

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("SigEPhoton computation time = ", F0.3 ," seconds.")', finish-start
        print *, 'Storing SigEphoton ...'
      endif
      start = finish

!      do i = 0, comm_size - 1
!        if (i == comm_rank) then
!          filename = 'eph_SigR_'
!          call write_spectrum_summed_over_k(filename,iter,sigma_r_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!          filename = 'eph_SigL'
!          call write_spectrum_summed_over_k(filename,iter,sigma_lesser_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!          filename = 'eph_SigG'
!          call write_spectrum_summed_over_k(filename,iter,sigma_greater_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        endif
!        call MPI_Barrier(MPI_COMM_WORLD, ierr)
!      enddo

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("SigEPhoton storage time = ", F0.3 ," seconds.")', finish-start
        print *, 'Computing Ephoton scattering rates ...'
      endif
      
      !!!! calculate e-photon scattering rates
            
      start = finish
      
      Ispec=czero
      Itot=czero
      do ik=1,nk
        call calc_block_collision(sigma_lesser_new_by_energies(:,:,:,:,ik),sigma_greater_new_by_energies(:,:,:,:,ik),&
            G_lesser_by_energies(:,:,:,:,ik),G_greater_by_energies(:,:,:,:,ik),local_NE,local_energies,spindeg,nm,nx,Itot_ik,Ispec_ik)
        Ispec=Ispec+Ispec_ik
        Itot=Itot+Itot_ik
      enddo
      do i = 0, comm_size - 1
        if (i == comm_rank) then
          filename = 'eph_Scat'
          call write_spectrum(filename,iter,Ispec,local_NE,local_energies,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/),append)    
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
      enddo                 
      
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("EPhoton scattering time = ", F0.3 ," seconds.")', finish-start        
      endif 

    endif
    
    deallocate(Ispec,Itot,Ispec_ik,Itot_ik) 
    
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      sigma_r_gw_by_energies(:,:,ix,:,:)=Sigma_r_gw_by_energies(:,:,3,:,:)
      sigma_lesser_gw_by_energies(:,:,ix,:,:)=Sigma_lesser_gw_by_energies(:,:,3,:,:)
      sigma_greater_gw_by_energies(:,:,ix,:,:)=Sigma_greater_gw_by_energies(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_r_gw_by_energies(:,:,nx-2,:,:)
      Sigma_lesser_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_lesser_gw_by_energies(:,:,nx-2,:,:)
      Sigma_greater_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_greater_gw_by_energies(:,:,nx-2,:,:)
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("Iteration ", I3.3, " time = ", F0.3 ," seconds.")', iter, finish-it_start
      print *, ''
    endif

  enddo
  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  start = MPI_Wtime()
  it_start = start

  if (comm_rank == 0) then
    print *, 'Computing G ...'
  endif      
      
  allocate(jdens_local(nb,nb,nx*ns,local_NE,nk))
  allocate(tot_cur_local(nb,nb,nx*ns,nk))
  allocate(tot_ecur_local(nb,nb,nx*ns,nk))
  allocate(tot_cur(nb,nb,nx*ns,nk))
  allocate(tot_ecur(nb,nb,nx*ns,nk))
      
  do ik=1,nk
    if (comm_rank == 0) then
      print *, ' ik=', ik,'/',nk
    endif
    !$omp parallel default(shared) private(ie)
    !$omp do 
    do ie=1, local_NE     
      !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_RGF_RS( &
        TEMP, nm, nx, local_energies(ie), mu, Hii(:,:,:,ik), H1i(:,:,:,ik), &
        sigma_lesser_gw_by_energies(:,:,:,ie,ik),sigma_greater_gw_by_energies(:,:,:,ie,ik), sigma_r_gw_by_energies(:,:,:,ie,ik), &
        g_lesser_by_energies(:,:,:,ie,ik), g_greater_by_energies(:,:,:,ie,ik), g_r_by_energies(:,:,:,ie,ik), &
        tre(ie,ik), tr(ie,ik), cur(:,:,:,ie,ik) ) 
    enddo
    !$omp end do 
    !$omp end parallel
    call calc_block_current(Hii(:,:,:,ik),g_lesser_by_energies(:,:,:,:,ik),cur(:,:,:,:,ik),local_NE,local_energies,spindeg,nb,ns,nm,nx,tot_cur_local(:,:,:,ik),tot_ecur_local(:,:,:,ik),jdens_local(:,:,:,:,ik))    
  enddo
  

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  finish = MPI_Wtime()

  if (comm_rank == 0) then
    print '("G computation time = ", F0.3 ," seconds.")', finish-start
    print *, 'Storing G ...'
  endif
  start = finish

  do i = 0, comm_size - 1
    if (i == comm_rank) then
      filename = 'gw_ldos'
      call write_spectrum_summed_over_k(filename,iter,g_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/), append)
      filename = 'gw_ndos'
      call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
      filename = 'gw_pdos'   
      call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
      filename = 'gw_Jdens'
      call write_current_spectrum_summed_over_kz(filename,iter,jdens_local,local_NE,local_energies,nx*NS,NB,Lx,nk, append)        
      filename = 'gw_trL'
      call write_transmission_spectrum_k(filename,iter,tr*spindeg,local_NE,local_energies,nk, append)
      filename = 'gw_trR'
      call write_transmission_spectrum_k(filename,iter,tre*spindeg,local_NE,local_energies,nk, append)
      call write_dos_summed_over_k('gw_dos',iter,G_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx, append)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo

  ! Reduce tr and tre
  local_sum_tr = sum(tr)
  local_sum_tre = sum(tre)
  call MPI_Reduce(local_sum_tr, global_sum_tr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(local_sum_tre, global_sum_tre, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  ! Reduce current and energy-current
  call MPI_Reduce(tot_cur_local, tot_cur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(tot_ecur_local, tot_ecur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (comm_rank == 0) then
    open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
    close(101)
    call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
    call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)      
  endif
  
  deallocate(tot_cur,tot_ecur,jdens_local,tot_cur_local,tot_ecur_local)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  finish = MPI_Wtime()

  if (comm_rank == 0) then
    print '("G storage time = ", F0.3 ," seconds.")', finish-start    
  endif

  deallocate( cur, tmp0, tmp1)
  ! deallocate(g_r_buf, g_lesser_extended_buf, g_greater_extended_buf)
  deallocate(g_r_buf, g_lesser_buf, g_greater_buf)
  deallocate(g_lesser_photon_buf, g_greater_photon_buf)
  deallocate(sigma_lesser_gw_buf, sigma_greater_gw_buf, sigma_r_gw_buf)
  deallocate(sigma_lesser_new_buf, sigma_greater_new_buf, sigma_r_new_buf)
  deallocate(P_retarded_buf, P_lesser_buf, P_greater_buf)
  deallocate(W_retarded_buf, W_lesser_buf, W_greater_buf)
  deallocate(Mii)

  deallocate(extended_local_energies)
end subroutine gw_ephoton_3d_ijs



! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic_ext(nm,length,nen,En,nop,Mii,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,pre_fact,intensity,hw)
  integer,intent(in)::nm,length,nen,nop
  real(8),intent(in)::en(nen),pre_fact,intensity,hw
  complex(8),intent(in),dimension(nm,nm,length)::Mii ! e-photon interaction matrix blocks
  complex(8),intent(in),dimension(nm,nm,length,nen+2*nop)::G_lesser,G_greater
  complex(8),intent(inout),dimension(nm,nm,length,nen)::Sig_retarded,Sig_lesser,Sig_greater
  !---------
  integer::ie,ix,ie2
  complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
    Sig_lesser=czero
    Sig_greater=czero
    Sig_retarded=czero
    ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
    !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
    !$omp parallel default(none) private(ie,A,B,ix,ie2) shared(length,nop,nen,nm,G_lesser,G_greater,Sig_lesser,Sig_greater,Mii)
    allocate(B(nm,nm))
    allocate(A(nm,nm))  
    !$omp do
    do ie=1,nen
      ie2 = ie+nop
      do ix=1,length
        ! Sig^<(E)
        A = czero
        ! if (ie-nop>=1) A =A+ G_lesser(:,:,ix,ie2-nop)
        ! if (ie+nop<=nen) A =A+ G_lesser(:,:,ix,ie2+nop)
        A =A+ G_lesser(:,:,ix,ie2-nop)
        A =A+ G_lesser(:,:,ix,ie2+nop)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_lesser(:,:,ix,ie) = A 
        ! Sig^>(E)
        A = czero
        ! if (ie-nop>=1) A =A+ G_greater(:,:,ix,ie2-nop)
        ! if (ie+nop<=nen) A =A+ G_greater(:,:,ix,ie2+nop)
        A =A+ G_greater(:,:,ix,ie2-nop)
        A =A+ G_greater(:,:,ix,ie2+nop)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_greater(:,:,ix,ie) = A
      enddo
    enddo  
    !$omp end do
    deallocate(A,B)
    !$omp end parallel
    Sig_greater=Sig_greater*pre_fact*intensity/hw**2
    Sig_lesser=Sig_lesser*pre_fact*intensity/hw**2
    Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
  end subroutine calc_sigma_ephoton_monochromatic_ext





!!!!!! need to merge with rgf_mod

! RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_RS(TEMP,nm,nx,E,mu,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,GRi1,GLi1,GGi1)     
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: ldos,ndens,pdens ! diag blocks of GFs    
    real(8),intent(inout),dimension(nm,nm,Nx) :: cur ! current density
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    ! ------- 
    COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    COMPLEX(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr,tmp1
    COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
    REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
    mul=mu(1)
    mur=mu(2)
    templ=temp(1)
    tempr=temp(2)
    !
    z=dcmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    S00=0.0d0
    do i=1,nm
        S00(i,i)=1.0d0
    enddo
    ! self energy on the left contact        
    H00(:,:)=Hii(:,:,1)+sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    tmp1=maxval(abs(transpose(G00) - G00)) / maxval(abs(G00))
    if (tmp1 > 2e-2) then
        print *, 'E=',E, '|G00^t-G00|=', tmp1
    endif
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)      
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+sigma_lesser_ph(:,:,1)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    Do l=2,nx-1                
        H00(:,:)=Hii(:,:,l)+sigma_r_ph(:,:,l)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)             
        C(:,:)=C(:,:)+sigma_lesser_ph(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact        
    H00(:,:)=Hii(:,:,nx)+sigma_r_ph(:,:,nx)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)!cmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01    
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+sigma_lesser_ph(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm        
      tim=tim-dble(B(i,i)-C(i,i))        
    enddo
    tr=tim
    ! transmission
    !-------------------------
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
        !-------------------------
!        A(:,:)=Gp(:,:)
!        call zgemm('c','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
!        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
!        !
!        Gp(:,:)= Glp(:,:,l) + C(:,:)
!        A(:,:)=Glp(:,:,l)
!        call zgemm('c','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
!        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
!        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
        pdens(:,:,l)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
      tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
    deallocate(Gl)
    deallocate(Gln)
    deallocate(Glp)
end subroutine green_RGF_RS


subroutine invert(A,nn)  
  integer :: info,lda,lwork,nn      
  integer, dimension(:), allocatable :: ipiv
  complex(8), dimension(nn,nn),intent(inout) :: A
  complex(8), dimension(:), allocatable :: work
  allocate(work(nn*nn))
  allocate(ipiv(nn))
  call zgetrf(nn,nn,A,nn,ipiv,info)
  if (info.ne.0) then
    print*,'SEVERE warning: zgetrf failed, info=',info
    A=czero
  else
    call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
    if (info.ne.0) then
      print*,'SEVERE warning: zgetri failed, info=',info
      A=czero
    endif
  endif
  deallocate(work)
  deallocate(ipiv)
end subroutine invert


Function ferm(a)
    Real (8) a,ferm
    ferm=1.0d0/(1.0d0+Exp(a))
End Function ferm


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)    
    complex(8), parameter :: alpha = dcmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = dcmplx(0.0d0,0.0d0)
    integer i,j,k,nm,nmax
    COMPLEX(8) :: z
    real(8) :: E,error
    REAL(8) :: TOL=1.0D-100  ! [eV]
    COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
    COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)
    COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:)
    COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
    !COMPLEX(8), ALLOCATABLE :: WORK(:)
    !COMPLEX(8), EXTERNAL :: ZLANGE
    Allocate( H_BB(nm,nm) )
    Allocate( H_SS(nm,nm) )
    Allocate( H_01(nm,nm) )
    Allocate( H_10(nm,nm) )
    Allocate( Id(nm,nm) )
    Allocate( A(nm,nm) )
    Allocate( B(nm,nm) )
    Allocate( C(nm,nm) )
    Allocate( tmp(nm,nm) )
    nmax=100
    z = dcmplx(E,1.0d-3)
    Id=0.0d0
    tmp=0.0d0
    do i=1,nm
        Id(i,i)=1.0d0
        tmp(i,i)=dcmplx(0.0d0,1.0d0)
    enddo
    H_BB = H00
    H_10 = H10
    H_01 = TRANSPOSE( CONJG( H_10 ) )
    H_SS = H00
    do i = 1, nmax
        A = z*S00 - H_BB
        !
        call invert(A,nm)
        !
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_SS = H_SS + C
        H_BB = H_BB + C
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)  
        H_10 = C    
        H_BB = H_BB + A
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_01 = C 
        ! NORM --> inspect the diagonal of A
        error=0.0d0
        DO k=1,nm
            DO j=1,nm
                error=error+sqrt(aimag(C(k,j))**2+Dble(C(k,j))**2)
            END DO
        END DO	
        tmp=H_SS
        IF ( abs(error) < TOL ) THEN	
            EXIT
        ELSE
        END IF
        IF (i .EQ. nmax) THEN
            write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error
            !call abort
            H_SS=H00
            H_BB=H00
        END IF
    enddo
    G00 = z*S00 - H_SS
    !
    call invert(G00,nm)
    !
    GBB = z*S00 - H_BB
    !
    call invert(GBB,nm)
    !
    Deallocate( tmp )
    Deallocate( A )
    Deallocate( B )
    Deallocate( C )
    Deallocate( H_BB )
    Deallocate( H_SS )
    Deallocate( H_01 )
    Deallocate( H_10 )
    Deallocate( Id )
end subroutine sancho




end module gw_mod
