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

use matrix_mod, only : cMatrix

implicit none

private

integer, parameter :: dp=8
complex(dp),allocatable :: Ham(:),S(:)
integer,allocatable :: row(:),col(:)
integer :: nnz

public :: deviceHam_load_COOmatrix, deviceHam_build_blocks, deviceHam_free



contains

subroutine deviceHam_free()
  if (allocated(Ham)) deallocate(Ham)
  if (allocated(S)) deallocate(S)
  if (allocated(row)) deallocate(row)
  if (allocated(col)) deallocate(col)
end subroutine deviceHam_free


subroutine deviceHam_load_COOmatrix(fname,overlap)
character(len=*), intent(in)        :: fname !! input text file name
logical, intent(in) :: overlap !! read overlap matrix

end subroutine deviceHam_load_COOmatrix



subroutine deviceHam_build_blocks


end subroutine deviceHam_build_blocks

end module deviceHam_mod
