!::Copyright (c) 2012, Kaspar Kirstein Nielsen (kaki@dtu.dk, Technical University of Denmark)
!::All rights reserved.
!::
!::Redistribution and use in source and binary forms, with or without modification, 
!::are permitted provided that the following conditions are met:
!::
!::Redistributions of source code must retain the above copyright notice, 
!::this list of conditions and the following disclaimer.
!::Redistributions in binary form must reproduce the above copyright notice, this 
!::list of conditions and the following disclaimer in the documentation and/or other materials 
!::provided with the distribution.
!::Neither the name of the Technical University of Denmark nor the names of its contributors 
!::may be used to endorse or promote products derived from this software without specific prior 
!::written permission.
!::THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
!::OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
!::AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
!::CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
!::DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
!::OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!:: STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
!::SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


include 'mkl_dss.f90'
module direct_solver
use mkl_dss

implicit none
contains
!::
!::This routine uses the Intel Math Kernel Library
!::to solve the equation Ax = B through a direct sparse
!::algorithm. A is assumed to be a sparse matrix on vector form,
!::i.e. length(A) = nNonZero and A is one-dimensional.
!::It represents the full matrix, which has a size of (n,n)
!::B is a column vector with size n and so is x
!::The results are returned via X and A and B are overwritten
!::cols has size nNonZero and contains the column indices
!::that correspond to the values in A
!::rows has size n+1 and contains the indices in A for the first
!::non-zero element in each row and the last element is equal to
!::nNonZero+1
!::i.e.: cols(i) = column index of the i'th value in A
!::rows(j) = index in A where the first non-zero element of the j'th row
!::is located
!::
subroutine solve_lin_direct( A, B, X, cols,rows,n, nNonZero )
real,intent(inout),dimension(nNonZero) :: A
real,intent(inout),dimension(n) :: B,X
integer,intent(in) :: n,nNonZero
integer,intent(in),dimension(n+1) :: rows
integer,intent(in),dimension(nNonZero) :: cols
integer,parameter :: nrhs=1 !::There is only a single column on the right hand side
integer :: error,i
type(mkl_dss_handle) :: handle ! allocate storage for the solver handle.
integer perm(1)

perm(1) = 0

! initialize the solver.
error = dss_create( handle, mkl_dss_defaults )
if (error /= mkl_dss_success) then
write(*,*) 'Error when calling dss_create',error
stop
endif

! define the non-zero structure of the matrix.
error = dss_define_structure( handle, mkl_dss_symmetric_structure, rows, n, n, cols, nNonZero )
if (error /= mkl_dss_success) then
write(*,*) 'Error when calling dss_define_structure',error
stop
endif

! reorder the matrix.
error = dss_reorder( handle, mkl_dss_defaults, perm )
if (error /= mkl_dss_success) then
write(*,*) 'Error when calling dss_reorder',error
stop
endif

! factor the matrix.
error = dss_factor_real( handle, mkl_dss_defaults, A )
if (error /= mkl_dss_success) then
write(*,*) 'Error when calling dss_factor_real'
stop
endif

error = dss_solve_real(handle, mkl_dss_defaults, B, 1, X )
if (error /= mkl_dss_success) then
write(*,*) 'Error when calling dss_solve_real',error
stop
endif

! deallocate solver storage and various local arrays.
error = dss_delete( handle, mkl_dss_defaults )
if (error /= mkl_dss_success ) then
write(*,*) 'Error calling dss_delete',error
stop
endif

end subroutine solve_lin_direct

end module direct_solver
