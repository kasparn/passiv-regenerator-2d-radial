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

module util

implicit none

contains

!::Performs linear interpolation in the x,y data set
!::For each of the values in xx (with size n)
subroutine interp1( x, y, xx, nx, n, yy )
real,dimension(nx),intent(in) :: x,y
real,dimension(n),intent(in) :: xx
integer,intent(in) :: nx,n
real,dimension(n),intent(out) :: yy
integer :: i,ind
real :: lin

do i=1,n    
    ind = locate( x, xx(i), nx )    
    lin = (xx(i) - x(ind))/(x(ind+1)-x(ind))
    yy(i) = (1-lin) * y(ind) + lin * y(ind+1)
enddo

end subroutine interp1

!::Based on the NR edition of this method
!::Simply finds the index i where:
!:: xx(i) <= x <= xx(i+1)
!::Returns 1 or n if out of lower or upper bounds, respectively
FUNCTION locate(xx,x,n)
	
	IMPLICIT NONE
	REAL, DIMENSION(n), INTENT(IN) :: xx
	REAL, INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,jl,jm,ju
	LOGICAL :: ascnd
	
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
END FUNCTION locate


end module util