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

module radialprofiles
use datatypes
use constants
use util
implicit none


contains

!::
!::Sets the radial velocity profile so far based on
!::a constant porosity but this could be changed
!::It is also remarked that the density of the fluid
!::is assumed constant - otherwise the normalization should
!::be re-done to take that into account in order to
!::conserve the fluid mass
!::
function getRadialProfile( op, geo )
type(opCond),intent(inout) :: op
type(geom),intent(in) :: geo
integer :: getRadialProfile
real,dimension(:),allocatable :: umean
real :: C1,kk
real,parameter :: c=0.05
integer :: i
!::Does nothing at the moment, should be implemented though.
!allocate(umean(geo%nx))
!
!umean = op%u(:,1)
!!::Normalization constant
!write(*,*) geo%dpar,-c*geo%dpar,-geo%Rad-c*geo%dpar
!C1 = 1/geo%Rad *( c*geo%dpar*( log(-c*geo%dpar) - log(-geo%Rad-c*geo%dpar) ) + geo%Rad)
!write(*,*) 'C1',C1
!stop
!do i=1,geo%nr_sf
!!::The profile that is constant and then drops rapidly
!!::to zero at R over a few sphere diameters
!    kk = 1 / ( 1 + c * geo%dpar / ( geo%Rad - geo%r(i) ) )
!    
!    op%u(:,i) = kk / C1 * umean
!enddo
!
!
!deallocate(umean)
!no error
getRadialProfile = 0
end function getRadialProfile


end module radialprofiles