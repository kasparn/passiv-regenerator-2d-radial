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

module grid
use datatypes
implicit none

contains

function gridConstant( geo )
type(geom),intent(inout) :: geo
integer :: nx,nr_sf,nr_w,Nr,i,gridConstant

  nx = geo%nx
  nr_sf = geo%nr_sf
  nr_w = geo%nr_w
 allocate(geo%dx(nx))
 allocate(geo%dr(nr_sf+nr_w),geo%r(nr_sf+nr_w))
 allocate(geo%rdn(nr_sf+nr_w),geo%rup(nr_sf+nr_w))
 !discretization in the x-direction
 !Can be made heterogeneous if wanted

 geo%dx = geo%L / geo%nx

 !discretization in the radial direction
 !is currently heterogeneous in the sense that
 !the solid/fluid domain has one resolution whereas the wall
 !domain has another
 
 geo%dr(1:nr_sf) = geo%Rad / nr_sf 
 geo%dr((nr_sf+1):(nr_sf+nr_w)) = geo%Hw / nr_w 

Nr = nr_sf + nr_w 

geo%r(1) = 0.5*geo%dr(1) 
geo%rup(1) = geo%dr(1) 
geo%rdn(1) = 0 

do i=2,Nr
    geo%r(i) = geo%r(i-1) + 0.5 * ( geo%dr(i-1) + geo%dr(i) ) 
    geo%rup(i) = geo%rup(i-1) + geo%dr(i) 
    geo%rdn(i) = geo%rdn(i-1) + geo%dr(i)      
enddo

!no error
gridConstant = 0

end function gridConstant

end module grid