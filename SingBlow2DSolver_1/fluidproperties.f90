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

module fluidproperties
use datatypes
use constants
use util
implicit none
contains

function getFluidProps( Tf, nx, nr, flIn, props )
type(fluidProps),intent(inout) :: flIn
real,intent(in),dimension(nx,nr) :: Tf
integer,intent(in) :: nx,nr
type(propSetup),intent(in) :: props
integer :: ret,getFluidProps

if ( props%fluidType .eq. domTypeInternal ) then
    !::Simply call the water properties
    !::This may be extended with more than just water
    ret = water( Tf, nx, nr, flIn )
elseif ( props%fluidType .eq. domTypeUserDef ) then
    !::Load user defined properties
    if ( props%kfNT .eq. 1 ) then
        flIn%kf = props%kf(2,1)
    else
        !::Make interpolation
        call interp1( props%kf(1,:), props%kf(2,:), Tf, props%kfNT, nx*nr, flIn%kf )
    endif
    
    if ( props%rhofNT .eq. 1 ) then
        flIn%rhof = props%rhof(2,1)
    else
        !::Make interpolation
        call interp1( props%rhof(1,:), props%rhof(2,:), Tf, props%rhofNT, nx*nr, flIn%rhof )
    endif
    
    if ( props%cfNT .eq. 1 ) then
        flIn%cf = props%cf(2,1)
    else
        !::Make interpolation
        call interp1( props%cf(1,:), props%cf(2,:), Tf, props%cfNT, nx*nr, flIn%cf )
    endif
    
    if ( props%mufNT .eq. 1 ) then
        flIn%muf = props%muf(2,1)
    else
        !::Make interpolation
        call interp1( props%muf(1,:), props%muf(2,:), Tf, props%mufNT, nx*nr, flIn%muf )
    endif
endif
!::No error
getFluidProps = 0

end function getFluidProps



!Returns the thermal properties of water as a function of the temperature
!In this implementation the properties are constant
function water( Tf, nx, nr, flIn )
type(fluidProps),intent(inout) :: flIn
real,intent(in),dimension(nx,nr) :: Tf
integer,intent(in) :: nx,nr
integer :: water

flIn%kf = 0.6;

flIn%cf = 4200;

flIn%rhof = 1000;

flIn%muf = 0.001;
!no error
water = 0
end function water

end module fluidproperties