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
module wallproperties
use datatypes
use constants
use util
implicit none

contains

function getWallProps( Tw, nx, nr, wlIn, props )
type(wallProps),intent(inout) :: wlIn
real,intent(in),dimension(nx,nr) :: Tw
integer,intent(in) :: nx,nr
type(propSetup),intent(in) :: props
integer :: ret,getWallProps

if ( props%wallType .eq. domTypeInternal ) then
    !::Simply call the Gd properties
    !::This may be extended with more than just water
    ret = PlasticSimple( Tw, nx, nr, wlIn )
elseif ( props%wallType .eq. domTypeUserDef ) then
    !::Load user defined properties
    if ( props%kwNT .eq. 1 ) then
        wlIn%kw = props%kw(2,1)
    else
        !::Make interpolation
        call interp1( props%kw(1,:), props%kw(2,:), Tw, props%kwNT, nx*nr, wlIn%kw )
    endif
    
    if ( props%rhowNT .eq. 1 ) then
        wlIn%rhow = props%rhow(2,1)
    else
        !::Make interpolation
        call interp1( props%rhow(1,:), props%rhow(2,:), Tw, props%rhowNT, nx*nr, wlIn%rhow )
    endif
    
    if ( props%cwNT .eq. 1 ) then
        wlIn%cw = props%cw(2,1)
    else
        !::Make interpolation
        call interp1( props%cw(1,:), props%cw(2,:), Tw, props%cwNT, nx*nr, wlIn%cw )
    endif    
endif
!::No error
getWallProps = 0

end function getWallProps

!returns the thermal properties of the regenerator
!housing.
!This implementation assumes constant properties of a plastic
!(similar to nylon)
function plasticsimple( Tw, nx, nr, wlIn )
type(wallProps),intent(inout) :: wlIn
real,dimension(nx,nr),intent(in) :: Tw
integer,intent(in) :: nx,nr
integer :: plasticSimple

wlIn%kw = 0.25;

wlIn%cw = 1500;

wlIn%rhow = 1000;
!no error
plasticSimple = 0

end function plasticsimple

!::Stainless steel (constant) properties
function StainlessSteel( Tw, nx, nr, wlIn )
type(wallProps),intent(inout) :: wlIn
real,dimension(nx,nr),intent(in) :: Tw
integer,intent(in) :: nx,nr
integer :: StainlessSteel

wlIn%kw(:,:) = 16;

wlIn%cw(:,:) = 450;

wlIn%rhow(:,:) = 7800;
!no error
StainlessSteel = 0

end function StainlessSteel

!::Aluminum (constant) properties
function Aluminum( Tw, nx, nr,wlIn )
type(wallProps),intent(inout) :: wlIn
real,dimension(nx,nr),intent(in) :: Tw
integer,intent(in) :: nx,nr
integer :: Aluminum

wlIn%kw(:,:) = 237;

wlIn%cw(:,:) = 900;

wlIn%rhow(:,:) = 2700;
!no error
Aluminum = 0
end function Aluminum

end module wallproperties