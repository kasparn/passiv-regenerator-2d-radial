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

module solidproperties
use datatypes
use constants
use util
implicit none

contains

function getSolidProps( Ts, nx, nr, slIn, props )
type(solidProps),intent(inout) :: slIn
real,intent(in),dimension(nx,nr) :: Ts
integer,intent(in) :: nx,nr
type(propSetup),intent(in) :: props
integer :: ret,getSolidProps

if ( props%solidType .eq. domTypeInternal ) then
    !::Simply call the Gd properties
    !::This may be extended with more than just water
    ret = GdConst( Ts, nx, nr, slIn )
elseif ( props%solidType .eq. domTypeUserDef ) then
    !::Load user defined properties
    if ( props%ksNT .eq. 1 ) then
        slIn%ks = props%ks(2,1)
    else
        !::Make interpolation
        call interp1( props%ks(1,:), props%ks(2,:), Ts, props%ksNT, nx*nr, slIn%ks )
    endif
    
    if ( props%rhosNT .eq. 1 ) then
        slIn%rhos = props%rhos(2,1)
    else
        !::Make interpolation
        call interp1( props%rhos(1,:), props%rhos(2,:), Ts, props%rhosNT, nx*nr, slIn%rhos )
    endif
    
    if ( props%csNT .eq. 1 ) then
        slIn%cs = props%cs(2,1)
    else
        !::Make interpolation
        call interp1( props%cs(1,:), props%cs(2,:), Ts, props%csNT, nx*nr, slIn%cs )
    endif    
endif
!::No error
getSolidProps = 0

end function getSolidProps

!Returns the solid properties as a function of T.
!Here, they are all assumed constant and Gd-like
function GdConst( Ts, nx, nr, slIn )
type(solidProps),intent(inout) :: slIn
real,intent(in),dimension(nx,nr) :: Ts
integer,intent(in) :: nx,nr
integer :: GdConst

slIn%ks = 10.5

slIn%cs = 300

slIn%rhos = 7900
!no error
GdConst = 0
end function GdConst

end module solidproperties