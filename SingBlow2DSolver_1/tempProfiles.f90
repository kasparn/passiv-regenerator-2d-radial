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

module tempprofiles
use datatypes
use constants
use util
implicit none

contains

function getTempProfile( dt, n_timesteps, prof, op, TCout,THout )
real,intent(in) :: dt
integer,intent(in) :: n_timesteps
type(profSetup),intent(in) :: prof
type(opCond),intent(in) :: op
integer :: getTempProfile,i
real,dimension(n_timesteps),intent(inout) :: TCout,THout
real,dimension(n_timesteps) :: t


if ( prof%TinProfile .ne. TinputFile ) then
    !Default use the constant mdot function
    getTempProfile = TConst( op%Tc, dt, n_timesteps, TCout )
    getTempProfile = TConst( op%Th, dt, n_timesteps, THout )
elseif ( prof%TinProfile .eq. TinputFile ) then
   
    t(1) = 0
    do i=2,n_timesteps
        t(i) = t(i-1) + dt
    enddo
    !Interpolate in a table provided externally
    call interp1( prof%Tin(1,:), prof%Tin(2,:), t, prof%TinNt, n_timesteps, TCout )    
    THout = TCout
endif
!No error
getTempProfile = 0

end function getTempProfile

!simply returns a constant Tinlet
function Tconst( T0, dt, n_timesteps, Tout )
integer :: Tconst
integer,intent(in) :: n_timesteps
real,intent(in) :: dt,T0
real,dimension(n_timesteps),intent(inout) :: Tout

Tout = T0
!no error
Tconst = 0

end function Tconst

end module tempprofiles