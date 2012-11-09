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

MODULE BEDCORRELATIONS
use datatypes

implicit none



contains

!::returns the Nusselt number, pressure gradient (dpdx), 
!::specific surface area (as) and the hydraulic diamter, dh.
!::This implementation is for packed spheres
function Spheres( geo, op, fl, sl, ti, corr )
type(geom),intent(in) :: geo
type(opCond),intent(in) :: op
type(fluidProps),intent(inout) :: fl
type(solidProps),intent(inout) :: sl
type(tim),intent(in) :: ti
type(corrParamsOut),intent(inout) :: corr
integer :: nx,nr_sf
real :: disp,por,dpar,tau,a0,f0,A,B,alp,bet,as,dh
real,dimension(:,:),allocatable :: Rep,Pr,Bi,Fo,Chi,kstat,kdispAx,kdispRad
integer :: Spheres

nx = geo%nx
nr_sf = geo%nr_sf

allocate( Rep(nx,nr_sf),Pr(nx,nr_sf),Bi(nx,nr_sf),Fo(nx,nr_sf),Chi(nx,nr_sf))
allocate( kstat(nx,nr_sf),kdispAx(nx,nr_sf), kdispRad(nx,nr_sf) )

tau = ti%t_end
!hydraulic diameter
corr%dh = 2./3. * geo%por / ( 1. - geo%por ) * geo%dpar 

!specific heat transfer surface area
corr%as = 6. * ( 1 - geo%por ) / geo%dpar 

!hydraulic diameter based Re-number
corr%Ref = abs( fl%rhof * op%u * geo%por * corr%dh / fl%muf )
!particle diameter based Re-number
Rep = abs( fl%rhof * op%u * geo%por * geo%dpar / fl%muf )
!fluid Prandtl number
Pr = fl%cf * fl%muf / fl%kf 

!from Wakao and Kaguei 1982, based on particle diameter
corr%Nu = 2. + 1.1 * Rep**0.6 * Pr**(1./3.) 
!This if for corrected for internal temperature gradients (see Engelbrecht
!et al 2006) (totally based on Engelbrecht's model implementation)
! Calculate Biot number and the heat transfer degradation factor
Bi = ( corr%Nu * fl%kf )/( 2 * sl%ks ) 
!the Fourier number
Fo = fl%kf * tau / ( fl%rhof * fl%cf * (geo%dpar/2)**2) 
chi = Fo*exp(0.246196-0.84878*log(Fo)-0.05639*(log(Fo))**2) 
corr%Nu = corr%Nu / ( 1 + chi * Bi/5)  !corrected heat transfer coefficient 

!Taken explicitly from the model of Engelbrecht
!Calculate static conductivity according to Hadley (1986) and axial
!dispersion from Kaviany (1995)
if ( geo%por .lt. 0.0827 ) then
    a0 = 10**(-4.898 * geo%por) 
elseif ( geo%por .lt. 0.298 ) then
    a0 = 10**(-.405 - 3.154 * (geo%por - 0.0827) ) 
elseif ( geo%por .lt. 0.58 ) then
    a0 = 10**( -1.084 - 6.778 * ( geo%por - 0.298) ) 
endif
f0 = 0.8 + 0.1 * geo%por 
kstat = fl%kf * ( (1-a0) * ( geo%por * f0 + sl%ks/fl%kf * (1-geo%por*f0) ) / &
        (1-geo%por*(1-f0) + sl%ks/fl%kf * geo%por * (1-f0)) + a0 * (2*(sl%ks/fl%kf)**2 * &
        (1-geo%por) + (1+2*geo%por) * sl%ks/fl%kf) / ((2+geo%por)*sl%ks/fl%kf + 1 - geo%por))  !W/m2-K, Hadley (1986) from Kaviany

!::Based on Kaviany 1995 AND Delgado (2006). The latter states that the radial
!::dispersion is about 1/5 of the axial when Ref > 10. When Ref<1 they are both
!::equal to kf and we assume linearity for 1 <= Ref < 10
where ( corr%Ref .ge. 10 )
    kdispAx = .75 * geo%por * Pr * corr%Ref  !W/m2-K, effective conduction in fluid due to axial dispersion
    kdispRad = kdispAx * 1./5.
endwhere

where( corr%Ref .gt. 1 .AND. corr%Ref .lt. 10 )    
    
    kdispAx = (0.75 * geo%por * Pr * 10. - 1) / ( 10. - 1.) * ( corr%Ref - 1 ) + 1;
    
    kdispRad = (1./5. * 0.75 * geo%por * Pr * 10. - 1) / ( 10. - 1.) * ( corr%Ref - 1 ) + 1;
endwhere

where( corr%Ref .le. 1 )
    kdispAx = 1.
    kdispRad = 1.
endwhere

kdispAx = kdispAx * fl%kf
kdispRad = kdispRad * fl%kf


!the effective solid thermal conductivity
!the disp parameters should be between 0 and 1 
!disp = 1 means the full dispersion formulation is taken into account
!disp = 0 means dispersion is neglected
sl%kstat = op%dispAx * kstat + (1-op%dispAx) * sl%ks
!fluid effective thermal conductivity (axially)
fl%kdispAx = op%dispAx * kdispAx + (1-op%dispAx) * fl%kf 
fl%kdispRad = op%dispRad * kdispRad + (1-op%dispRad) * fl%kf 
!convert it so that it is based on the hydraulic diameter
corr%Nu = corr%dh/geo%dpar * corr%Nu

!pressure gradient, the Ergun equation. It is important to realize the
!difference between u, the pore velocity, and u*por=v, the open area
!velocity
A = 180 !constants
B = 1.8 

alp = (1-geo%por)**2/geo%por**2 
bet = (1-geo%por)/geo%por**3 

corr%dpdx = A * alp * sum(fl%muf,2)/nr_sf / geo%dpar**2 * sum(abs(op%u),2)/nr_sf * geo%por + &
       B * bet * sum(fl%rhof,2)/nr_sf / geo%dpar * ( sum(abs(op%u),2)/nr_sf * geo%por )**2 

deallocate( Rep,Pr,Bi,Fo,Chi,kstat,kdispAx,kdispRad )


!no error
Spheres = 0

end function spheres


end module BEDCORRELATIONS