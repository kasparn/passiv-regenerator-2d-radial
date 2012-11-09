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

module solver
use datatypes
use wallproperties
use solidproperties
use grid
use fluidproperties
use solver_call_sparse
use constants
implicit none

contains

!::
!::Main call that takes care of setting up grid, properties etc and then
!::does the timestep loop
!::
subroutine SingBSolver( geo, ti, op, tR, flowFunc,tempFunc,&
bedFunc,gridFunc, fluidFunc,solidFunc,wallFunc, velocityFunc,props, prof )
type(Temp),intent(inout) :: tR
type(geom),intent(inout) :: geo
type(tim),intent(inout) :: ti
type(opCond),intent(inout) :: op
type(fluidProps) :: flProps
type(solidProps) :: slProps
type(wallProps) :: wlProps
type(corrParamsOut) :: bedCorr
type(propSetup) :: props
type(profSetup) :: prof
real,dimension(:,:),allocatable :: Qvisc,h
real,dimension(:),allocatable :: hw,mdot,Tc,Th
integer :: nx,nr_sf,nr_w,ret,n_timesteps,i,j,jj,cnt
integer,external :: bedFunc,gridFunc,flowFunc,tempFunc
integer,external :: solidFunc,wallFunc,fluidFunc,velocityFunc
real :: Ac,perc,vol,mdot1,mdot2
integer :: cnt1,cnt2


!!First setup the grid
ret = gridFunc( geo )
!dV = dx(1) * 2 * pi * r * dr 
nx = geo%nx  !grid points in the x-direction
nr_sf = geo%nr_sf !grid points in the radial direction in the fluid and solid domain
nr_w = geo%nr_w !grid points in the radial direction in the wall domain
n_timesteps = ti%n_timesteps
!duration of each timestep
ti%dt = ti%t_end/n_timesteps 
!temperatures of the three substances
allocate( tR%Tf(nx,nr_sf,n_timesteps+1) ,tR%Ts(nx,nr_sf,n_timesteps+1)  )
allocate( tR%Tw(nx,nr_w,n_timesteps+1) )
!average outlet temperature as a function of time
allocate( tR%TfHout(n_timesteps),tR%TfCout(n_timesteps) ) 
allocate( ti%t(n_timesteps+1) )

allocate( tR%qhot(op%maxIte),tR%qcold(op%maxIte) )
allocate( op%convErrHot(op%maxIte),op%convErrCold(op%maxIte) )
tR%qhot = 0
tR%qcold = 0
op%convErrHot = 0
op%convErrCold = 0

!initial temperature
if ( op%blowMode .eq. SingleBlowMode ) then
    tR%TfHout(1) = op%Tf0 
    tR%TfCout(1) = op%Tf0
    tR%Tf(:,:,1) = op%Tf0 
    tR%Ts(:,:,1) = op%Ts0 
    tR%Tw(:,:,1) = op%Tw0 
elseif ( op%blowMode .eq. PeriodicMode ) then
    !::Impose a linear temperature profile
    do i=1,geo%nx
        tR%Tf(i,:,1) = 1.*i/(1.*geo%nx) * ( op%Th - op%Tc ) + op%Tc
        tR%Ts(i,:,1) = tR%Tf(i,:,1)
        tR%Tw(i,:,1) = tR%Tf(i,:,1)
    enddo
    
endif
!Allocate the property arrays
allocate( flProps%kf(nx,nr_sf),flProps%kdispAx(nx,nr_sf),flProps%kdispRad(nx,nr_sf) )
allocate( flProps%cf(nx,nr_sf),flProps%rhof(nx,nr_sf),flProps%muf(nx,nr_sf) )
allocate( slProps%ks(nx,nr_sf),slProps%kstat(nx,nr_sf),slProps%cs(nx,nr_sf),slProps%rhos(nx,nr_sf) )
allocate( wlProps%kw(nx,nr_w),wlProps%cw(nx,nr_w),wlProps%rhow(nx,nr_w) )

allocate( mdot(n_timesteps), Tc(n_timesteps), Th(n_timesteps) )
!mass flow rate profile (as a function of time)
ret = flowFunc( ti%dt, n_timesteps, prof, op, mdot ) 

!Inlet temperatures as a function of time
ret = tempFunc( ti%dt, n_timesteps, prof, op, Tc, Th )

Ac = pi * geo%Rad**2 !bed cross sectional area

vol = geo%L * Ac 

!HT coefficient between wall and fluid
allocate( hw(nx) ) 
hw = op%hw 

allocate( op%u(nx,nr_sf) )
allocate( bedCorr%Nu(nx,nr_sf),bedCorr%Ref(nx,nr_sf) )
allocate( bedCorr%dpdx(nx) )

!Viscous dissipation in W at each i,j
allocate( Qvisc( nx, nr_sf ), h(nx, nr_sf) )
Qvisc(:,:) = 0.   

if ( op%blowMode .eq. SingleBlowMode ) then
    op%maxIte = 1
endif

!::Iteration loop
!::If the mode is set to SingleBlow only one iteration is performed
do jj=1,op%maxIte
    !::At all iterations except the first the final state is transfered to
    !::be the initial condition of the next iteration
    if ( jj .gt. 1 ) then
        tR%Tf(:,:,1) = tR%Tf(:,:,n_timesteps+1)
        tR%Ts(:,:,1) = tR%Ts(:,:,n_timesteps+1)
        tR%Tw(:,:,1) = tR%Tw(:,:,n_timesteps+1)
    endif
    
    cnt1 = 0
    cnt2 = 0
    mdot1 = 0
    mdot2 = 0
    
    perc = 0 
    cnt = 0 
    do i=1,n_timesteps
        if ( cnt .ge. n_timesteps/10 ) then
            perc = perc + 10 
            write(*,*) perc
            cnt = 0 
        endif    
        cnt = cnt + 1     
        !get the thermal properties which may be a function of T
        !thermal properties of the fluid
        ret = fluidFunc( tR%Tf(:,:,i),nx,nr_sf, flProps, props ) 
        !Thermal properties of the solid
        ret = solidFunc( tR%Ts(:,:,i), nx, nr_sf, slProps, props ) 
        
        !thermal properties of the wall    
        ret = wallFunc( tR%Tw(:,:,i), nx, nr_w, wlProps, props ) 

        !pore fluid velocity    
        op%u = mdot(i) / ( flProps%rhof * Ac * geo%por ) 
        !Radial velocity profile, assuming constant rhof !!!!
        !ret = velocityFunc( op, geo )
        !Nusselt number and pressure drop
        ret = bedFunc( geo, op, flProps, slProps, ti, bedCorr )
        bedCorr%Nu = bedCorr%Nu * op%Nuscl 
           
            
        do j=1,nr_sf
            !viscous dissipation at i and per m^3
            Qvisc(:,j) = op%viscDiss * abs( bedCorr%dpdx * geo%dx * mdot(i) / &
            flProps%rhof(:,1) ) / ( geo%dx * pi * geo%Rad**2 ) 
        enddo
        !heat transfer coefficient
        
        h = bedCorr%Nu * flProps%kf / bedCorr%dh
        !NTU = mean(mean( h * as * vol / ( mdot0 * cf ) ) ) 
        !Ref = mean(mean(Ref)) 
        !call the timestep
        
        call timestep_sp( geo%dx, geo%dr, geo%r, geo%rup, geo%rdn, &
                       geo%nx, geo%nr_sf, geo%nr_w,&
                       flProps%rhof,flProps%kdispAx,flProps%kdispRad,flProps%cf,&
                       slProps%rhos,slProps%kstat,slProps%cs,&
                       wlProps%rhow,wlProps%kw,wlProps%cw,&
                       h,bedCorr%as,ti%dt, geo%por, &
                       tR%Tf(:,:,i), tR%Ts(:,:,i), tR%Tw(:,:,i), &
                       op%u, hw, Qvisc, Tc(i), Th(i), &
                       op%axial, op%walls, op%wallsSol,tR%Tf(:,:,i+1), tR%Ts(:,:,i+1), tR%Tw(:,:,i+1)) 
       

    !get the <Tf(nx)> as a function of time
        tR%TfHout(i) = sum( tR%Tf(nx,:,i+1) * 2 * pi * geo%r(1:nr_sf) * geo%dr(1:nr_sf) ) / Ac  
        tR%TfCout(i) = sum( tR%Tf(1,:,i+1) * 2 * pi * geo%r(1:nr_sf) * geo%dr(1:nr_sf) ) / Ac  
        ti%t(i+1) = ti%t(i) + ti%dt
        if ( mdot(i) .gt. 0 ) then
            tR%qhot(jj) = tR%qhot(jj) + mdot(i) * (tR%TfHout(i) - op%Tc) * flProps%cf(nx,1)*ti%dt
            cnt1 = cnt1 + 1
            mdot1 = mdot1 + mdot(i)
        elseif ( mdot(i) .lt. 0 ) then
            tR%qcold(jj) = tR%qcold(jj) + mdot(i) * (op%Th - tR%TfCout(i)) * flProps%cf(1,1)*ti%dt
            cnt2 = cnt2 + 1
            mdot2 = mdot2 + mdot(i)
        endif
    enddo
    write(*,*) 'Flow balanced?',cnt1,cnt2
    write(*,*) mdot1,mdot2
    
    if ( jj .gt. 1 ) then
        op%convErrHot(jj) = abs( (tR%qhot(jj) - tR%qhot(jj-1))/tR%qhot(jj-1) )
        op%convErrCold(jj) = abs( (tR%qcold(jj) - tR%qcold(jj-1))/tR%qcold(jj-1) )
        
        !::The convergence is based on the relative change in qhot and qcold
        !::It is important to realize that they should not, in general, be 
        !::equal as long as viscous dissipation or other heat generating
        !::mechanisms are involved.
        write(*,*) 'Err',op%convErrHot(jj),op%convErrCold(jj)
        write(*,*) 'qhot,qcold',tR%qhot(jj),tR%qcold(jj)
        if ( op%convErrHot(jj) .le. op%convTol .AND. op%convErrCold(jj) .le. op%convTol ) then
            write(*,*) 'Convergence reached'
            exit
        elseif ( jj .eq. op%maxIte ) then
            write(*,*) 'Convergence was not reached after max. nr. of iterations'
        endif
        
    endif
    write(*,*) 'Iteration',jj,op%maxIte
    
enddo

open (14, file='Tf.dat',	&
			           status='unknown', form='unformatted',	&
			           access='direct', recl=2*geo%nx*geo%nr_sf)

do i=1,n_timesteps
    write(14,rec=i) tR%Tf(:,:,i)
enddo

close(14)
end subroutine SingBSolver
end module solver