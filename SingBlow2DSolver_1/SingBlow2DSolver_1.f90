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

program SingBlow2DSolver_1
use bedcorrelations
use constants
use datatypes
use flowprofiles
use fluidproperties
use grid
use solidproperties
use solver
use wallproperties
use tempprofiles
use radialprofiles
use io
implicit none

type(geom) :: geo
type(opCond) :: op
type(Temp) :: tR
type(tim) :: ti
type(propSetup) :: props
type(profSetup) :: prof
real :: CFL,rhof
real :: t1,t2
character(len=1000):: fileInput,fileOutputT,fileOutputQ
character(len=1000),parameter :: files='io.txt'

!::Get the file names
call getFileNames( files, fileInput, fileOutputT, fileOutputQ )

!::Load the model parameters
call read_input( fileInput, geo, op, ti, props, prof )

!::Load the thermal properties of the three domains
call loadProps( props )

!::Load the fluid flow profile
call loadFlowProfile( prof )

!::Load the inlet temperature profile
call loadTinProfile ( prof, op )

!Still a bit of a hack..
rhof = 1000
!::Determine the number of timesteps based on the CFL_max number
ti%n_timesteps = int(abs( prof%mdot0 * ti%t_end * geo%nx / &
                ( ti%CFL_max * rhof * pi * geo%Rad**2 * geo%por * geo%L ) )) + 1

!::Make sure there is an even number of timesteps
if ( mod( ti%n_timesteps, 2 ) .ne. 0 ) ti%n_timesteps = ti%n_timesteps + 1

write(*,*) 'n_timesteps',ti%n_timesteps

CFL = prof%mdot0 * ti%t_end * geo%nx / ( rhof * pi * geo%Rad**2 * geo%por * ti%n_timesteps * geo%L )

if ( CFL .gt. 0.5 ) then
    write(*,*) 'The Courant number, CFL, is >0.5 and you should be careful.'
    write(*,*) 'It is recommended to either decrease the number of mesh points or to increase the number of timesteps'
endif
ti%dt = ti%t_end/ti%n_timesteps


write(*,*) 'CFL',CFL,ti%dt

!::Call the solver
call cpu_time(t1)
call SingBSolver( geo, ti, op, tR, getFlowProfile, getTempProfile,Spheres,gridConstant, &
                  getFluidProps,getSolidProps,getWallProps, getRadialProfile,props, prof )
call cpu_time(t2)
write(*,*) 'total cpu time',t2-t1


call write_solution( fileOutputT, fileOutputQ, tR,geo,ti, op )
!::Do post-processing

end program SingBlow2DSolver_1

