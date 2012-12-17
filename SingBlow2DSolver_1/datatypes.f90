
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
module datatypes
implicit none

type fluidProps
real,dimension(:,:),allocatable :: kf,cf,rhof,muf,kdispAx,kdispRad
real,dimension(:),allocatable :: kfBdry
endtype fluidProps

type solidProps
real,dimension(:,:),allocatable :: ks,cs,rhos,kstat
endtype solidProps

type wallProps
real,dimension(:,:),allocatable :: kw,cw,rhow
endtype wallProps

type corrParamsOut
real,dimension(:,:),allocatable :: Nu,Ref
real,dimension(:),allocatable :: dpdx
real :: as, dh
endtype corrParamsOut

type geom
real,dimension(:),allocatable :: dx,dr,r,rup,rdn
integer :: nx,nr_sf,nr_w
real :: L, Rad, Hw,dpar,por
integer :: geomType
endtype geom

type Temp
real,dimension(:,:,:),allocatable :: Tf,Ts,Tw
real,dimension(:),allocatable ::TfHout,TfCout,qhot,qcold
endtype Temp

type tim
real :: dt,t_end,CFL_max
integer :: n_timesteps
real,dimension(:),allocatable :: t
endtype tim

type opCond
real :: axial,walls,wallsSol,dispAx,dispRad,Nuscl,viscDiss,Tc,Th
real :: Tf0,Ts0,Tw0,hw
real,dimension(:,:),allocatable :: u
integer :: blowMode,kfBdry
integer :: maxIte
real :: convTol
real,dimension(:),allocatable :: convErrHot,convErrCold
endtype opCond

type propSetup
integer :: wallType,fluidType,solidType
character(len=1000) :: fluidFile,solidFile,wallFile
integer :: kfNT,rhofNT,cfNT,mufNT
integer :: ksNT,rhosNT,csNT
integer :: kwNT,rhowNT,cwNT
real,dimension(:,:),allocatable :: kf,rhof,cf,muf
real,dimension(:,:),allocatable :: ks,rhos,cs
real,dimension(:,:),allocatable :: kw,rhow,cw
endtype propSetup

type profSetup
character(len=1000) :: profileFile,TinFile
integer :: flowProfile,TinProfile
real,dimension(:,:),allocatable :: mdotIn,Tin
integer :: nT,TinNt
real :: mdot0
endtype profSetup

end module datatypes

