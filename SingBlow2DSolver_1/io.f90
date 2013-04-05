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
module io
use datatypes
use constants
implicit none

contains

subroutine write_solution( fileT, fileQ,tR, geo, ti, op )
type(temp),intent(in) :: tR
type(geom),intent(in) :: geo
type(tim),intent(in) :: ti
type(opCond),intent(in) :: op
integer :: i
character(len=1000),intent(in) :: fileT,fileQ

!open (14, file='Tfout.dat',	&
!			           status='unknown', form='unformatted',	&
!			           access='direct', recl=recl*ti%n_timesteps)
!write(14,rec=1) tR%Tfout
!close(14)

open (15, file=fileT,   &
			status='replace', access='sequential',	&
			form='formatted',&
			action='write' )
			do i=1,ti%n_timesteps
			    write(15,*) ti%t(i),tR%TfHout(i),tR%TfCout(i)
			enddo
close(15)

open (15, file=fileQ,   &
			status='replace', access='sequential',	&
			form='formatted',&
			action='write' )
			do i=1,op%maxIte
			    write(15,*) i,tR%qhot(i),tR%qcold(i)
			enddo
close(15)

end subroutine write_solution

subroutine getFileNames( file, fileIn, fileOutT,fileOutQ )
character(len=1000),intent(in) :: file
character(len=1000),intent(inout) :: fileIn,fileOutT,fileOutQ

open(11,file=file,status='old',access='sequential',form='formatted',action='read')

read(11,*) fileIn
read(11,*) fileOutT
read(11,*) fileOutQ

close(11)

end subroutine getFileNames

subroutine loadTinProfile ( prof, op )
type(profSetup),intent(inout) :: prof
type(opCond),intent(inout) :: op

integer :: i
if ( prof%TinProfile .eq. TinputFile ) then
    allocate( prof%Tin(2,prof%TinNt) )
    open(11,file=prof%TinFile,status='old',access='sequential',&
    form='formatted',action='read')

    do i=1,prof%TinNt
        read(11,*) prof%Tin(1,i),prof%Tin(2,i)
    enddo
    close(11)
    op%Tf0 = prof%Tin(2,1)
    op%Ts0 = prof%Tin(2,1)
    op%Tw0 = prof%Tin(2,1)
    write(*,*) 'minmax Tin prof',minval(prof%Tin(2,:)),maxval(prof%Tin(2,:))
    
endif


end subroutine loadTinProfile

subroutine loadFlowProfile( prof )
type(profSetup),intent(inout) :: prof
integer :: i

if ( prof%flowProfile .eq. flowInputFile ) then
    allocate(prof%mdotIn(2,prof%nT))

    open(11,file=prof%profileFile,status='old',access='sequential',&
    form='formatted',action='read')

    do i=1,prof%nT
        read(11,*) prof%mdotIn(1,i),prof%mdotIn(2,i)
    enddo
    close(11)
    prof%mdot0 = maxval(prof%mdotIn(2,:))
    
    write(*,*) 'minmax mdot prof',minval(prof%mdotIn(2,:)),maxval(prof%mdotIn(2,:))
    
endif
end subroutine loadFlowProfile

!::Load any user-provided properties into tables
!::The structure of the input files is:
!:: T1 k[]1
!:: T2 k[]2
!:: .. ..
!:: Tn k[]n
!:: T1 rho[]1
!:: T2 rho[]2
!:: .. ..
!:: Tn rho[]n
!:: T1 c[]1
!:: T2 c[]2
!:: .. ..
!:: Tn c[]n
!:: (and for the fluid)
!:: T1 muf1
!:: T2 muf2
!:: .. ..
!:: Tn mufn
subroutine loadProps( props )
type(propSetup),intent(inout) :: props
integer :: i
if ( props%fluidType .eq. domTypeUserDef ) then
    !::Load the fluid properties
    allocate( props%kf(2,props%kfNT) )
    allocate( props%rhof(2,props%rhofNT) )
    allocate( props%cf(2,props%cfNT) )
    allocate( props%muf(2,props%mufNT) )
    open(11,file=props%fluidFile,status='old',access='sequential',&
    form='formatted',action='read')
    do i=1,props%kfNT
        read(11,*) props%kf(1,i),props%kf(2,i)
    enddo
    do i=1,props%rhofNT
        read(11,*) props%rhof(1,i),props%rhof(2,i)
    enddo
    do i=1,props%cfNT
        read(11,*) props%cf(1,i),props%cf(2,i)
    enddo
    do i=1,props%mufNT
        read(11,*) props%muf(1,i),props%muf(2,i)
    enddo
    close(11)
endif

if ( props%solidType .eq. domTypeUserDef ) then
    !::Load the solid properties
    allocate( props%ks(2,props%ksNT) )
    allocate( props%rhos(2,props%rhosNT) )
    allocate( props%cs(2,props%csNT) )
        
    open(11,file=props%solidFile,status='old',access='sequential',&
    form='formatted',action='read')
    do i=1,props%ksNT
        read(11,*) props%ks(1,i),props%ks(2,i)
    enddo
    do i=1,props%rhosNT
        read(11,*) props%rhos(1,i),props%rhos(2,i)
    enddo
    do i=1,props%csNT
        read(11,*) props%cs(1,i),props%cs(2,i)
    enddo    
    close(11)
endif

if ( props%wallType .eq. domTypeUserDef ) then
    !::Load the solid properties
    allocate( props%kw(2,props%kwNT) )
    allocate( props%rhow(2,props%rhowNT) )
    allocate( props%cw(2,props%cwNT) )
        
    open(11,file=props%wallFile,status='old',access='sequential',&
    form='formatted',action='read')
    do i=1,props%kwNT
        read(11,*) props%kw(1,i),props%kw(2,i)
    enddo
    do i=1,props%rhowNT
        read(11,*) props%rhow(1,i),props%rhow(2,i)
    enddo
    do i=1,props%cwNT
        read(11,*) props%cw(1,i),props%cw(2,i)
    enddo    
    close(11)

endif

end subroutine loadProps


!::Read the input from file
subroutine read_input( file, geo, op, ti, props,prof )
character(len=1000),intent(in) :: file
type(geom),intent(inout) :: geo
type(opCond),intent(inout) :: op
type(tim),intent(inout) :: ti
type(propSetup),intent(inout) :: props
type(profSetup),intent(inout) :: prof

open(11,file=file,status='old',access='sequential',form='formatted',action='read')

!::Read in the geometric parameters
!::Grid resolution
read(11,*) geo%nx,geo%nr_sf,geo%nr_w
write(*,*) 'nx,nr_sf,nr_w'
write(*,*) geo%nx,geo%nr_sf,geo%nr_w
!::Geometry
read(11,*) geo%L,geo%Rad,geo%Hw
write(*,*) 'L,Rad,Wt'
write(*,*) geo%L,geo%Rad,geo%Hw

!::Geometry type
read(11,*) geo%geomType
write(*,*) 'GeomType'
write(*,*) geo%geomType
!::Based on the geometry type load the relevant data
if ( geo%geomType .eq. geomTypeSpheres ) then
    read(11,*) geo%por,geo%dpar
    write(*,*) 'por,dpar'
    write(*,*) geo%por,geo%dpar
endif

!::Read in the operating conditions
read(11,*) prof%flowProfile
write(*,*) 'Flow profile (1=use the box profile,2=use input file)'
write(*,*) prof%flowProfile

if ( prof%flowProfile .eq. flowBoxProfile ) then
    read(11,*) prof%mdot0
    write(*,*) 'mdot0'
    write(*,*) prof%mdot0
elseif ( prof%flowProfile .eq. flowInputFile ) then
    read(11,*) prof%profileFile,prof%nT
    write(*,*) 'Flow profile file,nT'
    write(*,*) prof%profileFile,prof%nT
endif
read(11,*) op%axial,op%walls,op%wallsSol,op%dispAx,op%dispRad,op%kfBdry
write(*,*) 'axial,walls,wallsSol,dispAx,dispRad,kfBdry'
write(*,*) op%axial,op%walls,op%wallsSol,op%dispAx,op%dispRad,op%kfBdry

read(11,*) op%Nuscl,op%viscDiss
write(*,*) 'Nuscl,Viscous dissipation'
write(*,*) op%Nuscl,op%viscDiss

read(11,*) prof%TinProfile
write(*,*) 'Tin profile (1=use a constant,2=use input file)'
write(*,*) prof%TinProfile

if ( prof%TinProfile .eq. TinConstant ) then
    read(11,*) op%Tc,op%Th,op%Tf0,op%Ts0,op%Tw0
    write(*,*) 'Tc,Th,Tf0,Ts0,Tw0'
    write(*,*) op%Tc,op%Th,op%Tf0,op%Ts0,op%Tw0   
elseif ( prof%TinProfile .eq. TinputFile ) then
    read(11,*) prof%TinFile,prof%TinNt
    write(*,*) 'Tin file, Nt (Tin)'
    write(*,*) prof%TinFile, prof%TinNt    
endif


read(11,*) op%hw
write(*,*) 'Hw'
write(*,*) op%hw

!::Read in the timing variables
read(11,*) ti%t_end,ti%CFL_max
write(*,*) 't_end,CFL_max'
write(*,*) ti%t_end,ti%CFL_max

!::Properties (either use build-in fct or user provided data)
read(11,*) props%fluidType,props%solidType,props%wallType
write(*,*) 'fluidtype,solidtype,walltype'
write(*,*) props%fluidType,props%solidType,props%wallType

if ( props%fluidType .eq. domTypeUserDef ) then
    !::User specified fluid properties; ask for data file
    read(11,*) props%fluidFile,props%kfNT,props%rhofNT,props%cfNT,props%mufNT
    write(*,*) 'fluid file,kfNT,rhofNT,cfNT,mufNT'
    write(*,*) props%fluidFile,props%kfNT,props%rhofNT,props%cfNT,props%mufNT    
endif

if ( props%solidType .eq. domTypeUserDef ) then
    !::User specified solid properties; ask for data file
    read(11,*) props%solidFile,props%ksNT,props%rhosNT,props%csNT
    write(*,*) 'solid file,ksNT,rhosNT,csNT'
    write(*,*) props%solidFile,props%ksNT,props%rhosNT,props%csNT
endif

if ( props%wallType .eq. domTypeUserDef ) then
    !::User specified wall properties; ask for data file
    read(11,*) props%wallFile,props%kwNT,props%rhowNT,props%cwNT
    write(*,*) 'wall file,kwNT,rhowNT,cwNT'
    write(*,*) props%wallFile,props%kwNT,props%rhowNT,props%cwNT
endif

read(11,*) op%blowMode
write(*,*) 'SingleBlow or transient steady-state (1 if single,2 if transient)'
write(*,*) op%blowMode

if ( op%blowMode .eq. PeriodicMode ) then
    read(11,*) op%convTol,op%maxIte
    write(*,*) 'Convergence tolerance, max. nr. iterations'
    write(*,*) op%convTol,op%maxIte
else
    op%maxIte = 1
endif

close(11)

end subroutine read_input

end module io