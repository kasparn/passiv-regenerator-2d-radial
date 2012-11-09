module flowprofiles
use datatypes
use constants
use util
implicit none



contains

function getFlowProfile( dt, n_timesteps, prof, op,mdotOut )
real,intent(in) :: dt
integer,intent(in) :: n_timesteps
type(profSetup),intent(in) :: prof
type(opCond),intent(in) :: op
integer :: getFlowProfile,i
real,dimension(n_timesteps),intent(inout) :: mdotOut
real,dimension(n_timesteps) :: t


if ( op%blowMode .eq. SingleBlowMode ) then

    if ( prof%flowProfile .eq. flowBoxProfile ) then
        !Default use the constant mdot function
        getFlowProfile = flowBox( prof%mdot0, dt, n_timesteps, mdotOut )
    elseif ( prof%flowProfile .eq. flowInputFile ) then
       
        t(1) = 0
        do i=2,n_timesteps
            t(i) = t(i-1) + dt
        enddo
        !Interpolate in a table provided externally
        call interp1( prof%mdotIn(1,:), prof%mdotIn(2,:), t, prof%nT, n_timesteps, mdotOut )
    endif
    !No error
    getFlowProfile = 0
    return
elseif ( op%blowMode .eq. PeriodicMode ) then
    !Box profile is the only working implementation currently
    getFlowProfile = flowBox( prof%mdot0, dt, n_timesteps, mdotOut )
    !::Make half the cycle have opposite sign
    
    mdotOut((n_timesteps/2+1):n_timesteps) = -1.*mdotOut((n_timesteps/2+1):n_timesteps) 
    !mdotOut = mdotOut * -1
    !No error
    getFlowProfile = 0
    return
endif

!Error
getFlowProfile = 1

    
end function getFlowProfile

!simply returns a constant mdot as a function of time
function flowBox( mdot0, dt, n_timesteps, mdotOut )
integer :: flowBox
integer,intent(in) :: n_timesteps
real,intent(in) :: dt,mdot0
real,dimension(n_timesteps),intent(inout) :: mdotOut

mdotOut = mdot0
!no error
flowBox = 0

end function flowBox

end module flowprofiles