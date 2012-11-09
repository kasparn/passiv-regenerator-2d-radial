MODULE SOLVER_CALL
use matrix_inversion
!use solver_f90_test
implicit none
contains


!::
!::Progresses the unsteady heat transfer equations one timestep, dt
!::
subroutine timestep( dx, dr, r, rup, rdn, nx, nr_sf, nr_w,&
                        rhof,kf,cf,rhos,ks,cs,rhow,kw,cw,h,as,&
                        dt, por, Tf, Ts, Tw, u, hw, Qvisc, Tc, Th, &
                        axial, walls, Tfout, Tsout,Twout, A, B, Ainv, X)
real,intent(in) :: as,dt,por,axial,walls,Tc,Th
real,intent(in),dimension(nx,nr_sf) :: rhof,kf,cf,rhos,ks,cs,h,u,Qvisc
real,intent(in),dimension(nx,nr_sf) :: Tf,Ts
real,intent(inout),dimension(nx,nr_sf) :: Tfout,Tsout
real,intent(in),dimension(nx,nr_w) :: rhow,kw,cw
real,intent(in),dimension(nx,nr_w) :: Tw
real,intent(inout),dimension(nx,nr_w) :: Twout
real,intent(in),dimension(nx) :: hw,dx
real,intent(in),dimension(nr_sf+nr_w) :: dr,r,rup,rdn
real,dimension(nx*(2*nr_sf+nr_w),nx*(2*nr_sf+nr_w)) :: A,Ainv
real,dimension(nx*(2*nr_sf+nr_w)) :: B,X
integer,intent(in) :: nx,nr_sf,nr_w
integer :: i,j,rw,cl,jj,ipiv,info
real :: t1,t2

                    
call cpu_time( t1 )


!allocate( A( nx * ( 2 * nr_sf + nr_w ), nx * ( 2 * nr_sf + nr_w ) ) )
!allocate( B( nx * ( nr_sf * 2 + nr_w )) )
A(:,:) = 0
Ainv(:,:) = 0
B(:) = 0
X(:) = 0 

!loop through each node
do i=2,nx-1
   !The solid and fluid nodes

   !Lower line, i.e. j = 1
   rw = i 
   cl = i 
   !fluid part, central node            
   A( rw, cl ) = por * rhof(i,1) * cf(i,1) / dt + &!capacity term
                 2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/kf(i,1) + dr(1+1)/kf(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/kf(i,1) + dx(i+1)/kf(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kf(i,1) + dx(i-1)/kf(i-1,1) ) ) + &!conduction to left node
                 h(i,1) * as !heat transfer with the solid

   !fluid part, right node
   A( rw, cl+1 ) = axial * -2 / ( dx(i) * ( dx(i)/kf(i,1) + dx(i+1)/kf(i+1,1) ) ) !conduction to right node

   !fluid part, left node
   A( rw, cl-1 ) = axial * -2 / ( dx(i) * ( dx(i)/kf(i,1) + dx(i-1)/kf(i-1,1) ) ) !conduction to left node

   !fluid part, upper node
   A( rw, cl + nx ) = -2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/kf(i,1) + dr(1+1)/kf(i,1+1) ) ) !conduction to node r+dr

   !heat transfer between solid and fluid
   A( rw, cl + nx*nr_sf ) = - h(i,1) * as 

   !convection   
   A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(i,1) * cf(i,1) * u(i,1) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
   A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(i,1) * cf(i,1) * u(i,1) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 

   rw = i + nx * nr_sf 
   cl = i + nx * nr_sf 
   !solid part, central node
   A( rw, cl ) = (1 - por) * rhos(i,1) * cs(i,1) / dt + &!capacity term
                 2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/ks(i,1) + dr(1+1)/ks(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i+1)/ks(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i-1)/ks(i-1,1) ) ) + &!conduction to left node
                 h(i,1) * as !heat transfer between fluid and solid


   !solid part, right node
   A( rw, cl + 1 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i+1)/ks(i+1,1) ) ) !conduction to right node

   !solid part, left node
   A( rw, cl - 1 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i-1)/ks(i-1,1) ) ) !conduction to left node

   !solid part, upper node
   A( rw, cl + nx ) = -2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/ks(i,1) + dr(1+1)/ks(i,1+1) ) ) !conduction to node r+dr

   !HT with fluid
   A( rw, cl - nx * nr_sf ) = -h(i,1) * as !heat transfer between fluid and solid


   !lower line, i=1
   rw = i + 2 * nx * nr_sf 
   cl = i + 2 * nx * nr_sf 
   jj = nr_sf 
   !wall part, central node
   A( rw, cl ) = rhow(i,1) * cw(i,1) / dt + &!capacity term
                 2 / ( r(jj+1) * dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1+1)/kw(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i+1)/kw(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i-1)/kw(i-1,1) ) ) + &!conduction to left node
                 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/kf(i,nr_sf) + 1/hw(i) ) ) !HT with the fluid

   !wall part, right node
   A( rw, cl + 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i+1)/kw(i+1,1) ) ) !conduction to right node

   !wall part, left node
   A( rw, cl - 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i-1)/kw(i-1,1) ) ) !conduction to left node

   !wall part, upper node
   A( rw, cl + nx ) = -2 / ( r(jj+1) * dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1+1)/kw(i,1+1) ) ) !conduction to node r+dr

   !HT with the fluid
   A( rw, cl - nx * ( nr_sf + 1 ) ) = walls * -2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/kf(i,nr_sf) + 1/hw(i) ) ) 

   !Upper line, i.e. j = nr_sf
   rw = i + nx * (nr_sf-1) 
   cl = i + nx * (nr_sf-1) 
   !fluid part, central node
   A( rw, cl ) = por * rhof(i,nr_sf) * cf(i,nr_sf) / dt + &!capacity term
                 2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(i,nr_sf) + dr(nr_sf-1)/kf(i,nr_sf-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/kf(i,nr_sf) + dx(i+1)/kf(i+1,nr_sf) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kf(i,nr_sf) + dx(i-1)/kf(i-1,nr_sf) ) ) + &!conduction to left node
                 h(i,nr_sf) * as + &!HT btw. solid and fluid
                 walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/kf(i,nr_sf) + 1/hw(i) ) ) !HT with the wall

   !fluid part, right node
   A( rw, cl + 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kf(i,nr_sf) + dx(i+1)/kf(i+1,nr_sf) ) ) !conduction to right node

   !fluid part, left node
   A( rw, cl - 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kf(i,nr_sf) + dx(i-1)/kf(i-1,nr_sf) ) ) !conduction to left node

   !fluid part, lower node
   A( rw, cl - nx ) = &
       -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(i,nr_sf) + dr(nr_sf-1)/kf(i,nr_sf-1) ) ) !conduction to node r-dr

   A( rw, cl + nx * nr_sf ) = -h(i,nr_sf) * as !HT btw. solid and fluid

   !fluid part, HT with the wall
   A( rw, cl + nx * ( nr_sf + 1 ) ) =&
       walls * -2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/kf(i,nr_sf) + 1/hw(i) ) ) !HT with the wall

   !convection
   
   A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(i,nr_sf) * cf(i,nr_sf) * u(i,nr_sf) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
   A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(i,nr_sf) * cf(i,nr_sf) * u(i,nr_sf) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 

   rw = i + nx * nr_sf + nx * ( nr_sf - 1 ) 
   cl = i + nx * nr_sf + nx * ( nr_sf - 1 ) 
   !solid part, central node
   A( rw, cl ) = (1 - por) * rhos(i,nr_sf) * cs(i,nr_sf) / dt + &!capacity term
                 2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(i,nr_sf) + dr(nr_sf-1)/ks(i,nr_sf-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i+1)/ks(i+1,nr_sf) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i-1)/ks(i-1,nr_sf) ) ) + &!conduction to left node
                 h(i,nr_sf) * as !HT btw. solid and fluid

   !solid part, right node
   A( rw, cl + 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i+1)/ks(i+1,nr_sf) ) ) !conduction to right node

   !solid part, left node
   A( rw, cl - 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i-1)/ks(i-1,nr_sf) ) ) !conduction to left node

   !solid part, lower node
   A( rw, cl - nx ) = &
       -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(i,nr_sf) + dr(nr_sf-1)/ks(i,nr_sf-1) ) ) !conduction to node r-dr

   A( rw, cl - nx*nr_sf ) = -h(i,nr_sf) * as !HT btw. solid and fluid


   rw = i + 2 * nx * nr_sf + nx * ( nr_w - 1 ) 
   cl = i + 2 * nx * nr_sf + nx * ( nr_w - 1 ) 
   jj = nr_sf 
   !wall part, central node
   A( rw, cl ) = rhow(i,nr_w) * cw(i,nr_w) / dt + &!capacity term
                 2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(i,nr_w) + dr(jj+nr_w-1)/kw(i,nr_w-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i+1)/kw(i+1,nr_w) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i-1)/kw(i-1,nr_w) ) ) !conduction to left node

   !wall part, right node
   A( rw, cl + 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i+1)/kw(i+1,nr_w) ) ) !conduction to right node

   !wall part, left node
   A( rw, cl - 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i-1)/kw(i-1,nr_w) ) ) !conduction to left node

   !wall part, lower node
   A( rw, cl - nx ) = &
       -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(i,nr_w) + dr(jj+nr_w-1)/kw(i,nr_w-1) ) ) !conduction to node r-dr

   !inner nodes
   do j=2,nr_sf-1
      rw = nx*(j-1) + i 
      cl = nx*(j-1) + i 
      !fluid part, central nodes
      A( rw, cl ) =  por * rhof(i,j) * cf(i,j) / dt + &!capacity term
                     2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(i,j) + dr(j+1)/kf(i,j+1) )) + &!conduction to node r+dr
                     2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(i,j) + dr(j-1)/kf(i,j-1) ) ) + &!conduction to node r-dr
                     axial * 2 / ( dx(i) * ( dx(i)/kf(i,j) + dx(i+1)/kf(i+1,j) ) ) + &!conduction to right node
                     axial * 2 / ( dx(i) * ( dx(i)/kf(i,j) + dx(i-1)/kf(i-1,j) ) ) + &!conduction to left node
                     h(i,j) * as !HT btw. solid and fluid

     !fluid part, upper
     A( rw, cl + nx ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(i,j) + dr(j+1)/kf(i,j+1) )) !conduction to node r+dr
     !fluid part, lower
     A( rw, cl - nx ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(i,j) + dr(j-1)/kf(i,j-1) ) ) !conduction to node r-dr

     !fluid part, right
     A( rw, cl + 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kf(i,j) + dx(i+1)/kf(i+1,j) ) ) !conduction to right node

     !fluid part, left
     A( rw, cl - 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kf(i,j) + dx(i-1)/kf(i-1,j) ) ) !conduction to left node

     A( rw, cl + nx * nr_sf ) = -h(i,j) * as !HT btw. solid and fluid
     !convection
     A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(i,j) * cf(i,j) * u(i,j) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
     A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(i,j) * cf(i,j) * u(i,j) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
         
     rw = nx * ( nr_sf + (j-1) ) + i 
     cl = nx * ( nr_sf + (j-1) ) + i 
     ! solid part, inner nodes
     A( rw, cl ) =  (1 - por) * rhos(i,j) * cs(i,j) / dt + &!capacity term
                     2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) + &!conduction to node r+dr
                     2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) + &!conduction to node r-dr
                     axial * 2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i+1)/ks(i+1,j) ) ) + &!conduction to right node
                     axial * 2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i-1)/ks(i-1,j) ) ) + &!conduction to left node
                     h(i,j) * as !HT btw. solid and fluid

     !solid part, upper
     A( rw, cl + nx ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) )) !conduction to node r+dr
     !solid part, lower
     A( rw, cl - nx ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) !conduction to node r-dr

     !solid part, right
     A( rw, cl + 1 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i+1)/ks(i+1,j) ) ) !conduction to right node

     !solid part, left
     A( rw, cl - 1 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i-1)/ks(i-1,j) ) ) !conduction to left node

     A( rw, cl - nx*nr_sf ) = -h(i,j) * as !HT btw. solid and fluid


   enddo
   !wall, inner nodes (need their own loop since nr_sf ~= nr_w in
   !general
   do j=2,nr_w-1
     rw = 2 * nx * nr_sf + nx * (j-1) + i 
     cl = 2 * nx * nr_sf + nx * (j-1) + i 
     jj = nr_sf 
     ! wall part, inner nodes
     A( rw, cl ) =   rhow(i,j) * cw(i,j) / dt + &!capacity term
                     2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) + &!conduction to node r+dr
                     2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) + &!conduction to node r-dr
                     axial * 2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i+1)/kw(i+1,j) ) ) + &!conduction to right node
                     axial * 2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i-1)/kw(i-1,j) ) ) !conduction to left node                         

     !wall part, upper
     A( rw, cl + nx ) = -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) )) !conduction to node r+dr
     !wall part, lower
     A( rw, cl - nx ) = -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) !conduction to node r-dr

     !wall part, right
     A( rw, cl + 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i+1)/kw(i+1,j) ) ) !conduction to right node

     !wall part, left
     A( rw, cl - 1 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i-1)/kw(i-1,j) ) ) !conduction to left node

   enddo
enddo

!the two side-lines
do j=2,nr_sf-1
    !side line one (i = 1)
    rw = nx * ( j - 1 ) + 1 
    cl = nx * ( j - 1 ) + 1         
    !fluid part, central node
    A( rw, cl ) = por * rhof(1,j) * cf(1,j) / dt + &!capacity term
        2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(i,j) + dr(j+1)/kf(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(i,j) + dr(j-1)/kf(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/kf(1,j) + dx(1+1)/kf(1+1,j) ) ) + &!conduction to right node
        h(1,j) * as !HT btw. solid and fluid

    !fluid part, upper
    A( rw, cl + nx ) = &
        -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(i,j) + dr(j+1)/kf(i,j+1) ) ) !conduction to node r+dr

    !fluid part, lower
    A( rw, cl - nx ) = &
        -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(i,j) + dr(j-1)/kf(i,j-1) ) ) !conduction to node r-dr

    !fluid part, right node
    A( rw, cl + 1 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/kf(1,j) + dx(1+1)/kf(1+1,j) ) ) !conduction to right node

    A( rw, cl + nx * nr_sf ) = -h(1,j) * as !HT btw. solid and fluid
    !convection
    if ( u(1,j) .gt. 0 ) then
       !there will then be an explicit term in the B vector
       !A( rw, cl ) = A( rw, cl ) + rhof(1,j) * cf(1,j) * u(1,j) * por / dx(1) 
       A( rw, cl ) = A( rw, cl ) + rhof(1,j) * cf(1,j) * u(1,j) * por / ( 2 * dx(1) ) 
       
       A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,j) * cf(1,j) * u(1,j) * por / ( 2 * dx(1) ) 
       
    else if ( u(1,j) .lt. 0 ) then
       !outlet, use up-wind in this case
       A( rw, cl ) = A( rw, cl ) - rhof(1,j) * cf(1,j) * u(1,j) * por / dx(1) 
       A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,j) * cf(1,j) * u(1,j) * por / dx(1) 
    endif
 

    rw = nx * ( j - 1 ) + 1 + nx*nr_sf 
    cl = nx * ( j - 1 ) + 1 + nx*nr_sf 
    !solid part, central node
    A( rw, cl ) = (1 - por) * rhos(1,j) * cs(1,j) / dt + &!capacity term
        2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/ks(1,j) + dx(1+1)/ks(1+1,j) ) ) + &!conduction to right node
        h(1,j) * as 

    !solid part, upper
    A( rw, cl + nx ) = &
        -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) !conduction to node r+dr

    !solid part, lower
    A( rw, cl - nx ) = &
        -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) !conduction to node r-dr

    !solid part, right node
    A( rw, cl + 1 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/ks(1,j) + dx(1+1)/ks(1+1,j) ) ) !conduction to right node

    A( rw, cl - nx * nr_sf ) = -h(i,1) * as !HT btw. solid and fluid

    rw = nx * j 
    cl = nx * j 
    !sideline 2, i=nx
    !fluid part, central node
    A( rw, cl ) = por * rhof( nx, j ) * cf( nx, j ) / dt + &!capacity term
                  2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(nx,j) + dr(j+1)/kf(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(nx,j) + dr(j-1)/kf(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/kf(nx,j) + dx(nx-1)/kf(nx-1,j) ) ) + &!conduction to left node
                  h(nx,j) * as 


    !fluid part, upper
    A( rw, cl + nx ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kf(nx,j) + dr(j+1)/kf(nx,j+1) )) !conduction to node r+dr

    !fluid part, lower
    A( rw, cl - nx ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kf(nx,j) + dr(j-1)/kf(nx,j-1) )) !conduction to node r-dr

    !fluid part, left
    A( rw, cl - 1) = axial * -2 / ( dx(nx) * ( dx(nx)/kf(nx,j) + dx(nx-1)/kf(nx-1,j) ) ) !conduction to left node

    A( rw, cl + nx * nr_sf ) = -h(nx,j) * as !HT btw. solid and fluid
    !convection
    
    if ( u(nx,j) .gt. 0 ) then
        !outlet, use up-wind scheme
        !central node
        A( rw, cl ) = A( rw, cl ) + rhof(nx,j) * cf(nx,j) * u(nx,j) * por / dx(nx) 
        !left node
        A( rw, cl - 1) = A( rw, cl - 1 ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( dx(nx) ) 
    elseif ( u(nx,j) .lt. 0 ) then
        !inlet (then an explicit term is added to the B vector
        !central node
        A( rw, cl ) = A( rw, cl ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( 2 * dx(nx)  ) 
        
        A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( 2 * dx(nx)  ) 
        
    endif

    rw = nx * j + nx*nr_sf 
    cl = nx * j + nx*nr_sf 
    !solid part, central node
    A( rw, cl ) = ( 1 - por ) * rhos( nx, j ) * cs( nx, j ) / dt + &!capacity term
                  2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(nx,j) + dr(j+1)/ks(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(nx,j) + dr(j-1)/ks(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,j) + dx(nx-1)/ks(nx-1,j) ) ) + &!conduction to left node
                  h(nx, j) * as 

    !solid part, upper
    A( rw, cl + nx ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(nx,j) + dr(j+1)/ks(nx,j+1) )) !conduction to node r+dr

    !solid part, lower
    A( rw, cl - nx ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(nx,j) + dr(j-1)/ks(nx,j-1) )) !conduction to node r-dr

    !solid part, left
    A( rw, cl - 1) = axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,j) + dx(nx-1)/ks(nx-1,j) ) ) !conduction to left node

    A( rw, cl - nx * nr_sf ) = -h(nx,j) * as !HT btw. solid and fluid
enddo
!the two side lines for the wall (nr_w ~= nr_sf in general)
do j=2,nr_w-1
    !left most line, i=1
    rw = 2 * nx * nr_sf + nx * ( j - 1 ) + 1 
    cl = 2 * nx * nr_sf + nx * ( j - 1 ) + 1 
    jj = nr_sf 
    !wall part, central node
    A( rw, cl ) = rhow(1,j) * cw(1,j) / dt + &!capacity term
        2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/kw(1,j) + dx(1+1)/kw(1+1,j) ) ) !conduction to right node            

    !wall part, upper
    A( rw, cl + nx ) = &
        -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) !conduction to node r+dr

    !wall part, lower
    A( rw, cl - nx ) = &
        -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) !conduction to node r-dr

    !wall part, right node
    A( rw, cl + 1 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/kw(1,j) + dx(1+1)/kw(1+1,j) ) ) !conduction to right node

    !side line to the right (i=nx)
    rw = 2 * nx * nr_sf + nx * j 
    cl = 2 * nx * nr_sf + nx * j 
    !wall part, central node
    A( rw, cl ) = rhow( nx, j ) * cw( nx, j ) / dt + &!capacity term
                  2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j+1)/kw(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j-1)/kw(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,j) + dx(nx-1)/kw(nx-1,j) ) ) !conduction to left node                      

    !wall part, upper
    A( rw, cl + nx ) = -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j+1)/kw(nx,j+1) )) !conduction to node r+dr

    !wall part, lower
    A( rw, cl - nx ) = -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j-1)/kw(nx,j-1) )) !conduction to node r-dr

    !wall part, left
    A( rw, cl - 1) = axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,j) + dx(nx-1)/kw(nx-1,j) ) ) !conduction to left node

enddo

!!The lower left corner
rw = 1 
cl = 1 
!fluid part, central node
A( rw, cl ) =  por * rhof(1,1) * cf(1,1) / dt + &!capacity term
             2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kf(1,1) + dr(1+1)/kf(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/kf(1,1) + dx(1+1)/kf(1+1,1) ) ) + &!conduction to right node   
             h(1,1) * as 

!fluid part, upper
A( rw, cl + nx ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kf(1,1) + dr(1+1)/kf(1,1+1) )) !conduction to node r+dr

!fluid part, right
A( rw, cl + 1 ) = axial * -2 / ( dx(1) * ( dx(1)/kf(1,1) + dx(1+1)/kf(1+1,1) ) ) !conduction to right node

A( rw, cl + nx * nr_sf ) = -h(1,1) * as 
!convection
if ( u(1,1) .gt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( rw, cl ) = A( rw, cl ) + rhof(1,1) * cf(1,1) * u(1,1) * por / (2 * dx(1)) 
    
    A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,1) * cf(1,1) * u(1,1) * por / (2 * dx(1)) 
    !A( rw, cl ) = A( rw, cl ) + rhof(1,1) * cf(1,1) * u(1,1) * por / dx(1) 
elseif ( u(1,1) .lt. 0 ) then
    !outlet, use the up-wind scheme
    !central node
    A( rw, cl ) = A( rw, cl ) - rhof(1,1) * cf(1,1) * u(1,1) * por / dx(1) 
    !right node
    A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,1) * cf(1,1) * u(1,1) / dx(1) 
endif

rw = 1 + nx * nr_sf 
cl = 1 + nx * nr_sf 
!solid part, central node
A( rw, cl ) =  ( 1 - por ) * rhos(1,1) * cs(1,1) / dt + &!capacity term
             2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(1,1) + dr(1+1)/ks(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/ks(1,1) + dx(1+1)/ks(1+1,1) ) ) + &!conduction to right node                
             h(1,1) * as 

!solid part, upper
A( rw, cl + nx ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(1,1) + dr(1+1)/ks(1,1+1) )) !conduction to node r+dr

!solid part, right
A( rw, cl + 1 ) = axial * -2 / ( dx(1) * ( dx(1)/ks(1,1) + dx(1+1)/ks(1+1,1) ) ) !conduction to right node

A( rw, cl - nx * nr_sf ) = -h(1,1) * as 

!wall
rw = 1 + 2 * nx * nr_sf 
cl = 1 + 2 * nx * nr_sf 
jj = nr_sf 
!wall part, central node
A( rw, cl ) = rhow(1,1) * cw(1,1) / dt + &!capacity term
             2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(1,1) + dr(jj+1+1)/kw(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/kw(1,1) + dx(1+1)/kw(1+1,1) ) ) + &!conduction to right node                 
             walls * 2 / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/kf(1,nr_sf) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) !HT with the fluid

!wall part, upper
A( rw, cl + nx ) = -2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(1,1) + dr(jj+1+1)/kw(1,1+1) )) !conduction to node r+dr

!wall part, right
A( rw, cl + 1 ) = axial * -2 / ( dx(1) * ( dx(1)/kw(1,1) + dx(1+1)/kw(1+1,1) ) ) !conduction to right node

!HT with the fluid
A( rw, cl - nx * ( nr_sf + 1 ) ) = walls * -2 / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/kf(1,nr_sf) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) 

rw = nx 
cl = nx 
!!Lower right corner
!fluid part, central node
A( rw, cl ) = por * rhof( nx, 1 ) * cf( nx, 1 ) / dt + &!capacity term
              2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kf(nx,1) + dr(1+1)/kf(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/kf(nx,1) + dx(nx-1)/kf(nx-1,1) ) ) + &!conduction to left node
              h(nx, 1 ) * as 

!fluid part, upper node
A( rw, cl + nx ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kf(nx,1) + dr(1+1)/kf(nx,1+1) )) !conduction to node r+dr

!fluid part, left node
A( rw, cl - 1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/kf(nx,1) + dx(nx-1)/kf(nx-1,1) ) ) !conduction to the left node

A( rw, nx * nr_sf + nx ) = -h(nx,1) * as 
if ( u(nx,1) .gt. 0 ) then
    !outlet, use the upwind scheme
    !central node
    A( rw, cl ) = A( rw, cl ) + rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( dx(nx) ) 
    !left node
    A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( dx(nx) ) 
elseif ( u(nx,1) .lt. 0 ) then
    !inlet, (then an explicit term is added to the B vector)
    !central node
    A( rw, cl ) = A( rw, cl ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( 2 * dx(nx) ) 
    
    A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( 2 * dx(nx) ) 
endif

rw = nx  + nx * nr_sf 
cl = nx  + nx * nr_sf 
!solid part, central node
A( rw, cl ) = ( 1 - por ) * rhos( nx, 1 ) * cs( nx, 1 ) / dt + &!capacity term
              2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(nx,1) + dr(1+1)/ks(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,1) + dx(nx-1)/ks(nx-1,1) ) ) + & !conduction to left node   
              h(nx,1) * as 

!solid part, upper node
A( rw, cl + nx ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(nx,1) + dr(1+1)/ks(nx,1+1) )) !conduction to node r+dr

!solid part, left node
A( rw, cl - 1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,1) + dx(nx-1)/ks(nx-1,1) ) ) !conduction to the left node
!HT with fluid
A( rw, cl - nx * nr_sf ) = -h(nx,1) * as 


!wall part, lower right corner
rw = nx  + 2 * nx * nr_sf 
cl = nx  + 2 * nx * nr_sf 
jj = nr_sf 
!wall part, central node
A( rw, cl ) = rhow( nx, 1 ) * cw( nx, 1 ) / dt + &!capacity term
              2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(nx,1) + dr(jj+1+1)/kw(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,1) + dx(nx-1)/kw(nx-1,1) ) ) + &!conduction to left node
              walls * 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/kf(nx,nr_sf) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) !HT with the fluid

!wall part, upper node
A( rw, cl + nx ) = -2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(nx,1) + dr(jj+1+1)/kw(nx,1+1) )) !conduction to node r+dr

!wall part, left node
A( rw, cl - 1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,1) + dx(nx-1)/kw(nx-1,1) ) ) !conduction to the left node

!HT with the fluid
A( rw, cl - nx * ( nr_sf + 1 ) ) = walls * -2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/kf(nx,nr_sf) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) 

rw = nx * ( nr_sf - 1 ) + 1 
cl = nx * ( nr_sf - 1 ) + 1 
!!Upper left corner
!fluid part, central node
A( rw, cl ) = por * rhof( 1, nr_sf ) * cf( 1, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(1,nr_sf) + dr(nr_sf-1)/kf(1,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/kf(1,nr_sf) + dx(1+1)/kf(1+1,nr_sf) ) ) + &!conduction to right node
                h(1, nr_sf) * as + &!HT with the solid
                walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/kf(1,nr_sf) + 1/hw(1) ) ) !HT with the wall

!fluid part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(1,nr_sf) + dr(nr_sf-1)/kf(1,nr_sf-1) ) ) !conduction to node r-dr

!fluid part, right node
A( rw, cl + 1 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/kf(1,nr_sf) + dx(1+1)/kf(1+1,nr_sf) ) ) !conduction to right node

!fluid part, HT with the wall
A( rw, cl + nx * ( nr_sf + 1 ) ) = walls * -2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/kf(1,nr_sf) + 1/hw(1) ) ) 

!HT with the solid
A( rw, cl + nx * nr_sf ) = -h(1, nr_sf ) * as 

!convection 
 
if ( u(1,nr_sf) .gt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( rw, cl ) = A( rw, cl ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / ( 2*dx(1) ) 
    A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / ( 2*dx(1) ) 
    !A( rw, cl ) = A( rw, cl ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / dx(1) 

elseif ( u(1,nr_sf) .lt. 0 ) then
    !outlet, use the upwind scheme
    !central node
    A( rw, cl ) = A( rw, cl ) - rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / dx(1) 
    !right node
    A( rw, cl + 1 ) = A( rw, cl + 1 ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / dx(1) 
endif



rw = nx * ( nr_sf - 1 ) + 1 + nx * nr_sf 
cl = nx * ( nr_sf - 1 ) + 1 + nx * nr_sf 
!solid part, central node
A( rw, cl ) = (1 - por) * rhos( 1, nr_sf ) * cs( 1, nr_sf ) / dt + &!capacity term                    
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(1,nr_sf) + dr(nr_sf-1)/ks(1,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/ks(1,nr_sf) + dx(1+1)/ks(1+1,nr_sf) ) ) + &!conduction to right node
                h(1, nr_sf ) * as 

!solid part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(1,nr_sf) + dr(nr_sf-1)/ks(1,nr_sf-1) ) ) !conduction to node r-dr

!solid part, right node
A( rw, cl + 1 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/ks(1,nr_sf) + dx(1+1)/ks(1+1,nr_sf) ) ) !conduction to right node

A( rw, cl - nx * nr_sf ) = -h(1, nr_sf ) * as 

!upper left corner, wall
rw = nx * ( nr_w - 1 ) + 1 + 2 * nx * nr_sf 
cl = nx * ( nr_w - 1 ) + 1 + 2 * nx * nr_sf 
jj = nr_sf 
!wall part, central node
A( rw, cl ) =   rhow( 1, nr_w ) * cw( 1, nr_w ) / dt + &!capacity term                    
                2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(1,nr_w) + dr(jj+nr_w-1)/kw(1,nr_w-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/kw(1,nr_w) + dx(1+1)/kw(1+1,nr_w) ) ) !conduction to right node                    

!wall part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(1,nr_w) + dr(jj+nr_w-1)/kw(1,nr_w-1) ) ) !conduction to node r-dr

!wall part, right node
A( rw, cl + 1 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/kw(1,nr_w) + dx(1+1)/kw(1+1,nr_w) ) ) !conduction to right node


rw = nx * nr_sf 
cl = nx * nr_sf 
!!Upper right corner
!fluid part, central node
A( rw, cl ) = por * rhof( nx, nr_sf ) * cf( nx, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(nx,nr_sf) + dr(nr_sf-1)/kf(nx,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/kf(nx,nr_sf) + dx(nx-1)/kf(nx-1,nr_sf) ) ) + &!conduction to left node
                h( nx, nr_sf ) * as + &!HT with the solid
                walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/kf(nx,nr_sf) + 1/hw(nx) ) ) !HT with the wall

!fluid part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kf(nx,nr_sf) + dr(nr_sf-1)/kf(nx,nr_sf-1) ) ) !conduction to node r-dr

!fluid part, left node
A( rw, cl - 1 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/kf(nx,nr_sf) + dx(nx-1)/kf(nx-1,nr_sf) ) ) !conduction to left node

!HT with the solid
A( rw, cl + nx * nr_sf ) = -h( nx, nr_sf ) * as 

!HT with the wall
A( rw, cl + nx * ( nr_sf + 1 ) ) = walls * -2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/kf(nx,nr_sf) + 1/hw(nx) ) ) 
!convection 
if ( u(nx,nr_sf) .gt. 0 ) then
    !outlet, use the upwind scheme
    !central node
    A( rw, cl ) = A( rw, cl ) + rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / dx(nx) 
    !left node
    A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / dx(nx) 
elseif ( u(nx,nr_sf) .lt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( rw, cl ) = A( rw, cl ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / ( 2 * dx(nx) ) 
    A( rw, cl - 1 ) = A( rw, cl - 1 ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx-1,nr_sf) * por / ( 2 * dx(nx)) 
endif

rw = nx * nr_sf + nx * nr_sf 
cl = nx * nr_sf + nx * nr_sf 
!solid part, central node
A( rw, cl ) = (1 - por) * rhos( nx, nr_sf ) * cs( nx, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(nx,nr_sf) + dr(nr_sf-1)/ks(nx,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,nr_sf) + dx(nx-1)/ks(nx-1,nr_sf) ) ) + &!conduction to left node
                h( nx, nr_sf ) * as 

!solid part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(nx,nr_sf) + dr(nr_sf-1)/ks(nx,nr_sf-1) ) ) !conduction to node r-dr

!solid part, left node
A( rw, cl - 1 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,nr_sf) + dx(nx-1)/ks(nx-1,nr_sf) ) ) !conduction to left node

A( rw, cl - nx * nr_sf ) = -h( nx, nr_sf ) * as 

!upper right corner, wall
rw = nx * nr_w + 2 * nx * nr_sf 
cl = nx * nr_w + 2 * nx * nr_sf 
jj = nr_sf 
!wall part, central node
A( rw, cl ) =   rhow( nx, nr_w ) * cw( nx, nr_w ) / dt + &!capacity term
                2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(nx,nr_w) + dr(jj+nr_w-1)/kw(nx,nr_w-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,nr_w) + dx(nx-1)/kw(nx-1,nr_w) ) ) !conduction to left node                    

!wall part, lower node
A( rw, cl - nx ) = &
    -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(nx,nr_w) + dr(jj+nr_w-1)/kw(nx,nr_w-1) ) ) !conduction to node r-dr

!wall part, left node
A( rw, cl - 1 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,nr_w) + dx(nx-1)/kw(nx-1,nr_w) ) ) !conduction to left node


do i=1,nx
    do j=1,nr_sf
        !fluid nodes
        B( (j-1)*nx + i ) = por * rhof(i,j) * cf(i,j) * Tf(i,j) / dt + &!capacity term
                            Qvisc(i,j) !viscous dissipation term 
        !solid nodes
        B( (j-1)*nx + i + nx * nr_sf ) = ( 1 - por ) * rhos(i,j) * cs(i,j) * Ts(i,j) / dt 
    enddo
    do j=1,nr_w
        !wall nodes
        B( (j-1)*nx + i + 2 * nx * nr_sf ) = rhow(i,j) * cw(i,j) * Tw(i,j) / dt 
    enddo
enddo
!convection 
if ( u(1,1) .gt. 0 ) then
    !inlet from cold side (then an explicit term is added to the B vector)
    B( 1:nx:(nx*(nr_sf-1)+1) ) = B( 1:nx:(nx*(nr_sf-1)+1) ) + (rhof(1,:) * cf(1,:) * u(1,:) * por * Tc / dx(1))
endif

if ( u(nx,1) .lt. 0 ) then
    !inlet from hot side (then an explicit term is added to the B vector)
    B( nx:nx:(nx*nr_sf) ) = B( nx:nx:(nx*nr_sf) ) - (rhof(nx,:) * cf(nx,:) * u(nx,:) * por * Th / dx(nx))
endif
call cpu_time( t2 )
write(*,*) 'Bef matrix inversion'
write(*,*) t2-t1
!::Solve Ax = B
!call invMatrixLUD( A, Ainv, B, X, nx * ( 2*nr_sf + nr_w ) )
!call solve_lin_direct()
write(*,*) 'ipiv',ipiv
write(*,*) 'info',info
stop

do j=1,nr_sf
    Tfout(:,j) = X( nx*(j-1)+1:nx*j ) 
    Tsout(:,j) = X( (nx*nr_sf+nx*(j-1)+1):(nx*nr_sf+nx*j) ) 
enddo
do j=1,nr_w
    Twout(:,j) = X( (2*nx*nr_sf+nx*(j-1)+1):(2*nx*nr_sf+nx*j) ) 
enddo
write(*,*) Tfout
stop
end subroutine timestep

end module SOLVER_CALL