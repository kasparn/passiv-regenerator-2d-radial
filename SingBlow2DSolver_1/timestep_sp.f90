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


MODULE SOLVER_CALL_SPARSE
use direct_solver
implicit none
contains


!::
!::Progresses the unsteady heat transfer equations one timestep, dt
!::The matrix assembled is sparse
!::
subroutine timestep_sp( dx, dr, r, rup, rdn, nx, nr_sf, nr_w,&
                        rhof,kfAx,kfRad,kfBdry,cf,rhos,ks,cs,rhow,kw,cw,h,as,&
                        dt, por, Tf, Ts, Tw, u, hw, Qvisc, Tc, Th, &
                        axial, walls, wallsSol,Tfout, Tsout,Twout )
real,intent(in) :: as,dt,por,axial,walls,wallsSol,Tc,Th
real,intent(in),dimension(nx,nr_sf) :: rhof,kfAx,kfRad,cf,rhos,ks,cs,h,u,Qvisc
real,intent(in),dimension(nx,nr_sf) :: Tf,Ts
real,intent(inout),dimension(nx,nr_sf) :: Tfout,Tsout
real,intent(in),dimension(nx,nr_w) :: rhow,kw,cw
real,intent(in),dimension(nx,nr_w) :: Tw
real,intent(inout),dimension(nx,nr_w) :: Twout
real,intent(in),dimension(nx) :: hw,dx,kfBdry
real,intent(in),dimension(nr_sf+nr_w) :: dr,r,rup,rdn
real,dimension(:),allocatable :: A,B,X
integer,dimension(:),allocatable :: cols,rows
integer,intent(in) :: nx,nr_sf,nr_w
integer :: i,j,jj,nNonZero,ind,indFluidOffSet,indFluidSolidOffSet
integer :: rw_in
real :: t1,t2,t3
integer,dimension(:),allocatable :: rowIndex,columns
real,dimension(:),allocatable :: values,rhs,sol

                    
!call cpu_time( t1 )

!::Find the number of non-zero elements
nNonZero = 2 * ( (nx-2)*(nr_sf-2)*6 + 2*(nx-2)*5 + 2*(nr_sf-2)*5 + 4 * 4 ) + & !solid and fluid
           (nx-2)*(nr_w-2)*5 + &!wall inner
           2*(nx-2)*4 + &!buttom lines
           2*(nr_w-2)*4 + &!side lines
           4*3  + &!wall corners         
           2 * nx + &!HT betw. wall and fluid
           2 * nx !HT betw. wall and solid
!write(*,*) '# non-zero elements',nNonZero

!::Allocate sparse array and cols and rows
allocate( A(nNonZero) )
allocate( B(nx*(2*nr_sf+nr_w)), X(nx*(2*nr_sf+nr_w)) )
allocate( cols(nNonZero), rows(nx*(2*nr_sf+nr_w)+1) )


A(:) = 0
B(:) = 0
X(:) = 0 
cols(:) = 0
rows(:) = 0

rows(nx*(2*nr_sf+nr_w)+1) = nNonZero+1


!::The first row is the equation for the lower left fluid corner
rows(1) = 1
!central node
A( 1 ) =  por * rhof(1,1) * cf(1,1) / dt + &!capacity term
             2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kfRad(1,1) + dr(1+1)/kfRad(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/kfAx(1,1) + dx(1+1)/kfAx(1+1,1) ) ) + &!conduction to right node   
             h(1,1) * as !HT with the solid
cols( 1 ) = 1

!fluid part, right
A( 2 ) = axial * -2 / ( dx(1) * ( dx(1)/kfAx(1,1) + dx(1+1)/kfAx(1+1,1) ) ) !conduction to right node
cols( 2 ) = 2

!upper node
A( 3 ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kfRad(1,1) + dr(1+1)/kfRad(1,1+1) )) !conduction to node r+dr
cols(3) = 1 + nx


A( 4 ) = -h(1,1) * as !HT with the solid
cols(4) = 1 + nx * nr_sf

!convection (cols has already been set!!!)
if ( u(1,1) .gt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( 1 ) = A( 1 ) + rhof(1,1) * cf(1,1) * u(1,1) * por / (2 * dx(1)) 
    
    A( 2 ) = A( 2 ) + rhof(1,1) * cf(1,1) * u(1,1) * por / (2 * dx(1)) 
    
elseif ( u(1,1) .lt. 0 ) then
    !outlet
    !central node
    A( 1 ) = A( 1 ) - rhof(1,1) * cf(1,1) * u(1,1) * por / dx(1) 
    !right node
    A( 2 ) = A( 2 ) + rhof(1,1) * cf(1,1) * u(1,1) * por / dx(1) 
endif


!::The lower right corner
ind = 4 + (nx-2) * 5 + 1 !The lower left corner plus the line-nodes between the corners
rows(nx) = ind
!fluid part, central node
A( ind+1 ) = por * rhof( nx, 1 ) * cf( nx, 1 ) / dt + &!capacity term
              2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kfRad(nx,1) + dr(1+1)/kfRad(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/kfAx(nx,1) + dx(nx-1)/kfAx(nx-1,1) ) ) + &!conduction to left node
              h(nx, 1 ) * as 
cols(ind+1) = nx

!fluid part, left node
A( ind ) = axial * -2 / ( dx(nx) * ( dx(nx)/kfAx(nx,1) + dx(nx-1)/kfAx(nx-1,1) ) ) !conduction to the left node
cols(ind) = cols(ind+1)-1

!fluid part, upper node
A( ind+2 ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/kfRad(nx,1) + dr(1+1)/kfRad(nx,1+1) )) !conduction to node r+dr
cols(ind+2) = cols(ind+1)+nx

A( ind+3 ) = -h(nx,1) * as  !HT with the fluid
cols(ind+3) = cols(ind+1) + nx * nr_sf

if ( u(nx,1) .gt. 0 ) then
    !outlet
    !central node
    A( ind+1 ) = A( ind+1 ) + rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( dx(nx) ) 
    !left node
    A( ind ) = A( ind ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( dx(nx) ) 
elseif ( u(nx,1) .lt. 0 ) then
    !inlet, (then an explicit term is added to the B vector)
    !central node
    A( ind+1 ) = A( ind+1 ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( 2 * dx(nx) ) 
    !left node
    A( ind ) = A( ind ) - rhof(nx,1) * cf(nx,1) * u(nx,1) * por / ( 2 * dx(nx) ) 
endif

do i=2,nx-1
!::The lower line w/o corners, i.e. j=1, 2 <= i <= nx-1   
   !fluid part, central node   
   ind = 4 + 5 * (i-2) + 1 !lower left corner off set plus five for each of the elements in the line
   rows(i) = ind
   
   A( ind+1 ) = por * rhof(i,1) * cf(i,1) / dt + &!capacity term
                 2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/kfRad(i,1) + dr(1+1)/kfRad(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,1) + dx(i+1)/kfAx(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,1) + dx(i-1)/kfAx(i-1,1) ) ) + &!conduction to left node
                 h(i,1) * as !heat transfer with the solid
   cols( ind+1 ) = i

   !fluid part, left node
   A( ind ) = axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,1) + dx(i-1)/kfAx(i-1,1) ) ) !conduction to left node
   cols(ind) = cols(ind+1)-1

   !fluid part, right node
   A( ind+2 ) = axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,1) + dx(i+1)/kfAx(i+1,1) ) ) !conduction to right node
   cols(ind+2) = cols(ind+1)+1
   
   !fluid part, upper node
   A( ind + 3 ) = -2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/kfRad(i,1) + dr(1+1)/kfRad(i,1+1) ) ) !conduction to node r+dr
   cols(ind+3) = cols(ind+1) + nx

   !heat transfer between solid and fluid
   A( ind+4 ) = - h(i,1) * as 
   cols(ind+4) = cols(ind+1) + nx * nr_sf
   
   !convection   
   !left node
   A( ind ) = A( ind ) - rhof(i,1) * cf(i,1) * u(i,1) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
   !right node
   A( ind + 2 ) = A( ind + 2 ) + rhof(i,1) * cf(i,1) * u(i,1) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 


!:: Inner fluid nodes (still looping over i!)
    do j=2,nr_sf-1
        
        !::lower left and right corners offset (each has four) 
        !::plus lower line (each element contributes with five)
        !::plus the first elements of the left lines (five)
        !::plus each inner line (contributing with (nx-2)*6 each)
        !::Plus the trail of the current line (each with 6)
        !::plus the right line (contributing with five for each line)
        ind = 2*4 + 5 * (nx-2) + 5 * (j-1) + (nx-2)*6*(j-2) + (i-2)*6 + 5 * (j-2) + 1
        
        rows((j-1)*nx+i) = ind
        !fluid part, central nodes
        A( ind+2 ) =  por * rhof(i,j) * cf(i,j) / dt + &!capacity term
                 2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(i,j) + dr(j+1)/kfRad(i,j+1) )) + &!conduction to node r+dr
                 2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(i,j) + dr(j-1)/kfRad(i,j-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,j) + dx(i+1)/kfAx(i+1,j) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,j) + dx(i-1)/kfAx(i-1,j) ) ) + &!conduction to left node
                 h(i,j) * as !HT btw. solid and fluid
        cols(ind+2) = (j-1) * nx + i
        
        !fluid part, lower
        A( ind ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(i,j) + dr(j-1)/kfRad(i,j-1) ) ) !conduction to node r-dr
        cols(ind) = cols(ind+2) - nx
        
        !fluid part, left
        A( ind+1 ) = axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,j) + dx(i-1)/kfAx(i-1,j) ) ) !conduction to left node
        cols(ind+1) = cols(ind+2)-1
        
        !fluid part, right
        A( ind + 3 ) = axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,j) + dx(i+1)/kfAx(i+1,j) ) ) !conduction to right node
        cols(ind+3) = cols(ind+2)+1
        
        !fluid part, upper
        A( ind+4 ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(i,j) + dr(j+1)/kfRad(i,j+1) )) !conduction to node r+dr
        cols(ind+4) = cols(ind+2) + nx       
        

        A( ind+5 ) = -h(i,j) * as !HT btw. solid and fluid
        cols(ind+5) = cols(ind+2) + nx*nr_sf
        !convection
        !left node
        A( ind+1 ) = A( ind+1 ) - rhof(i,j) * cf(i,j) * u(i,j) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
        !right node
        A( ind+3 ) = A( ind+3 ) + rhof(i,j) * cf(i,j) * u(i,j) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
    enddo !end of inner fluid nodes
    
    !::The upper line, i.e. j=nr_sf,2<=i<=nx-1
    !::The two lower corners (each w. 4)
    !::Upper left corners (5)
    !::plus the lower line (each w. 5)
    !::plus the two side lines (each w. 5)
    !::plus inner nodes (each w. 6)
    !::The trailing part of the upper line (each w. 6)
    ind = 2 * 4 + 5 + 5 * (nx-2) + (nr_sf-2) * 2 * 5 + (nr_sf-2)*(nx-2)*6 + (i-2)*6 + 1
    rows(nx*(nr_sf-1)+i) = ind
   !fluid part, central node
    A( ind+2 ) = por * rhof(i,nr_sf) * cf(i,nr_sf) / dt + &!capacity term
                 2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(i,nr_sf) + dr(nr_sf-1)/kfRad(i,nr_sf-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,nr_sf) + dx(i+1)/kfAx(i+1,nr_sf) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kfAx(i,nr_sf) + dx(i-1)/kfAx(i-1,nr_sf) ) ) + &!conduction to left node
                 h(i,nr_sf) * as + &!HT btw. solid and fluid
                 walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/kfBdry(i) + 1/hw(i) ) ) !HT with the wall
                 
    cols(ind+2) = nx * (nr_sf-1) + i
    
    !fluid part, lower node
    A( ind ) = &
       -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(i,nr_sf) + dr(nr_sf-1)/kfRad(i,nr_sf-1) ) ) !conduction to node r-dr
    cols(ind) = cols(ind+2) - nx   
    
    !fluid part, left node
    A( ind + 1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,nr_sf) + dx(i-1)/kfAx(i-1,nr_sf) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2) - 1
   
   !fluid part, right node
    A( ind+3 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kfAx(i,nr_sf) + dx(i+1)/kfAx(i+1,nr_sf) ) ) !conduction to right node
    cols(ind+3) = cols(ind+2) + 1
   
     !HT btw. solid and fluid
    A( ind+4 ) = -h(i,nr_sf) * as
    cols(ind+4) = cols(ind+2) + nx * nr_sf
    
    !fluid part, HT with the wall
    A( ind+5 ) =&
        -2 * walls / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/kfBdry(i) + 1/hw(i) ) ) !HT with the wall
    cols(ind+5) = cols(ind+2) + nx*nr_sf + nx
   !convection
   !Left node
   A( ind+1 ) = A( ind + 1 ) - rhof(i,nr_sf) * cf(i,nr_sf) * u(i,nr_sf) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 
   !right node
   A( ind+3 ) = A( ind + 3 ) + rhof(i,nr_sf) * cf(i,nr_sf) * u(i,nr_sf) * por / ( dx(i) + 0.5*(dx(i-1)+dx(i+1)) ) 

    
enddo !end of i-loop for the fluid nodes

!!Upper left corner
!::The two lower corners (each w. 4)
!::plus the lower line (each w. 5)
!::plus the two side lines (each w. 5)
!::plus inner nodes (each w. 6)
ind = 2 * 4 + 5 * (nx-2) + (nr_sf-2) * 2 * 5 + (nr_sf-2)*(nx-2)*6 + 1
rows( nx * (nr_sf-1) + 1 ) = ind
!fluid part, central node
A( ind+1 ) = por * rhof( 1, nr_sf ) * cf( 1, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(1,nr_sf) + dr(nr_sf-1)/kfRad(1,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/kfAx(1,nr_sf) + dx(1+1)/kfAx(1+1,nr_sf) ) ) + &!conduction to right node
                h(1, nr_sf) * as + &!HT with the solid                
                walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/kfBdry(1) + 1/hw(1) ) ) !HT with the wall
cols(ind+1) = nx * (nr_sf-1) + 1
!fluid part, lower node
A( ind ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(1,nr_sf) + dr(nr_sf-1)/kfRad(1,nr_sf-1) ) ) !conduction to node r-dr
cols(ind) = cols(ind+1) - nx
!fluid part, right node
A( ind+2 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/kfAx(1,nr_sf) + dx(1+1)/kfAx(1+1,nr_sf) ) ) !conduction to right node
cols(ind+2) = cols(ind+1) + 1

!HT with the solid
A( ind+3 ) = -h(1, nr_sf ) * as 
cols(ind+3) = cols(ind+1) + nx*nr_sf

!fluid part, HT with the wall
A( ind+4 ) = -2 * walls / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/kfBdry(1) + 1/hw(1) ) ) 
cols(ind+4) = cols(ind+1) + nx*nr_sf+nx


!convection 
 
if ( u(1,nr_sf) .gt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( ind+1 ) = A( ind+1 ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / ( 2*dx(1) ) 
    !right node
    A( ind+2 ) = A( ind+2 ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / ( 2*dx(1) )   

elseif ( u(1,nr_sf) .lt. 0 ) then
    !outlet
    !central node
    A( ind+1 ) = A( ind+1 ) - rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / dx(1) 
    !right node
    A( ind+2 ) = A( ind+2 ) + rhof(1,nr_sf) * cf(1,nr_sf) * u(1,nr_sf) * por / dx(1) 
endif


!!Upper right corner
!::The two lower corners (each w. 4)
!::The upper left corner (5)
!::plus the lower line (each w. 5)
!::plus the two side lines (each w. 5)
!::plus inner nodes (each w. 6)
!::plus the upper line (each w. 6)
ind = 2 * 4 + 5 + 5 * (nx-2) + (nr_sf-2) * 2 * 5 + (nr_sf-2)*(nx-2)*6 + (nx-2)*6 +1
rows( nx * nr_sf ) = ind
!fluid part, central node
A( ind+2 ) = por * rhof( nx, nr_sf ) * cf( nx, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(nx,nr_sf) + dr(nr_sf-1)/kfRad(nx,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/kfAx(nx,nr_sf) + dx(nx-1)/kfAx(nx-1,nr_sf) ) ) + &!conduction to left node
                h( nx, nr_sf ) * as + &!HT with the solid
                walls * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/kfBdry(nx) + 1/hw(nx) ) ) !HT with the wall
cols(ind+2) = nx * nr_sf
!fluid part, lower node
A( ind ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/kfRad(nx,nr_sf) + dr(nr_sf-1)/kfRad(nx,nr_sf-1) ) ) !conduction to node r-dr
cols(ind) = cols(ind+2) - nx
!fluid part, left node
A( ind+1 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/kfAx(nx,nr_sf) + dx(nx-1)/kfAx(nx-1,nr_sf) ) ) !conduction to left node
cols(ind+1) = cols(ind+2) - 1
!HT with the solid
A( ind+3 ) = -h( nx, nr_sf ) * as 
cols(ind+3) = cols(ind+2) + nx*nr_sf
!HT with the wall
A( ind+4 ) = -2 * walls / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/kfBdry(nx) + 1/hw(nx) ) ) 
cols(ind+4) = cols(ind+2) + nx*nr_sf+nx
!convection 
if ( u(nx,nr_sf) .gt. 0 ) then
    !outlet
    !central node
    A( ind+2 ) = A( ind+2 ) + rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / dx(nx) 
    !left node
    A( ind+1 ) = A( ind+1 ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / dx(nx) 
elseif ( u(nx,nr_sf) .lt. 0 ) then
    !inlet (then an explicit term is added to the B vector)
    !central node
    A( ind+2 ) = A( ind+2 ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx,nr_sf) * por / ( 2 * dx(nx) ) 
    !left node
    A( ind+1 ) = A( ind+1 ) - rhof(nx,nr_sf) * cf(nx,nr_sf) * u(nx-1,nr_sf) * por / ( 2 * dx(nx)) 
endif

!Fluid left and right lines (i=1 or nx and 2<=j<=nr_sf-1
do j=2,nr_sf-1
    !side line one (i = 1)
    !::The two lower corners (each w. 4)
    !::The lower line (each ele. w. 5)
    !::The left line (each 5)
    !::The inner lines (each 6)
    !::The right line (each 5)
    ind = 2 * 4 + (nx-2) * 5 + (j-2)*5 + (j-2)*(nx-2)*6 + (j-2)*5 + 1
    rows( nx*(j-1)+1) = ind
    !fluid part, central node
    A( ind+1 ) = por * rhof(1,j) * cf(1,j) / dt + &!capacity term
        2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(i,j) + dr(j+1)/kfRad(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(i,j) + dr(j-1)/kfRad(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/kfAx(1,j) + dx(1+1)/kfAx(1+1,j) ) ) + &!conduction to right node
        h(1,j) * as !HT btw. solid and fluid
    cols(ind+1) = nx*(j-1)+1
    
    !fluid part, lower
    A( ind ) = &
        -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(i,j) + dr(j-1)/kfRad(i,j-1) ) ) !conduction to node r-dr
    cols(ind) = cols(ind+1) - nx
    
    !fluid part, right node
    A( ind+2 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/kfAx(1,j) + dx(1+1)/kfAx(1+1,j) ) ) !conduction to right node
    cols(ind+2) = cols(ind+1) + 1
    
    !fluid part, upper
    A( ind+3 ) = &
        -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(i,j) + dr(j+1)/kfRad(i,j+1) ) ) !conduction to node r+dr
    cols(ind+3) = cols(ind+1) + nx
    
    A( ind+4 ) = -h(1,j) * as !HT btw. solid and fluid
    cols(ind+4) = cols(ind+1) + nx*nr_sf
    
    !convection
    if ( u(1,j) .gt. 0 ) then
       !there will then be an explicit term in the B vector
       !central node   
       A( ind+1 ) = A( ind+1 ) + rhof(1,j) * cf(1,j) * u(1,j) * por / ( 2 * dx(1) ) 
       !right node
       A( ind+2 ) = A( ind+2 ) + rhof(1,j) * cf(1,j) * u(1,j) * por / ( 2 * dx(1) ) 
       
    else if ( u(1,j) .lt. 0 ) then
       !outlet
       !central node
       A( ind+1 ) = A( ind+1 ) - rhof(1,j) * cf(1,j) * u(1,j) * por / dx(1) 
       !right node
       A( ind+2 ) = A( ind+2 ) + rhof(1,j) * cf(1,j) * u(1,j) * por / dx(1) 
    endif
    
    !Right line, i=nx
    !::The two lower corners (each w. 4)
    !::The lower line (each ele. w. 5)
    !::The left line (each 5)
    !::The inner lines (each 6)
    !::The right line (each 5)
    ind = 2 * 4 + (nx-2) * 5 + (j-1)*5 + (j-1)*(nx-2)*6 + (j-2)*5 + 1
    rows( nx*j ) = ind  
   !fluid part, central node
    A( ind+2 ) = por * rhof( nx, j ) * cf( nx, j ) / dt + &!capacity term
                  2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(nx,j) + dr(j+1)/kfRad(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(nx,j) + dr(j-1)/kfRad(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/kfAx(nx,j) + dx(nx-1)/kfAx(nx-1,j) ) ) + &!conduction to left node
                  h(nx,j) * as 
    cols(ind+2) = nx*j

    !fluid part, lower
    A( ind ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/kfRad(nx,j) + dr(j-1)/kfRad(nx,j-1) )) !conduction to node r-dr
    cols(ind) = cols(ind+2)-nx

    !fluid part, left
    A( ind+1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/kfAx(nx,j) + dx(nx-1)/kfAx(nx-1,j) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2)-1

    !fluid part, upper
    A( ind+3 ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/kfRad(nx,j) + dr(j+1)/kfRad(nx,j+1) )) !conduction to node r+dr
    cols(ind+3) = cols(ind+2)+nx
       
    A( ind+4 ) = -h(nx,j) * as !HT btw. solid and fluid
    cols(ind+4) = cols(ind+2) + nx*nr_sf
    
    !convection    
    if ( u(nx,j) .gt. 0 ) then
        !outlet
        !central node
        A( ind+2 ) = A( ind+2 ) + rhof(nx,j) * cf(nx,j) * u(nx,j) * por / dx(nx) 
        !left node
        A( ind+1 ) = A( ind+1 ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( dx(nx) ) 
    elseif ( u(nx,j) .lt. 0 ) then
        !inlet (then an explicit term is added to the B vector
        !central node
        A( ind+2 ) = A( ind+2 ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( 2 * dx(nx)  ) 
        
        A( ind+1 ) = A( ind+1 ) - rhof(nx,j) * cf(nx,j) * u(nx,j) * por / ( 2 * dx(nx)  ) 
        
    endif
enddo

!::This is to keep things a bit simpler (the solid indices are offset with this value)
!::The two lower corners (each w. 4)
!::The two upper corners (5)
!::plus the lower line (each w. 5)
!::plus the two side lines (each w. 5)
!::plus inner nodes (each w. 6)
!::plus the upper line (each w. 6)
indFluidOffSet = 2 * 4 + 2 * 5 + 5 * (nx-2) + (nr_sf-2) * 2 * 5 + (nr_sf-2)*(nx-2)*6 + (nx-2)*6
!::Fluid part ends here
!::--------------------------------::!
!::Solid part begins here

!::Lower left corner
ind = indFluidOffSet + 1
rows( nx*nr_sf+1 ) = ind

!solid part, central node
A( ind+1 ) =  ( 1 - por ) * rhos(1,1) * cs(1,1) / dt + &!capacity term
             2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(1,1) + dr(1+1)/ks(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/ks(1,1) + dx(1+1)/ks(1+1,1) ) ) + &!conduction to right node                
             h(1,1) * as 
cols(ind+1) = nx * nr_sf + 1

!HT with the fluid
A( ind ) = -h(1,1) * as 
cols(ind) = cols(ind+1) - nx*nr_sf

!solid part, right
A( ind+2 ) = axial * -2 / ( dx(1) * ( dx(1)/ks(1,1) + dx(1+1)/ks(1+1,1) ) ) !conduction to right node
cols(ind+2) = cols(ind+1) + 1

!solid part, upper
A( ind+3 ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(1,1) + dr(1+1)/ks(1,1+1) )) !conduction to node r+dr
cols(ind+3) = cols(ind+1) + nx

!::Lower right corner
!::The offset
!::The lower left corner (4)
!::The lower line (each having 5)
ind = indFluidOffset + 4 + (nx-2)*5 + 1
rows(nx*nr_sf+nx) = ind

!solid part, central node
A( ind+2 ) = ( 1 - por ) * rhos( nx, 1 ) * cs( nx, 1 ) / dt + &!capacity term
              2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(nx,1) + dr(1+1)/ks(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,1) + dx(nx-1)/ks(nx-1,1) ) ) + & !conduction to left node   
              h(nx,1) * as 
cols(ind+2) = nx*nr_sf + nx

!HT with fluid
A( ind ) = -h(nx,1) * as 
cols(ind) = cols(ind+2) - nx*nr_sf

!solid part, left node
A( ind+1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,1) + dx(nx-1)/ks(nx-1,1) ) ) !conduction to the left node
cols(ind+1) = cols(ind+2) - 1

!solid part, upper node
A( ind+3 ) = -2 / ( r(1)*dr(1) / rup(1) * ( dr(1)/ks(nx,1) + dr(1+1)/ks(nx,1+1) )) !conduction to node r+dr
cols(ind+3) = cols(ind+2) + nx


do i=2,nx-1

    !::Lower line,i.e. j=1
    !::The offset
    !::The lower left corner (4)
    !::The trail of the lower line
    ind = indFluidOffset + 4 + (i-2)*5 + 1
    rows(nx*nr_sf + i) = ind
    !solid part, central node
    A( ind+2 ) = (1 - por) * rhos(i,1) * cs(i,1) / dt + &!capacity term
                 2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/ks(i,1) + dr(1+1)/ks(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i+1)/ks(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i-1)/ks(i-1,1) ) ) + &!conduction to left node
                 h(i,1) * as !heat transfer between fluid and solid
    cols(ind+2) = nx*nr_sf+i

    !HT with fluid
    A( ind ) = -h(i,1) * as !heat transfer between fluid and solid
    cols(ind) = cols(ind+2) - nx*nr_sf

    !solid part, left node
    A( ind+1 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i-1)/ks(i-1,1) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2) - 1

    !solid part, right node
    A( ind+3 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,1) + dx(i+1)/ks(i+1,1) ) ) !conduction to right node
    cols(ind+3) = cols(ind+2) + 1

    !solid part, upper node
    A( ind+4 ) = -2 / ( r(1) * dr(1) / rup(1) * ( dr(1)/ks(i,1) + dr(1+1)/ks(i,1+1) ) ) !conduction to node r+dr
    cols(ind+4) = cols(ind+2) + nx

    !inner nodes
    do j=2,nr_sf-1
        !::The offset
        !::Lower corners (each 4)     
        !::Lower line (each 5)
        !::left line
        !::Right line
        !::Inner nodes
        !::Current inner line
        ind = indFluidOffset + 4*2 + (nx-2)*5 + (j-1)*5 + (j-2) * 5 + (j-2)*(nx-2)*6 + (i-2)*6 + 1        
        rows(nx*nr_sf + (j-1)*nx + i ) = ind
        ! solid part, inner nodes
        A( ind+3 ) =  (1 - por) * rhos(i,j) * cs(i,j) / dt + &!capacity term
                 2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) + &!conduction to node r+dr
                 2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i+1)/ks(i+1,j) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i-1)/ks(i-1,j) ) ) + &!conduction to left node
                 h(i,j) * as !HT btw. solid and fluid
        cols(ind+3) = nx*nr_sf + (j-1)*nx + i

        A( ind ) = -h(i,j) * as !HT btw. solid and fluid
        cols(ind) = cols(ind+3) - nx*nr_sf

        !solid part, lower
        A( ind+1 ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) !conduction to node r-dr
        cols(ind+1) = cols(ind+3) - nx

        !solid part, left
        A( ind+2 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i-1)/ks(i-1,j) ) ) !conduction to left node
        cols(ind+2) = cols(ind+3)-1

        !solid part, right
        A( ind+4 ) = axial * -2 / ( dx(i) * ( dx(i)/ks(i,j) + dx(i+1)/ks(i+1,j) ) ) !conduction to right node
        cols(ind+4) = cols(ind+3) + 1

        !solid part, upper
        A( ind+5 ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) )) !conduction to node r+dr
        cols(ind+5) = cols(ind+3) + nx
     
    enddo
   
    !::The upper line, i.e. j=nr_sf
    !::The offset    
    !::Lower corners (each 4)
    !::Upper left corner (5)
    !::Lower line (each 5)
    !::Side lines (two lines each w. 5)
    !::Central nodes (each 6)
    !::Trail of upper line (each 6)
    ind = indFluidOffset + 2*4 + 5 + (nx-2)*5 + (nr_sf-2) * 2 * 5 + (nr_sf-2)*(nx-2)*6 + (i-2)*6 + 1
    rows(nx*nr_sf + nx*(nr_sf-1)+i) = ind
    
    !solid part, central node
    A( ind+3 ) = (1 - por) * rhos(i,nr_sf) * cs(i,nr_sf) / dt + &!capacity term
                 2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(i,nr_sf) + dr(nr_sf-1)/ks(i,nr_sf-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i+1)/ks(i+1,nr_sf) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i-1)/ks(i-1,nr_sf) ) ) + &!conduction to left node
                 h(i,nr_sf) * as + &!HT btw. solid and fluid
                 wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/ks(i,nr_sf) + 1/hw(i) ) ) !HT with the wall
    cols(ind+3) = nx*nr_sf + nx*(nr_sf-1)+i
    
    !HT btw. solid and fluid
    A( ind ) = -h(i,nr_sf) * as 
    cols(ind) = cols(ind+3) - nx*nr_sf
   
    !solid part, lower node
    A( ind+1 ) = &
       -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(i,nr_sf) + dr(nr_sf-1)/ks(i,nr_sf-1) ) ) !conduction to node r-dr
    cols(ind+1) = cols(ind+3) - nx

    !solid part, left node
    A( ind+2 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i-1)/ks(i-1,nr_sf) ) ) !conduction to left node
    cols(ind+2) = cols(ind+3) - 1
    
    !solid part, right node
    A( ind+4 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/ks(i,nr_sf) + dx(i+1)/ks(i+1,nr_sf) ) ) !conduction to right node
    cols(ind+4) = cols(ind+3) + 1 
   
    !HT with the wall
    A( ind+5 ) = &
    -wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(i,1) + dr(nr_sf)/ks(i,nr_sf) + 1/hw(i) ) ) !HT with the wall  
    cols(ind+5) = cols(ind+3) + nx
enddo

!::The two sidelines
do j=2,nr_sf-1
    !::Left line, i.e. i=1   
    !::The offset    
    !::Lower corners (each 4)        
    !::Lower line (each 5)
    !::Side lines (two lines each w. 5)
    !::Central nodes (each 6)
    ind = indFluidOffset + 2*4 + (nx-2)*5 + (j-2)*5 + (j-2)*(nx-2)*6 + (j-2)*5 + 1
    rows( nx*nr_sf + nx * (j-1) + 1 ) = ind    
    !solid part, central node
    A( ind+2 ) = (1 - por) * rhos(1,j) * cs(1,j) / dt + &!capacity term
        2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/ks(1,j) + dx(1+1)/ks(1+1,j) ) ) + &!conduction to right node
        h(1,j) * as 
    cols(ind+2) = nx*nr_sf + nx * (j-1) + 1
    
    !HT btw. solid and fluid
    A( ind ) = -h(i,1) * as
    cols(ind) = cols(ind+2) - nx*nr_sf
    
    !solid part, lower
    A( ind+1 ) = &
        -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(i,j) + dr(j-1)/ks(i,j-1) ) ) !conduction to node r-dr
    cols(ind+1) = cols(ind+2) - nx
    
    !solid part, right node
    A( ind+3 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/ks(1,j) + dx(1+1)/ks(1+1,j) ) ) !conduction to right node
    cols(ind+3) = cols(ind+2) + 1
    
    !solid part, upper
    A( ind+4 ) = &
        -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(i,j) + dr(j+1)/ks(i,j+1) ) ) !conduction to node r+dr
    cols(ind+4) = cols(ind+2) + nx    
    
    
    !::Right line, i.e. i=nx    
    !::The offset    
    !::Lower corners (each 4)        
    !::Lower line (each 5)
    !::Left line (each w. 5)
    !::Central nodes (each 6)
    !::Right line (each w. 5)
    ind = indFluidOffset + 2*4 + (nx-2)*5 + (j-1)*5 + (j-1)*(nx-2)*6 + (j-2)*5 + 1
    rows( nx*nr_sf + nx * j ) = ind       
    !solid part, central node
    A( ind+3 ) = ( 1 - por ) * rhos( nx, j ) * cs( nx, j ) / dt + &!capacity term
                  2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(nx,j) + dr(j+1)/ks(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(nx,j) + dr(j-1)/ks(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,j) + dx(nx-1)/ks(nx-1,j) ) ) + &!conduction to left node
                  h(nx, j) * as 
    cols(ind+3) = nx * nr_sf + nx * j
    
    !HT btw. solid and fluid
    A( ind ) = -h(nx,j) * as 
    cols(ind) = cols(ind+3) - nx*nr_sf
    
    !solid part, lower
    A( ind+1 ) = -2 / ( r(j)*dr(j) / rdn(j) * ( dr(j)/ks(nx,j) + dr(j-1)/ks(nx,j-1) )) !conduction to node r-dr
    cols(ind+1) = cols(ind+3) - nx
    
    !solid part, left
    A( ind+2 ) = axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,j) + dx(nx-1)/ks(nx-1,j) ) ) !conduction to left node
    cols(ind+2) = cols(ind+3) - 1
    
    !solid part, upper
    A( ind+4 ) = -2 / ( r(j)*dr(j) / rup(j) * ( dr(j)/ks(nx,j) + dr(j+1)/ks(nx,j+1) )) !conduction to node r+dr
    cols(ind+4) = cols(ind+3) + nx
enddo


!::Upper left corner
!::The offset    
!::Lower corners (each 4)        
!::Lower line (each 5)
!::Left and right lines (two each w. 5)
!::Central nodes (each 6)
ind = indFluidOffset + 2*4 + (nx-2)*5 + (nr_sf-2)*5*2 + (nr_sf-2)*(nx-2)*6 + 1
rows( nx*nr_sf + nx * (nr_sf-1) + 1 ) = ind
!solid part, central node
A( ind+2 ) = (1 - por) * rhos( 1, nr_sf ) * cs( 1, nr_sf ) / dt + &!capacity term                    
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(1,nr_sf) + dr(nr_sf-1)/ks(1,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/ks(1,nr_sf) + dx(1+1)/ks(1+1,nr_sf) ) ) + &!conduction to right node
                h(1, nr_sf ) * as + & !HT with the fluid
                wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/ks(1,nr_sf) + 1/hw(1) ) ) !HT with the wall
cols(ind+2) = nx*nr_sf + nx*(nr_sf-1) + 1

!HT with the fluid
A( ind ) = -h(1, nr_sf ) * as 
cols(ind) = cols(ind+2) - nx*nr_sf

!solid part, lower node
A( ind+1 ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(1,nr_sf) + dr(nr_sf-1)/ks(1,nr_sf-1) ) ) !conduction to node r-dr
cols(ind+1) = cols(ind+2) - nx

!solid part, right node
A( ind+3 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/ks(1,nr_sf) + dx(1+1)/ks(1+1,nr_sf) ) ) !conduction to right node
cols(ind+3) = cols(ind+2) + 1

!HT with the wall
A( ind+4 ) = -wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(1,1) + dr(nr_sf)/ks(1,nr_sf) + 1/hw(1) ) ) !HT with the wall
cols(ind+4) = cols(ind+2) + nx
!::Upper right corner
!::The offset    
!::Lower corners (each 4)
!::Upper left corner (5)        
!::Lower line (each 5)
!::Upper line (each 6)
!::Left and right lines (two each w. 5)
!::Central nodes (each 6)
ind = indFluidOffset + 2*4 + 5 + (nx-2)*5 + (nx-2)*6 + (nr_sf-2)*5*2 + (nr_sf-2)*(nx-2)*6 + 1
rows( nx*nr_sf + nx * nr_sf ) = ind

!solid part, central node
A( ind+3 ) = (1 - por) * rhos( nx, nr_sf ) * cs( nx, nr_sf ) / dt + &!capacity term
                2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(nx,nr_sf) + dr(nr_sf-1)/ks(nx,nr_sf-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/ks(nx,nr_sf) + dx(nx-1)/ks(nx-1,nr_sf) ) ) + &!conduction to left node
                h( nx, nr_sf ) * as + & !HT with the fluid
                wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/ks(nx,nr_sf) + 1/hw(nx) ) ) !HT with the wall
cols(ind+3) = nx * nr_sf + nx * nr_sf

!HT with the fluid
A( ind ) = -h( nx, nr_sf ) * as 
cols(ind) = cols(ind+3) - nx * nr_sf

!solid part, lower node
A( ind+1 ) = &
    -2 / ( r(nr_sf) * dr(nr_sf) / rdn(nr_sf) * ( dr(nr_sf)/ks(nx,nr_sf) + dr(nr_sf-1)/ks(nx,nr_sf-1) ) ) !conduction to node r-dr
cols(ind+1) = cols(ind+3) - nx

!solid part, left node
A( ind+2 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/ks(nx,nr_sf) + dx(nx-1)/ks(nx-1,nr_sf) ) ) !conduction to left node
cols(ind+2) = cols(ind+3) - 1

!HT with the wall
A( ind+4 ) = &
-wallsSol * 2 / ( r(nr_sf) * dr(nr_sf) / rup(nr_sf) * ( dr(nr_sf+1)/kw(nx,1) + dr(nr_sf)/ks(nx,nr_sf) + 1/hw(nx) ) ) !HT with the wall
cols(ind+4) = cols(ind+3) + nx
!::End of solid part::!
!::-------------------------::!
!::Begin of wall part::!
!::Index offset for the wall
!::The fluid off set
!::The four corners in the solid
!::The lower and upper lines in the solid
!::The left and right lines in the solid
!::The inner nodes in the solid
indFluidSolidOffSet = indFluidOffSet + 2*4 + 2*5 + (nx-2)*5 + (nx-2)*6 + (nr_sf-2)*5*2 + (nr_sf-2)*(nx-2)*6
!::Lower left corner of the wall

ind = indFluidSolidOffSet + 1
rows( 2 * nx*nr_sf + 1 ) = ind
jj = nr_sf 
!wall part, central node
A( ind+2 ) = rhow(1,1) * cw(1,1) / dt + &!capacity term
             2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(1,1) + dr(jj+1+1)/kw(1,1+1) )) + &!conduction to node r+dr           
             axial * 2 / ( dx(1) * ( dx(1)/kw(1,1) + dx(1+1)/kw(1+1,1) ) ) + &!conduction to right node                              
             walls * 2 / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/kfBdry(1) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) + &!HT with the fluid
             wallsSol * 2 / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/ks(1,nr_sf) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) !HT with the solid
cols(ind+2) = 2 * nx*nr_sf + 1

!HT with the fluid
A( ind ) = -2 * walls / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/kfBdry(1) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) 
cols(ind) = cols(ind+2) - nx * nr_sf - nx

!HT with the solid
A( ind+1 ) = -2 * wallsSol / ( r(jj+1)*dr(jj+1) / rdn(jj+1) * ( dr(jj+1-1)/ks(1,nr_sf) + dr(jj+1)/kw(1,1) + 1/hw(1) ) ) 
cols(ind+1) = cols(ind+2) - nx

!wall part, right
A( ind+3 ) = axial * -2 / ( dx(1) * ( dx(1)/kw(1,1) + dx(1+1)/kw(1+1,1) ) ) !conduction to right node
cols(ind+3) = cols(ind+2) + 1

!wall part, upper
A( ind+4 ) = -2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(1,1) + dr(jj+1+1)/kw(1,1+1) )) !conduction to node r+dr
cols(ind+4) = cols(ind+2) + nx

!::Lower right corner, wall
!::The offset
!::The lower left corner
!::The lower line (each w. 6)
ind = indFluidSolidOffSet + 5 + 6 * (nx-2) + 1
rows( 2*nx*nr_sf+nx ) = ind
jj = nr_sf 
!wall part, central node
A( ind+3 ) = rhow( nx, 1 ) * cw( nx, 1 ) / dt + &!capacity term
              2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(nx,1) + dr(jj+1+1)/kw(nx,1+1) )) + &!conduction to node r+dr           
              axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,1) + dx(nx-1)/kw(nx-1,1) ) ) + &!conduction to left node
              walls * 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/kfBdry(nx) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) + & !HT with the fluid
              wallsSol * 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/ks(nx,nr_sf) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) !HT with the fluid
cols(ind+3) = 2*nx*nr_sf + nx

!HT with the fluid
A( ind ) = -2 * walls / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/kfBdry(nx) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) 
cols(ind) = cols(ind+3) - nx*nr_sf - nx

!HT with the solid
A( ind+1 ) = -2 * wallsSol / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj)/ks(nx,nr_sf) + dr(jj+1)/kw(nx,1) + 1/hw(nx)) ) 
cols(ind+1) = cols(ind+3) - nx

!wall part, left node
A( ind+2 ) = axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,1) + dx(nx-1)/kw(nx-1,1) ) ) !conduction to the left node
cols(ind+2) = cols(ind+3) - 1

!wall part, upper node
A( ind+4 ) = -2 / ( r(jj+1)*dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(nx,1) + dr(jj+1+1)/kw(nx,1+1) )) !conduction to node r+dr
cols(ind+4) = cols(ind+3) + nx


!::Lower and upper lines and inner nodes
do i=2,nx-1

   !::Lower line, i.e. j=1
   !::The offset
   !::The lower left corner
   !::The lower line (each w. 6)
    ind = indfluidSolidOffSet + 5 + 6 * (i-2) + 1
    rows( 2*nx*nr_sf + i ) = ind
   
    jj = nr_sf 
   !wall part, central node
    A( ind+3 ) = rhow(i,1) * cw(i,1) / dt + &!capacity term
                 2 / ( r(jj+1) * dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1+1)/kw(i,1+1) ) ) + &!conduction to node r+dr                 
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i+1)/kw(i+1,1) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i-1)/kw(i-1,1) ) ) + &!conduction to left node
                 walls * 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/kfBdry(i) + 1/hw(i) ) ) + &!HT with the fluid
                 wallsSol * 2 / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/ks(i,nr_sf) + 1/hw(i) ) ) !HT with the fluid
    cols(ind+3) = 2*nx*nr_sf + i
   
   !HT with the fluid
    A( ind ) = -2 * walls / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/kfBdry(i) + 1/hw(i) ) ) 
    cols(ind) = cols(ind+3) - nx * nr_sf - nx
    
    !HT with the solid
    A( ind+1 ) = -2 * wallsSol / ( r(jj+1) * dr(jj+1) / rdn(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1-1)/ks(i,nr_sf) + 1/hw(i) ) ) 
    cols(ind+1) = cols(ind+3) - nx
   
   !wall part, left node
    A( ind+2 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i-1)/kw(i-1,1) ) ) !conduction to left node
    cols(ind+2) = cols(ind+3) - 1
   
   !wall part, right node
    A( ind+4 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,1) + dx(i+1)/kw(i+1,1) ) ) !conduction to right node
    cols(ind+4) = cols(ind+3) + 1
    
    !wall part, upper node
    A( ind+5 ) = -2 / ( r(jj+1) * dr(jj+1) / rup(jj+1) * ( dr(jj+1)/kw(i,1) + dr(jj+1+1)/kw(i,1+1) ) ) !conduction to node r+dr
    cols(ind+5) = cols(ind+3) + nx
   
   
    !::Upper line, i.e. j=nr_w
    !::The offset
    !::The lower corners (each w. 5)
    !::The lower line (each w. 6)
    !::The two sidelines (i=1,nx respectively, each w. 4)
    !::The inner nodes, each w. 5)
    !::Upper left corner (3)
    !::Trail of upper line
    ind = indFluidSolidOffSet + 5*2 + 6 * (nx-2) + (nr_w-2)*4*2 + (nx-2)*(nr_w-2)*5 + 3 + (i-2)*4 + 1
    rows( 2*nx*nr_sf + nx*(nr_w-1)+i ) = ind
   
    jj = nr_sf 
    !wall part, central node
    A( ind+2 ) = rhow(i,nr_w) * cw(i,nr_w) / dt + &!capacity term
                 2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(i,nr_w) + dr(jj+nr_w-1)/kw(i,nr_w-1) ) ) + &!conduction to node r-dr
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i+1)/kw(i+1,nr_w) ) ) + &!conduction to right node
                 axial * 2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i-1)/kw(i-1,nr_w) ) ) !conduction to left node
    cols(ind+2) = 2*nx*nr_sf + nx*(nr_w-1) + i
    
    !wall part, lower node
    A( ind ) = &
        -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(i,nr_w) + dr(jj+nr_w-1)/kw(i,nr_w-1) ) ) !conduction to node r-dr
    cols(ind) = cols(ind+2) - nx
    
    !wall part, left node
    A( ind+1 ) = &
       axial * -2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i-1)/kw(i-1,nr_w) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2) - 1
    
    !wall part, right node
    A( ind+3 ) = &
        axial * -2 / ( dx(i) * ( dx(i)/kw(i,nr_w) + dx(i+1)/kw(i+1,nr_w) ) ) !conduction to right node
    cols(ind+3) = cols(ind+2) + 1
    
   !::wall, inner nodes
   do j=2,nr_w-1
    !::Lower line, i.e. j=1
    !::The offset
    !::The lower corners (each w. 5)
    !::The lower line (each w. 6)
    !::The left line (each w. 4)
    !::Pass inner lines (each w. 5)
    !::Trailing, current, inner line (each w. 5)
    !::The right line (each w. 4)
    ind = indFluidSolidOffset + 2*5 + (nx-2) * 6 + (j-1)*4 + 5 * (j-2)*(nx-2) + (i-2)*5 + (j-2)*4 + 1
    rows( 2*nx*nr_sf + (j-1)*nx + i ) = ind
    jj = nr_sf 
    ! wall part, inner nodes
    A( ind+2 ) =   rhow(i,j) * cw(i,j) / dt + &!capacity term
                     2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) + &!conduction to node r+dr
                     2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) + &!conduction to node r-dr
                     axial * 2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i+1)/kw(i+1,j) ) ) + &!conduction to right node
                     axial * 2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i-1)/kw(i-1,j) ) ) !conduction to left node                         
    cols(ind+2) = 2*nx*nr_sf + (j-1)*nx + i

    !wall part, lower
    A( ind ) = -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) !conduction to node r-dr
    cols(ind) = cols(ind+2) - nx
    
    !wall part, left
    A( ind+1 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i-1)/kw(i-1,j) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2) - 1
    
    !wall part, right
    A( ind+3 ) = axial * -2 / ( dx(i) * ( dx(i)/kw(i,j) + dx(i+1)/kw(i+1,j) ) ) !conduction to right node
    cols(ind+3) = cols(ind+2) + 1
    
    !wall part, upper
    A( ind+4 ) = -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) )) !conduction to node r+dr
    cols(ind+4) = cols(ind+2) + nx
   enddo
enddo


!the two side lines for the wall
do j=2,nr_w-1
    !::Left most line, i=1
    !::The offset
    !::The two lower corners (each w. 5)
    !::The lower line (each w. 6)
    !::The two side lines (each w. 4)
    !::The inner nodes (each w. 5)
    ind = indFluidSolidOffSet + 2*5 + (nx-2) * 6 + (j-2)*4*2 + (nx-2)*(j-2)*5 + 1
    rows( 2*nx*nr_sf + (j-1)*nx + 1) = ind
    jj = nr_sf 
    !wall part, central node
    A( ind+1 ) = rhow(1,j) * cw(1,j) / dt + &!capacity term
        2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) + &!conduction to node r+dr
        2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) + &!conduction to node r-dr
        axial * 2 / ( dx(1) * ( dx(1)/kw(1,j) + dx(1+1)/kw(1+1,j) ) ) !conduction to right node            
    cols(ind+1) = 2*nx*nr_sf + (j-1)*nx + 1
    
    !wall part, lower
    A( ind ) = &
        -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j-1)/kw(i,j-1) ) ) !conduction to node r-dr
    cols(ind) = cols(ind+1) - nx
    
    !wall part, right node
    A( ind+2 ) = &
        axial * -2 / ( dx(1) * ( dx(1)/kw(1,j) + dx(1+1)/kw(1+1,j) ) ) !conduction to right node
    cols(ind+2) = cols(ind+1) + 1
    
    !wall part, upper
    A( ind+3 ) = &
        -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(i,j) + dr(jj+j+1)/kw(i,j+1) ) ) !conduction to node r+dr
    cols(ind+3) = cols(ind+1) + nx
    
    
    !::Right side line (i=nx)
    !::The offset
    !::The two lower corners (each w. 5)
    !::The lower line (each w. 6)
    !::The left side line (each w. 4)
    !::The right side line (each w. 4)
    !::The inner nodes (each w. 5)
    ind = indFluidSolidOffSet + 2*5 + (nx-2) * 6 + (j-1)*4 + (j-2)*4 + (nx-2)*(j-1)*5 + 1
    rows( 2*nx*nr_sf + j*nx ) = ind
    
    !wall part, central node
    A( ind+2 ) = rhow( nx, j ) * cw( nx, j ) / dt + &!capacity term
                  2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j+1)/kw(nx,j+1) )) + &!conduction to node r+dr
                  2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j-1)/kw(nx,j-1) )) + &!conduction to node r-dr
                  axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,j) + dx(nx-1)/kw(nx-1,j) ) ) !conduction to left node                      
    cols(ind+2) = 2*nx*nr_sf + j*nx
    
    !wall part, lower
    A( ind ) = -2 / ( r(jj+j)*dr(jj+j) / rdn(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j-1)/kw(nx,j-1) )) !conduction to node r-dr
    cols(ind) = cols(ind+2) - nx
    
    !wall part, left
    A( ind+1 ) = axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,j) + dx(nx-1)/kw(nx-1,j) ) ) !conduction to left node
    cols(ind+1) = cols(ind+2) - 1
    
    !wall part, upper
    A( ind+3 ) = -2 / ( r(jj+j)*dr(jj+j) / rup(jj+j) * ( dr(jj+j)/kw(nx,j) + dr(jj+j+1)/kw(nx,j+1) )) !conduction to node r+dr
    cols(ind+3) = cols(ind+2) + nx
    
enddo


!::Upper left corner, wall
!::Offset
!::Two lower corners (each 5)
!::Lower line (each 6)
!::Left and right lines (each 4)
!::Inner nodes (each 5)
ind = indFluidSolidOffSet + 2*5 + (nx-2)*6 +(nr_w-2)*2*4 + (nr_w-2)*(nx-2)*5 + 1
rows( 2*nx*nr_sf + nx*(nr_w-1)+1 ) = ind
jj = nr_sf 
!wall part, central node
A( ind+1 ) =   rhow( 1, nr_w ) * cw( 1, nr_w ) / dt + &!capacity term                    
                2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(1,nr_w) + dr(jj+nr_w-1)/kw(1,nr_w-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(1) * ( dx(1)/kw(1,nr_w) + dx(1+1)/kw(1+1,nr_w) ) ) !conduction to right node                    
cols(ind+1) = 2*nx*nr_sf + nx*(nr_w-1) + 1
!wall part, lower node
A( ind ) = &
    -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(1,nr_w) + dr(jj+nr_w-1)/kw(1,nr_w-1) ) ) !conduction to node r-dr
cols(ind) = cols(ind+1) - nx

!wall part, right node
A( ind+2 ) = &
    axial * -2 / ( dx(1) * ( dx(1)/kw(1,nr_w) + dx(1+1)/kw(1+1,nr_w) ) ) !conduction to right node
cols(ind+2) = cols(ind+1) + 1


!::Upper right corner, wall
!::Offset
!::Two lower corners (each 5)
!::Lower line (each 6)
!::Left and right lines (each 4)
!::Inner nodes (each 5)
!::Upper left corner (3)
!::Upper line (each 4)
ind = indFluidSolidOffSet + 2*5 + (nx-2)*6 +(nr_w-2)*2*4 + (nr_w-2)*(nx-2)*5 + 3 + (nx-2)*4 + 1
rows( 2*nx*nr_sf + nx*nr_w ) = ind

jj = nr_sf 
!wall part, central node
A( ind+2 ) =   rhow( nx, nr_w ) * cw( nx, nr_w ) / dt + &!capacity term
                2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(nx,nr_w) + dr(jj+nr_w-1)/kw(nx,nr_w-1) ) ) + &!conduction to node r-dr
                axial * 2 / ( dx(nx) * ( dx(nx)/kw(nx,nr_w) + dx(nx-1)/kw(nx-1,nr_w) ) ) !conduction to left node                    
cols(ind+2) = 2*nx*nr_sf + nx*nr_w
!wall part, lower node
A( ind ) = &
    -2 / ( r(jj+nr_w) * dr(jj+nr_w) / rdn(jj+nr_w) * ( dr(jj+nr_w)/kw(nx,nr_w) + dr(jj+nr_w-1)/kw(nx,nr_w-1) ) ) !conduction to node r-dr
cols(ind) = cols(ind+2) - nx
!wall part, left node
A( ind+1 ) = &
    axial * -2 / ( dx(nx) * ( dx(nx)/kw(nx,nr_w) + dx(nx-1)/kw(nx-1,nr_w) ) ) !conduction to left node
cols(ind+1) = cols(ind+2) - 1


!::End of the wall nodes::!
!::------------------------------::!
!::Begin of assembling the B vector (rhs of the equation system)

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
    do j=1,nr_sf
        B( (j-1)*nx + 1 ) = B( (j-1)*nx + 1 ) + (rhof(1,j) * cf(1,j) * u(1,j) * por * Tc / dx(1))        
    enddo
endif

if ( u(nx,1) .lt. 0 ) then
    !inlet from hot side (then an explicit term is added to the B vector)
    do j=1,nr_sf
        B( (j-1)*nx + nx ) = B( (j-1)*nx + nx ) - (rhof(nx,j) * cf(nx,j) * u(nx,j) * por * Th / dx(nx))
    enddo
endif
!call cpu_time( t2 )
!write(*,*) 'Bef matrix inversion'
!write(*,*) t2-t1

!open (12, file='A_vals.dat',	&
!           status='unknown', form='unformatted',	&
!           access='direct', recl=2*nNonZero)
!write(12,rec=1) A
!close(12)
!
!open (12, file='col_vals.dat',	&
!           status='unknown', form='unformatted',	&
!           access='direct', recl=1*nNonZero)
!write(12,rec=1) cols
!close(12)
!
!open (12, file='row_vals.dat',	&
!           status='unknown', form='unformatted',	&
!           access='direct', recl=1*(2*nx*nr_sf + nx*nr_w)+1)
!write(12,rec=1) rows
!close(12)

!::Solve Ax = B
call solve_lin_direct( A, B, X, cols,rows,nx*(2*nr_sf+nr_w), nNonZero )



!call cpu_time( t3 )
!write(*,*) 'after matrix inversion'
!write(*,*) t3-t2


do j=1,nr_sf
    Tfout(:,j) = X( nx*(j-1)+1:nx*j ) 
    Tsout(:,j) = X( (nx*nr_sf+nx*(j-1)+1):(nx*nr_sf+nx*j) ) 
enddo
do j=1,nr_w
    Twout(:,j) = X( (2*nx*nr_sf+nx*(j-1)+1):(2*nx*nr_sf+nx*j) ) 
enddo

deallocate( A,B, X,cols,rows )


end subroutine timestep_sp

end module SOLVER_CALL_SPARSE