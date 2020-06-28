module rep_ker
implicit none

contains

function drker30(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker30, xl, xs

xl = x
xs = xi
if (x .le. xi) then
  xl = xi
  xs = x
end if

drker30=3d0/xl - 3d0/2d0 * xs/xl**2 + 3d0/10d0 * xs**2/xl**3

end function drker30

function ddrker30(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker30, xl, xs

if (x .le. xi) then
  ddrker30 =  3d0/5d0 * x/xi**3 - 3d0/(2d0*xi**2)
else
  ddrker30  = -3d0/x**2 + 3d0*xi/x**3 - 9d0/10d0 * xi**2/x**4
end if

end function ddrker30

function drker31(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker31, xl, xs

xl = x
xs = xi
if (x .le. xi) then
  xl = xi
  xs = x
end if

drker31 = 3d0/(4d0*xl**2) - 3d0/5d0 * xs/xl**3 + 3d0/20d0 * xs**2/xl**4
end function drker31

function ddrker31(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker31, xl, xs

if (x .le. xi) then
  ddrker31 = 3d0/10d0 * x/xi**4 - 3d0/(5d0*xi**3)
else
  ddrker31  = -3d0/(2d0*x**3) + 9d0/5d0 * xi/x**4 - 3d0/5d0 * xi**2/x**5
end if

!ddrker31 = 0.0d0
end function ddrker31

function drker35(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker35, xl, xs

xl = x
xs = xi
if (x .le. xi) then
  xl = xi
  xs = x
end if

drker35 = 3d0/(56d0*xl**6) - 1d0/14d0 * xs/xl**7 + 1d0/40d0 * xs**2/xl**8

end function drker35


function ddrker35(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker35, xl, xs

if (x .le. xi) then
  ddrker35 = 1d0/20d0 * x/xi**8 - 1d0/(14d0*xi**7)
else
  ddrker35  = -1d0/5d0 * xi**2/x**9 + 0.5d0 * xi/x**8 - 9d0/(28d0*x**7)
end if
end function ddrker35

end module

module arrays
implicit none
real*8,allocatable,dimension(:):: coeff
real*8,allocatable, dimension(:,:)::dataarray
real*8,parameter:: angtoau=1.0d0/0.529177249d0
integer, parameter :: na=4, nr=6, nker=10
real*8, parameter :: zero=0.0d0, &
zero1=-16.085339160030447d0,ecut=10.5d0
integer,parameter:: nframe=2400, ndim=3*na, ndata=3200
end module

subroutine eval(r1,ener,der)
use arrays
use rep_ker
implicit none
real*8, dimension (:) :: r1(nr), r2(nr), der(nr)
real*8, dimension (:) :: kerval(3,nker), derval(3,nker)
real*8 :: ener
integer :: ii

ener=0.0d0
der=0.0d0
do ii = 1, ndata
  r2(:)=dataarray(ii,:)
kerval(1,1)=drker35(r1(1),r2(1))
kerval(1,2)=drker35(r1(2),r2(2))
kerval(1,3)=drker35(r1(2),r2(3))
kerval(1,4)=drker35(r1(3),r2(2))
kerval(1,5)=drker35(r1(3),r2(3))
kerval(1,6)=drker35(r1(4),r2(4))
kerval(1,7)=drker35(r1(4),r2(5))
kerval(1,8)=drker35(r1(5),r2(4))
kerval(1,9)=drker35(r1(5),r2(5))
kerval(1,10)=drker35(r1(6),r2(6))

kerval(2,1)=drker31(r1(1),r2(1))
kerval(2,2)=drker31(r1(2),r2(2))
kerval(2,3)=drker31(r1(2),r2(3))
kerval(2,4)=drker31(r1(3),r2(2))
kerval(2,5)=drker31(r1(3),r2(3))
kerval(2,6)=drker31(r1(4),r2(4))
kerval(2,7)=drker31(r1(4),r2(5))
kerval(2,8)=drker31(r1(5),r2(4))
kerval(2,9)=drker31(r1(5),r2(5))
kerval(2,10)=drker31(r1(6),r2(6))

kerval(3,1)=drker30(r1(1),r2(1))
kerval(3,2)=drker30(r1(2),r2(2))
kerval(3,3)=drker30(r1(2),r2(3))
kerval(3,4)=drker30(r1(3),r2(2))
kerval(3,5)=drker30(r1(3),r2(3))
kerval(3,6)=drker30(r1(4),r2(4))
kerval(3,7)=drker30(r1(4),r2(5))
kerval(3,8)=drker30(r1(5),r2(4))
kerval(3,9)=drker30(r1(5),r2(5))
kerval(3,10)=drker30(r1(6),r2(6))

derval(1,1)=ddrker35(r1(1),r2(1))
derval(1,2)=ddrker35(r1(2),r2(2))
derval(1,3)=ddrker35(r1(2),r2(3))
derval(1,4)=ddrker35(r1(3),r2(2))
derval(1,5)=ddrker35(r1(3),r2(3))
derval(1,6)=ddrker35(r1(4),r2(4))
derval(1,7)=ddrker35(r1(4),r2(5))
derval(1,8)=ddrker35(r1(5),r2(4))
derval(1,9)=ddrker35(r1(5),r2(5))
derval(1,10)=ddrker35(r1(6),r2(6))

derval(2,1)=ddrker31(r1(1),r2(1))
derval(2,2)=ddrker31(r1(2),r2(2))
derval(2,3)=ddrker31(r1(2),r2(3))
derval(2,4)=ddrker31(r1(3),r2(2))
derval(2,5)=ddrker31(r1(3),r2(3))
derval(2,6)=ddrker31(r1(4),r2(4))
derval(2,7)=ddrker31(r1(4),r2(5))
derval(2,8)=ddrker31(r1(5),r2(4))
derval(2,9)=ddrker31(r1(5),r2(5))
derval(2,10)=ddrker31(r1(6),r2(6))

derval(3,1)=ddrker30(r1(1),r2(1))
derval(3,2)=ddrker30(r1(2),r2(2))
derval(3,3)=ddrker30(r1(2),r2(3))
derval(3,4)=ddrker30(r1(3),r2(2))
derval(3,5)=ddrker30(r1(3),r2(3))
derval(3,6)=ddrker30(r1(4),r2(4))
derval(3,7)=ddrker30(r1(4),r2(5))
derval(3,8)=ddrker30(r1(5),r2(4))
derval(3,9)=ddrker30(r1(5),r2(5))
derval(3,10)=ddrker30(r1(6),r2(6))

ener = ener + (&
2.0d0*kerval(1,1)+&
kerval(1,2)+&
kerval(1,3)+&
kerval(1,4)+&
kerval(1,5)+&
kerval(1,6)+&
kerval(1,7)+&
kerval(1,8)+&
kerval(1,9)+&
2.0d0*kerval(1,10)+&
kerval(2,1)*kerval(2,2)*kerval(2,6)+&
kerval(2,1)*kerval(2,3)*kerval(2,7)+&
kerval(2,1)*kerval(2,4)*kerval(2,8)+&
kerval(2,1)*kerval(2,5)*kerval(2,9)+&
kerval(2,2)*kerval(2,5)*kerval(2,10)+&
kerval(2,3)*kerval(2,4)*kerval(2,10)+&
kerval(2,6)*kerval(2,9)*kerval(2,10)+&
kerval(2,7)*kerval(2,8)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10))*coeff(ii)
      
der( 1 )=der(1) + (&
2.0d0*derval(1,1)+&
derval(2,1)*kerval(2,2)*kerval(2,6)+&
derval(2,1)*kerval(2,3)*kerval(2,7)+&
derval(2,1)*kerval(2,4)*kerval(2,8)+&
derval(2,1)*kerval(2,5)*kerval(2,9)+&
derval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
derval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10))*coeff(ii)

der( 2 )= der(2) + (&
derval(1,2)+&
derval(1,3)+&
kerval(2,1)*derval(2,2)*kerval(2,6)+&
kerval(2,1)*derval(2,3)*kerval(2,7)+&
derval(2,2)*kerval(2,5)*kerval(2,10)+&
derval(2,3)*kerval(2,4)*kerval(2,10)+&
kerval(3,1)*derval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*derval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10))*coeff(ii)

der( 3 )= der(3) + (&
derval(1,4)+&
derval(1,5)+&
kerval(2,1)*derval(2,4)*kerval(2,8)+&
kerval(2,1)*derval(2,5)*kerval(2,9)+&
kerval(2,2)*derval(2,5)*kerval(2,10)+&
kerval(2,3)*derval(2,4)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*derval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*derval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10))*coeff(ii)

der( 4 )=der(4) + (&
derval(1,6)+&
derval(1,7)+&
kerval(2,1)*kerval(2,2)*derval(2,6)+&
kerval(2,1)*kerval(2,3)*derval(2,7)+&
derval(2,6)*kerval(2,9)*kerval(2,10)+&
derval(2,7)*kerval(2,8)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*derval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*derval(3,7)*kerval(3,8)*kerval(3,10))*coeff(ii)

der( 5 )=der(5) + (&
derval(1,8)+&
derval(1,9)+&
kerval(2,1)*kerval(2,4)*derval(2,8)+&
kerval(2,1)*kerval(2,5)*derval(2,9)+&
kerval(2,6)*derval(2,9)*kerval(2,10)+&
kerval(2,7)*derval(2,8)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*derval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*derval(3,8)*kerval(3,10))*coeff(ii)

der( 6 )=der(6) + (&
2.0d0*derval(1,10)+&
kerval(2,2)*kerval(2,5)*derval(2,10)+&
kerval(2,3)*kerval(2,4)*derval(2,10)+&
kerval(2,6)*kerval(2,9)*derval(2,10)+&
kerval(2,7)*kerval(2,8)*derval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*derval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*derval(3,10))*coeff(ii)

end do

end subroutine

program peswrapper
use arrays
implicit none
real*8,parameter::pi=dacos(-1.0d0),const=pi/180.0d0, dx=0.005d0
real*8:: ener, d0h, d02h
real*8,dimension(:)::dvdr(nr), rdist(nr), r(nr), dvdx(4)
real*8,dimension(:,:)::xyz(na,3), dvdxyz(na,3), txyz(na,3)
character (len=6) :: symbol
integer :: ii, id, i1, i2

allocate(dataarray(ndata,nr),coeff(ndata))
open(unit=11,file="coeff.dat")
do ii=1,ndata
  read(11,*)dataarray(ii,:),coeff(ii)
end do

open(unit=10,file="inp.xyz", status="old")
read(10,*) 
read(10,*)
do ii = 1, na 
  read(10,*) symbol, xyz(ii,:)
end do

!xyz=xyz/0.529177249d0
rdist=0.0d0
id=0
do i1 = 1, na
  do i2=i1+1,na
    id = id+1
    rdist(id)=sqrt((xyz(i1,1)-xyz(i2,1))**2+(xyz(i1,2)-xyz(i2,2))**2+&
    (xyz(i1,3)-xyz(i2,3))**2)
  end do
end do

call eval(rdist,ener,dvdr)

dvdxyz=0.0d0
id=0
do i1 = 1, na
  do i2=i1+1,na
    id=id+1
    dvdxyz(i1,1) = dvdxyz(i1,1) + (xyz(i1,1)-xyz(i2,1))/rdist(id)*dvdr(id)
    dvdxyz(i2,1) = dvdxyz(i2,1) + (xyz(i2,1)-xyz(i1,1))/rdist(id)*dvdr(id)

    dvdxyz(i1,2) = dvdxyz(i1,2) + (xyz(i1,2)-xyz(i2,2))/rdist(id)*dvdr(id)
    dvdxyz(i2,2) = dvdxyz(i2,2) + (xyz(i2,2)-xyz(i1,2))/rdist(id)*dvdr(id)

    dvdxyz(i1,3) = dvdxyz(i1,3) + (xyz(i1,3)-xyz(i2,3))/rdist(id)*dvdr(id)
    dvdxyz(i2,3) = dvdxyz(i2,3) + (xyz(i2,3)-xyz(i1,3))/rdist(id)*dvdr(id)
  end do
end do

!do ii=1,na
!  do i2=1,3
!  txyz=xyz
!  txyz(ii,i2)=xyz(ii,i2)-2.0d0*dx
!  call rcalc(txyz,r)
!  call eval(r,ener,dvdr)
!  dvdx(1)=ener
!
!  txyz=xyz
!  txyz(ii,i2)=xyz(ii,i2)-dx
!  call rcalc(txyz,r)
!  call eval(r,ener,dvdr)
!  dvdx(2)=ener
!
!  txyz=xyz
!  txyz(ii,i2)=xyz(ii,i2)+dx
!  call rcalc(txyz,r)
!  call eval(r,ener,dvdr)
!  dvdx(3)=ener
!
!  txyz=xyz
!  txyz(ii,i2)=xyz(ii,i2)+2.0d0*dx
!  call rcalc(txyz,r)
!  call eval(r,ener,dvdr)
!  dvdx(4)=ener
!
!  d0h=(dvdx(3)-dvdx(2))/2.0d0/dx
!  d02h=(dvdx(4)-dvdx(1))/4.0d0/dx
!  dvdxyz(ii,i2)=(4.0d0*d0h-d02h)/3.0d0
!  end do
!end do
!
!call eval(rdist,ener,dvdr)

!ener=ener*27.21138d0

!dvdxyz=dvdxyz/0.529177249d0*27.21138d0

open(unit=11,file="ener.out")
open(unit=12,file="grad.out")

write(11,*)"FINAL SINGLE POINT ENERGY =", ener, "eV (zero is set to the sum of atomic enetgies)"
write(12,*)"# The current gradient in eV/angstrom"
do ii = 1, na
  write(12,*)dvdxyz(ii,:)
end do
write(12,*)"# The end"

close(11)
close(12)

end program
