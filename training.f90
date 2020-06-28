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

!ddrker30=0.0d0

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
real*8,allocatable,dimension(:)::vpot, coeff, energy, wt
real*8,allocatable, dimension(:,:)::dataarray, fdisp
real*8,parameter:: angtoau=1.0d0/0.529177249d0
integer, parameter :: na=4, nr=6, nker=10
real*8, parameter :: zero=0.0d0, & !zero is the sum of atomic energies
zero1=-16.085339160030447d0, ecut=10.5d0 !zero1 is the minimum energy
integer:: ndata, ngrad
end module

program potential
use arrays
implicit none
integer,parameter:: n=4, ndim=3*n, nframe=4001, ntrain=3200
real*8,dimension(:) :: r(nr)
real*8,allocatable, dimension(:,:)::disp, grad
integer::ii,jj,fno, id, kk, ll, idb, i1, i2, nl, i3
real*8 :: ener, gmax, gav, emax
character(len=10)::dummy

open(unit=10,file="training.dat")
open(unit=11,file="coeff.dat")

allocate(grad(ndim,nframe), disp(ndim,nframe), energy(nframe))

idb=0
do ii=1,nframe
  read(10,*)fno
  do jj=1,ndim,3
    read(10,*)dummy, disp(jj:jj+2,ii)
  end do

  read(10,*)energy(ii)

  if (energy(ii)-zero1 < ecut) idb=idb+1
  if (energy(ii)-zero1 > ecut) print*,idb,energy(ii),energy(ii)-zero1 ,ecut

  do jj=1,ndim,3
    read(10,*)grad(jj:jj+2,ii)
  end do

end do

print*,"idb =", idb
grad=-1.0d0*grad
!disp=disp*angtoau
gmax=max(maxval(grad),abs(minval(grad)))
emax=maxval(energy)-zero1

allocate(fdisp(ndim,idb), dataarray(idb,nr), vpot(idb*(ndim+1)), wt(idb*(ndim+1)))

i1=0
i2=0
do nl=1, nframe
  fno=nl
  id = 0
  r=0.0d0
  do ii = 1, na
    kk=(ii-1)*3+1
    do jj=ii+1,na
      ll = (jj-1)*3+1
      id = id+1
      r(id)=sqrt(sum((disp(kk:kk+2,fno)-disp(ll:ll+2,fno))**2))
    end do
  end do

  if (energy(fno)-zero1 < ecut) then
    i2=i2+1
    dataarray(i2,:)=r
    fdisp(:,i2)=disp(:,fno)
    i1=i1+1
    i3=i1
    gav=0.0d0
    vpot(i1)=energy(fno)-zero
    wt(i1)=emax/abs(emax+energy(fno)-zero1)
    do ii = 1, na
      i1=i1+1
      vpot(i1)=grad((ii-1)*3+1,fno)
      wt(i1)=gmax/(gmax+abs(vpot(i1)))!*4.0d0
      i1=i1+1
      vpot(i1)=grad((ii-1)*3+2,fno)
      wt(i1)=gmax/(gmax+abs(vpot(i1)))!*4.0d0
      i1=i1+1
      vpot(i1)=grad((ii-1)*3+3,fno)
      wt(i1)=gmax/(gmax+abs(vpot(i1)))!*4.0d0
    end do
  end if
  if (i2==ntrain) exit
end do

ndata=i2
print*,ndata

ngrad=i1-ndata
print*,ndata,ngrad,i1

deallocate(disp,grad,energy)

call trik()

!Test
do ii = 1, ndata
  r=dataarray(ii,:)
  call eval(r,ener)
  write(91,*)vpot((ii-1)*(ndim+1)+1), ener
end do

end program

subroutine eval(r1,ener)
use arrays
use rep_ker
implicit none
real*8, dimension (:) :: r1(nr), r2(nr)
real*8, dimension (:,:) :: kerval(3,nker)
real*8 :: ener
integer :: ii

ener=0.0d0
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
      
end do

end subroutine 

subroutine trik()
use arrays
use rep_ker
implicit none
real*8, allocatable, dimension(:,:) :: qm, dqm
integer :: ii, jj, kk, id, id1, id2, k1
integer :: n, nrhs, rank
integer :: lda, ldb, info
real*8, parameter :: piby180 = acos(-1.0d0)/180.0d0
real*8, dimension (:) :: r1(nr), r2(nr), der(nr)
real*8, dimension (:,:) :: kerval(3,nker), derval(3,nker)
real*8, dimension (:) :: tmpx(na), tmpy(na), tmpz(na), tmpdx(na), tmpdy(na), tmpdz(na)
real*8, allocatable, dimension(:) :: bi, s
real*8, dimension(:) :: work(150000)

allocate(qm(ndata+ngrad,ndata), coeff(ndata),s(ndata), bi(ndata+ngrad), dqm(ndata+ngrad,ndata))

n = ndata+ngrad
nrhs = 1
lda = ndata
ldb = ndata

print*,ndata,ngrad, n

qm=0.0d0

k1=0
do ii = 1, ndata
  r1(:)=dataarray(ii,:)
  do jj = 1, ndata
    kk=k1
    r2(:)=dataarray(jj,:)

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

kk=kk+1

qm(kk, jj) = &
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
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10)

der( 1 )=&
2.0d0*derval(1,1)+&
derval(2,1)*kerval(2,2)*kerval(2,6)+&
derval(2,1)*kerval(2,3)*kerval(2,7)+&
derval(2,1)*kerval(2,4)*kerval(2,8)+&
derval(2,1)*kerval(2,5)*kerval(2,9)+&
derval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
derval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10)

der( 2 )=&
derval(1,2)+&
derval(1,3)+&
kerval(2,1)*derval(2,2)*kerval(2,6)+&
kerval(2,1)*derval(2,3)*kerval(2,7)+&
derval(2,2)*kerval(2,5)*kerval(2,10)+&
derval(2,3)*kerval(2,4)*kerval(2,10)+&
kerval(3,1)*derval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*derval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10)

der( 3 )=&
derval(1,4)+&
derval(1,5)+&
kerval(2,1)*derval(2,4)*kerval(2,8)+&
kerval(2,1)*derval(2,5)*kerval(2,9)+&
kerval(2,2)*derval(2,5)*kerval(2,10)+&
kerval(2,3)*derval(2,4)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*derval(3,5)*kerval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*derval(3,4)*kerval(3,7)*kerval(3,8)*kerval(3,10)

der( 4 )=&
derval(1,6)+&
derval(1,7)+&
kerval(2,1)*kerval(2,2)*derval(2,6)+&
kerval(2,1)*kerval(2,3)*derval(2,7)+&
derval(2,6)*kerval(2,9)*kerval(2,10)+&
derval(2,7)*kerval(2,8)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*derval(3,6)*kerval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*derval(3,7)*kerval(3,8)*kerval(3,10)

der( 5 )=&
derval(1,8)+&
derval(1,9)+&
kerval(2,1)*kerval(2,4)*derval(2,8)+&
kerval(2,1)*kerval(2,5)*derval(2,9)+&
kerval(2,6)*derval(2,9)*kerval(2,10)+&
kerval(2,7)*derval(2,8)*kerval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*derval(3,9)*kerval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*derval(3,8)*kerval(3,10)

der( 6 )=&
2.0d0*derval(1,10)+&
kerval(2,2)*kerval(2,5)*derval(2,10)+&
kerval(2,3)*kerval(2,4)*derval(2,10)+&
kerval(2,6)*kerval(2,9)*derval(2,10)+&
kerval(2,7)*kerval(2,8)*derval(2,10)+&
kerval(3,1)*kerval(3,2)*kerval(3,5)*kerval(3,6)*kerval(3,9)*derval(3,10)+&
kerval(3,1)*kerval(3,3)*kerval(3,4)*kerval(3,7)*kerval(3,8)*derval(3,10)
      
  do id1 = 1, na
    tmpx(id1)=fdisp((id1-1)*3+1,ii)!*angtoau
    tmpy(id1)=fdisp((id1-1)*3+2,ii)!*angtoau
    tmpz(id1)=fdisp((id1-1)*3+3,ii)!*angtoau
  end do
  tmpdx=0.0d0
  tmpdy=0.0d0
  tmpdz=0.0d0
  id=0
  do id1 = 1, na
    do id2=id1+1,na
      id=id+1
      tmpdx(id1) = tmpdx(id1) + (tmpx(id1)-tmpx(id2))/r1(id)*der(id)
      tmpdx(id2) = tmpdx(id2) + (tmpx(id2)-tmpx(id1))/r1(id)*der(id)

      tmpdy(id1) = tmpdy(id1) + (tmpy(id1)-tmpy(id2))/r1(id)*der(id)
      tmpdy(id2) = tmpdy(id2) + (tmpy(id2)-tmpy(id1))/r1(id)*der(id)

      tmpdz(id1) = tmpdz(id1) + (tmpz(id1)-tmpz(id2))/r1(id)*der(id)
      tmpdz(id2) = tmpdz(id2) + (tmpz(id2)-tmpz(id1))/r1(id)*der(id)
    end do
  end do

  do id1=1,na
    kk=kk+1
    qm(kk, jj)=tmpdx(id1)
    kk=kk+1
    qm(kk, jj)=tmpdy(id1)
    kk=kk+1
    qm(kk, jj)=tmpdz(id1)
  end do

  end do
  k1=kk

end do
print*,kk,k1

do ii = 1, ndata+ngrad
  qm(ii,:) = qm(ii,:)*wt(ii)
  bi(ii)=vpot(ii)*wt(ii)
end do

dqm=qm
print*,"TOT QMATRIX DONE"
call system('date')
!--------------lapack---------------
!   solve the equations a*x = b.
call dgelss(n,ndata,1,qm,n,bi,n,s,1.0d-15,rank,work,150000,info)
coeff(:)=bi(1:ndata)

!     check for the exact singularity.
if( info.gt.0 ) then
   write(*,*)'the leading minor of order ',info,' is not positive'
   write(*,*)'definite; the solution could not be computed.'
   stop
end if
!------------------------------------
print*,"TOT COEFF DONE"

do ii=1,ndata
  write(11,*)dataarray(ii,:),coeff(ii)
end do

do ii=1,n
  write(92,*)vpot(ii)*wt(ii),dot_product(dqm(ii,:),coeff(:)), wt(ii)
end do

call system('date')
deallocate(qm,dqm)

end subroutine
