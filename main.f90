Program main
implicit none
integer,parameter  :: dp=8
!real(dp),parameter :: beta=6._dp
!integer,parameter :: ibeta=2
integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,j
integer :: b1,b2,b3,b4,b11,b12,b13,b14
integer :: c1,c2,c3,c4,c11,c12,c13,c14
integer :: j1,j2,jj1,jj2,jjj1,jjj2

real(dp):: z5(2,2,2), z6(2,2,2,2),z78(8,2,8,2)
real(dp)::twopatch(2,8,2,2,8,2)
real(dp)::tripatch(8,2,2,2,8,2)
real(dp)::e_H

!real(dp),parameter:: HJ(2,2)=reshape(  &
! (/ exp(-beta),exp(beta) ,&
!    exp(beta),exp(-beta) /) , (/2,2/) )

!real(dp),parameter:: HJ(2,2)=reshape(  &
! (/ 0.1_dp**ibeta,10._dp**ibeta ,&
!    10._dp**ibeta,0.1_dp**ibeta /) , (/2,2/) )

integer  :: seq
real(dp) :: beta, HJ(2,2)

 write(6,*) "       J","   Z=Sum e^-H"
do seq=1,10,9
beta=real(seq,dp)
HJ(1,1)=exp(-beta)
HJ(2,2)=exp(-beta)
HJ(1,2)=exp(beta)
HJ(2,1)=exp(beta)
!write(6,"(8F12.6)") HJ(:,:)

z5(:,:,:)=0._dp
do i1=1,2
do i6=1,2
do i11=1,2
   z5(i11,i6,i1)=fz5(i11,i6,i1)
enddo
enddo
enddo
!write(6,*) "z5"
!write(6,"(8F12.6)") z5(:,:,:)

z6 (:,:,:,:)=0._dp
do i1=1,2
do i7=1,2
do i14=1,2
do i11=1,2
   z6(i11,i14,i7,i1)=z5(i11,1,i1)*HJ(1,i14)*HJ(1,i7) &
                    +z5(i11,2,i1)*HJ(2,i14)*HJ(2,i7)
enddo
enddo
enddo
enddo
!write(6,*) "z6"
!write(6,"(8F12.6)") z6(:,:,:,:)


z78(:,:,:,:)=0._dp
do i1=1,2
j1=0
do i2=1,2
do i3=1,2
do i4=1,2
j1=j1+1
do i11=1,2
j2=0
do i12=1,2
do i13=1,2
do i14=1,2
j2=j2+1
   do i7=1,2
   do i8=1,2
   z78(j2,i11,j1,i1)=z78(j2,i11,j1,i1) &
     +z6(i11,i14,i7,i1)*HJ(i2,i1)*HJ(I7,I2) &
     *z6(i12,i13,i8,i4)*HJ(i3,i4)*HJ(I8,I3) &
     *HJ(I8,I7)
   enddo
   enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
!write(6,*) "z78"
!write(6,"(8E12.4)") z78(:,:,:,:)

twopatch(:,:,:,:,:,:)=0._dp
do i1=1,2
j1=0;do i2=1,2;do i3=1,2;do i4=1,2;j1=j1+1

  do i11=1,2
! j2=0;do i12=1,2;do i13=1,2;do i14=1,2;j2=j2+1

do b1=1,2
!jj1=0;do b2=1,2;do b3=1,2;do b4=1,2;jj1=jj1+1

  do b11=1,2
  jj2=0;do b12=1,2;do b13=1,2;do b14=1,2;jj2=jj2+1

  j2=0;do i12=1,2;do i13=1,2;do i14=1,2;j2=j2+1

  twopatch(b11,jj2,b1,i11,j1,i1)=twopatch(b11,jj2,b1,i11,j1,i1) &
    +z78(j2,i11,j1,i1)*z78(jj2,b11,j2,b1) 
   !+z78(j2,i11,j1,i1)*z78(jj2,b11,jj1,b1) !(jj1=j2)
  enddo
  enddo
  enddo

enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
!write(6,*) "twopatch"
!write(6,"(8E12.4)") twopatch(:,:,:,:,:,:)


 tripatch(:,:,:,:,:,:)=0._dp
 do i1=1,2
 j1=0;do i2=1,2;do i3=1,2;do i4=1,2;j1=j1+1
 do i11=1,2 !loop
 !j2=0;do i12=1,2;do i13=1,2;do i14=1,2;j2=j2+1
 
 do b1=1,2
 !jj1=0;do b2=1,2;do b3=1,2;do b4=1,2;jj1=jj1+1
 do b11=1,2
 !jj2=0;do b12=1,2;do b13=1,2;do b14=1,2;jj2=jj2+1
 
 
  !do c1=1,2
   jjj1=0;do c2=1,2;do c3=1,2;do c4=1,2;jjj1=jjj1+1 !loop
   do c11=1,2
   jjj2=0;do c12=1,2;do c13=1,2;do c14=1,2;jjj2=jjj2+1
 
 ! print*,"wzj"
 
   tripatch(jjj2,c11,b11,b1,j1,i1)= tripatch(jjj2,c11,b11,b1,j1,i1) &
    +twopatch(b11,jjj1,b1,i11,j1,i1)*z78(jjj2,c11,jjj1,i11)
!!  +twopatch(b11,jj2 ,b1,i11,j1,i1)*z78(jjj2,c11,jjj1,c1)
!!  i11=c1; jj2=jjj1
!
  !enddo
  !enddo
  !enddo
  !enddo
   enddo
   enddo
   enddo
   enddo
 
 
 enddo
 enddo
 enddo
 enddo
 enddo
     
 enddo
 enddo
 enddo
 enddo
 enddo
!write(6,*) "tripatch"
!write(6,"(8E12.4)") tripatch(:,:,:,:,:,:)




e_H=0._dp
 do i1=1,2
 j1=0;do i2=1,2;do i3=1,2;do i4=1,2;j1=j1+1
 
 do b1=1,2
 do b11=1,2

 do c11=1,2
 jjj2=0;do c12=1,2;do c13=1,2;do c14=1,2;jjj2=jjj2+1
 
 !print*,"wzj"
 e_H=e_H+tripatch(jjj2,c11,b11,b1,j1,i1)*tripatch(j1,b1,i1,c11,jjj2,b11)
 enddo
 enddo
 enddo
 enddo
 enddo
     
 enddo
 enddo
 enddo
 enddo
 enddo
 write(6,"(F12.6,E16.6,3F16.8)") beta,e_H,log(e_H)/60._dp, e_H*exp(-66._dp*beta),-(log(e_H)-log(1.6d4))/beta
enddo
 write(6,*) "                               Degeneracy","   lowest Energy"

contains
 real(dp) function  fz5(i,j,k)
 integer :: i,j,k
 fz5=0._dp
 fz5= HJ(i,1)*HJ(j,1)*HJ(k,1) + HJ(i,2)*HJ(j,2)*HJ(k,2)
 end function
end Program main
