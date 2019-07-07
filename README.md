# truss

program truss

implicit none
real :: L0,P,Pe,fs,fl,fc,ft,ffp,ffc,ffs
real :: updown,weight_min,rax,ray,rlx,lem
real :: L(4,6),N(4,6),rivet(17,2),element(31,6),An(17,31),lambda(17,31),omega(17,31),fk(17,31),diameter(4,6),anc(9,6)
real :: Rs1(17),Rs2(17),Rl(17),tmin(17),weight(17,31),e1(17),e2(17),rivet_pitch(17)
real :: lk(17),g1(17),g2(17),length(4,6),width(4,6),thick(6,6)
integer :: number(17),num(4,6),a_number(2)
real :: alpha,pi
integer :: i,j,k,m

fs = 120D0 !MPa
fl = 300D0 !MPa
fc = 160D0 !MPa
ft = 160D0 !MPa
ffp = 160D0 !MPa
ffc = 160D0 !MPa
ffs = 120D0 !MPa

write(*,*) "Enter l"
read(*,*) L0
write(*,*) "Enter p"
read(*,*) P
Pe = P*1.5 !安全率を考慮

alpha = ATAN(1/L0)
pi = 4 * ATAN(1.0)

write(*,*) "right_up(1) or right_down(2)?"
read(*,*) updown

if(updown == 1) then
call right_up(L0,Pe,N,L,alpha)
else if(updown == 2) then
call right_down(L0,Pe,N,L,alpha)
else
write(*,*) "error"
end if

open(10,file="rivet.dat")
do i = 1,17
read(10,*) rivet(i,1:2)
e1(i) = 5D0 / 3D0 * (rivet(i,1) + 5D0)
e2(i) = 5D0 / 3D0 * (rivet(i,1) - 1D0)
rivet_pitch(i) = 2.5 * rivet(i,1)
end do
close(10)

do i = 1,17
Rs1(i) = fs * pi * (rivet(i,2)**2) / 4D0 !一面剪断
Rs2(i) = fs * 2D0 * pi * (rivet(i,2)**2) / 4D0 !二面剪断
end do

open(11,file="element.dat")
do i=1,31
read(11,*) element(i,1:6)
end do
close(11)

!Deside number of rivet
do i = 1,4
do j = 1,6

weight_min = 10000000

if(N(i,j) > 0) then !引張

do k = 1,17

if( i == 1 ) then
number(k) = ceiling(abs(N(i,j)*1000)/Rs2(k)) !二面剪断
else if( i == 4 ) then
number(k) = ceiling(abs(N(i,j)*1000)/Rs2(k)) !二面剪断
else
number(k) = ceiling(abs(N(i,j)*1000)/Rs1(k)) !一面剪断
end if

g1(k) = (e1(k) + (rivet_pitch(k) * (number(k) - 1))) * 2D0 !ガセットプレートが覆う長さ
!g2(k) = (e1(k) + (rivet_pitch(k) * (ceiling((number(k) / 2D0)) - 1))) * 2
!lk(k) = L(i,j) - 2D0 * (e2(k) + rivet_pitch(k) * (number(k) - 1)) !(部材の長さ)-(ガセットプレートが覆う長さ)
Rl(k) = N(i,j)/dble(number(k)) !リベット1本あたりの剪断力
tmin(k) = Rl(k) / (fl * rivet(i,2)) !最小厚さ
do m = 1,31
An(k,m) = element(m,3) - rivet(k,2) * element(m,2) !有効断面積
if(N(i,j)/An(k,m) <= ft & !有効断面積が満たすべき条件
 .and. element(m,2) >= tmin(k) & !部材の厚さが満たすべき条件
 .and. rivet(k,1) <= element(m,6) & !リベット径が満たすべき条件 
 .and. g1(k) < (L(i,j) * 0.21 * 1000D0) & !ガセットプレートが覆う長さが満たす範囲
 .and. number(k) <= 6  & !リベットの本数が満たすべき条件
 .and. 2 * element(m,2) <= 10 * rivet(k,2)) then !板の総厚が満たすべき条件

weight(k,m) = element(m,5) * L(i,j) !部材の重量

if(weight(k,m) < weight_min) then
weight_min = weight(k,m)
width(i,j) = element(m,1)
thick(i,j) = element(m,2)
num(i,j) = number(k)
diameter(i,j) = rivet(k,1)
length(i,j) = g1(k)
end if

end if

end do 
end do

end if

if(N(i,j) < 0D0) then !圧縮
do k = 1,17

if( i == 1 ) then
number(k) = ceiling(abs(N(i,j)*1000)/Rs2(k))
else if( i == 4 ) then
number(k) = ceiling(abs(N(i,j)*1000)/Rs2(k))
else
number(k) = ceiling(abs(N(i,j)*1000)/Rs1(k))
end if

g1(k) = (e1(k) + (rivet_pitch(k) * (number(k) - 1))) * 2
g2(k) = (e1(k) + (rivet_pitch(k) * (ceiling((number(k) / 2D0)) - 1))) * 2
lk(k) = L(i,j) - 2D0 * (e2(k) + rivet_pitch(k) * (number(k) - 1))
Rl(k) = abs(N(i,j))/dble(number(k))
tmin(k) = Rl(k)/rivet(i,2)
do m = 1,31

lambda(k,m) = lk(k) / element(m,4)

!座屈応力を求める
if(lambda(k,m) < 30) then
fk(k,m) = fc
else if(30 <= lambda(k,m) .and. lambda(k,m) <= 100) then
fk(k,m) = fc * (1 - 0.400 * (lambda(k,m) / 100)**2)
else if(lambda(k,m) > 100) then
fk(k,m) = 0.600 * fc / ((lambda(k,m) / 100)**2)
else
fk(k,m) = 0D0
end if

omega(k,m) = fc / fk(k,m)

if(abs(N(i,j)) * omega(k,m) / element(m,3) <= fc & !座屈しないための条件
 .and. element(m,2) >= tmin(k) & !部材の厚さが満たすべき条件
 .and. rivet(k,1) <= element(m,6) & !リベット径が満たすべき条件 
 .and. g1(k) < (L(i,j) * 0.2 * 1000) & !ガセットプレートが覆う長さが満たす範囲
 .and. number(k) <= 6 & !リベットの本数が満たすべき条件
 .and. 2 * element(m,2) <= 10 * rivet(k,2)) then !板の総厚が満たすべき条件

weight(k,m) = element(m,5) * L(i,j) !部材の重量

if(weight(k,m) < weight_min) then
weight_min = weight(k,m)
width(i,j) = element(m,1)
thick(i,j) = element(m,2) 
num = number(k)
diameter(i,j) = rivet(k,1)
length(i,j) = g1(k)
end if

end if

end do 
end do

end if

if(j < 6) then
write(*,*) "L(",i,",",j,") = ", L(i,j) * 1000
write(*,*) "N(",i,",",j,") = ", N(i,j)
write(*,*) "Number of rivet = ",num(i,j)
write(*,*) "Diameter of rivet = ",diameter(i,j)
write(*,*) "Length covered by gusset plate : ",length(i,j)
write(*,*) "Width of element : ",width(i,j)
write(*,*) "Thickness of element",thick(i,j)
write(*,*) "Percentage of length",length(i,j) / (L(i,j)*10) 
else if(i == 2) then
write(*,*) "L(",i,",",j,") = ", L(i,j) * 1000
write(*,*) "N(",i,",",j,") = ", N(i,j)
write(*,*) "Number of rivet = '",num(i,j)
write(*,*) "Diameter of rivet = ",diameter(i,j)
write(*,*) "Length covered by gusset plate : ",length(i,j)
write(*,*) "Width of element : ",width(i,j)
write(*,*) "Thickness of element",thick(i,j)
write(*,*) "Percentage of length",length(i,j) / (L(i,j)*10) 
end if

end do
end do

!ガセットプレート厚さ
do i = 1,6
thick(5,i) = max(thick(1,i),thick(1,i+1),thick(2,i),thick(3,i))
thick(6,i) = max(thick(2,i),thick(3,i-1),thick(4,i-1),thick(4,i))
if(thick(5,i) < 3.2D0) then
thick(5,i) = 3.2D0
end if
if(thick(6,i) < 3.2D0) then
thick(6,i) = 3.2D0
end if
if(i < 6) then
write(*,*) "Thickness of gusset plate (1,",i,") = ", thick(5,i)
end if
write(*,*) "Thickness of gusset plate (2,",i,") = ", thick(6,i)
end do

!反力
rax = -L0*Pe/4.0
ray = Pe/2
rlx = -rax

write(*,*) "Power of anchor_bolt"
write(*,*) "Rax = ",rax,"Ray = ",ray,"Rlx = ",rlx

!アンカーボルト
open(12,file="anchor_bolt.dat")
do i = 1,9
read(12,*) anc(i,1:4)
anc(i,5) = (anc(i,1) ** 2) * pi / 4D0 !断面積
anc(i,6) = anc(i,1) * pi !円周
if(ceiling(abs(rax) * 1000 /(anc(i,5) * ffp)) > ceiling(abs(ray) * 1000/(anc(i,5) * ffs))) then
a_number(1) = ceiling(abs(rax) * 1000 /(anc(i,5) * ffp))
else 
a_number(1) = ceiling(abs(ray) * 1000 /(anc(i,5) * ffs))
end if
a_number(2) = ceiling(abs(rlx) * 1000 /(anc(i,5) * ffp))
lem = (abs(rax) * 1000 /(0.7D0 * a_number(1))-((pi*pi) * anc(i,1) * (anc(i,4) + anc(i,1) / 2D0))) / anc(i,6) &
- anc(i,3) + anc(i,2) + (2 * (anc(i,1) + anc(i,4))) !埋め込み長さ

write(*,*) "Number of anchor bolt",a_number(1),a_number(2)
write(*,*) "Length of anchor bolt",lem

end do
close(12)

end program

subroutine right_up(L0,Pe,N,L,alpha)

implicit none
real :: L0,Pe,alpha
real :: N(4,6),L(4,6)
real :: s(5),c(5)
integer :: k
real :: rax,ray,rlx

L(2,1) = 2.0
do k = 1,5
L(1,k) = L0/5
L(2,k+1) = 2-(1.0/5*k)
L(3,k) = SQRT((L(1,k)**2)+(L(2,k)**2))
L(4,k) = 0.2/SIN(alpha)
end do

do k = 1,5
s(k) = L(1,6-k)/L(3,6-k)
c(k) = L(2,6-k)/L(3,6-k)
end do

rax = -L0*Pe/4.0
ray = Pe/2
rlx = -rax

do k = 1,5
if(k == 1) then
N(4,5) = 0
N(2,6) = Pe/2
N(3,5) = -N(2,6)/c(k)
N(1,5) = -N(3,5)*s(k)
else 
N(4,6-k) = ((N(4,7-k)*COS(alpha))+(N(3,7-k)*s(k-1)))/COS(alpha)
N(2,7-k) = (N(4,6-k)*SIN(alpha))-(N(3,7-k)*c(k-1))-(N(4,7-k)*SIN(alpha))
N(3,6-k) = -N(2,7-k)/c(k)
N(1,6-k) = N(1,7-k)-(N(3,6-k)*s(k))
end if
end do
N(2,1) = ray

end subroutine

subroutine right_down(L0,Pe,N,L,alpha)

implicit none
real :: L0,Pe,alpha
real :: N(4,6),L(4,6)
real :: s(5),c(5)
integer :: k
real :: rax,ray,rlx

do k = 1,5
L(1,k) = L0/5
L(2,k+1) = 2-(1.0/5*k)
L(3,k) = SQRT((L(1,k)**2)+(L(2,k+1)**2))
L(4,k) = 0.2/SIN(alpha)
end do
L(2,1) = 2.0

do k = 1,5
s(k) = L(1,6-k)/L(3,6-k)
c(k) = L(2,7-k)/L(3,6-k)
end do

rax = -L0*Pe/4.0
ray = pe/2
rlx = -rax

N(1,5) = 0
N(2,6) = 0
do k = 1,5
if(k == 1) then
N(3,6-k) = Pe*COS(alpha)/(2*((s(k)*SIN(alpha))+(COS(alpha)*c(k))))
N(2,6-k) = -c(k)*N(3,6-k)
N(4,6-k) = -(s(k)/COS(alpha))*N(3,6-k)
else if(k <  5) then
N(3,6-k) = -N(2,7-k)*COS(alpha)/((s(k)*SIN(alpha))+(COS(alpha)*c(k)))
N(2,6-k) = -c(k)*N(3,6-k)
N(4,6-k) = (N(2,7-k)+(N(4,7-k)*SIN(alpha))+(N(3,6-k)*c(k)))/SIN(alpha)
N(1,6-k) = N(3,7-k)*s(k-1)+N(1,7-k)
else 
N(1,6-k) = N(3,7-k)*s(k-1) + N(1,7-k)
N(3,6-k) = -(N(1,6-k)+rax)/s(5)
N(2,6-k) = ray - (N(3,6-k)*c(5))
N(4,6-k) = -rlx/COS(alpha)
end if
end do

end subroutine
