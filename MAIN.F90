program Julia2B
!
! This code is an example code for 2-d plane strain explicit time integration 
! method with iEFGM which is suggested by Saigal
! Also modified NBNM for weight factor is adopted
!
! 2000. 9. 3. To implement multi-domain/property treatment and perfectly plastic
! behavior. Also contact option will be inserted
! 
! 2000. 9. 16. Disable imsl library.
!
! 2000. 12. 27 TCC algorthm implemented
!
! 2000. 12. 30 hyperelastic behavior implemented
!
! 2001. 1.19 Modify internal force term to use 1st Piola-Kirhoffstress
!
!=============================================================================
implicit none
!=============================================================================
! 일단 5종류 설정.
real*8::MaterialProperty(3,5)
real*8::imsistress(3)
common/material/MaterialProperty
integer, parameter::N20 = 40
integer, dimension(:), allocatable::NumberofDOI, NumberOfGDOI, NodeExtF, &
     NodeInitVel, NodeOfMaster, NodeOfSlave, NodeFixed, NumberOfDOI2
integer, dimension(:,:), allocatable::NodeOfDOI, NodeOfCell, GaussOfDOI, &
     NodeOfDOI2
real*8, dimension(:), allocatable::MassOfNode, Jacobian
real*8, dimension(:,:), allocatable::InitCoordOfNode, CoordOfNode, &
     DisplOfNode, VelOfNode, PreVelOfNode, AccelOfNode, Phi4Node, dxPhi4Node,&
     dyPhi4Node, Phi4Gauss, dxPhi4Gauss, dyPhi4Gauss, ExtForce, &
     InitVelocity, CoordOfGauss, dxPhi, dyPhi
real*8::TimeIncrement, TimeIncrBack, TimeAccumulated, TimeToStop, TimeInterval
real::t0, t1, t2
character(len=15)::dummy, filename
character(len=10)::number
integer::NumberOfTotalNode, NumberOfCell, NumberOfExtForce, i, j, &
     AllocateStatus, NumberOfMaterial, NT, NumberOfVelCond, &
     NumberOfMasterNode, NumberOfSlaveNode, NumberOfFixedBC
t0 = 0.0
t1 = secnds(t0)
!==============================================================================
! open input/output file
!
open(unit=11,file="finemesh.inp")
open(unit=15,file="output.dat")
!
data number/'0123456789'/
filename = "output00.dis"
!=============================================================================
read(11,*) dummy
read(11,*) NumberOfTotalNode, NumberOfCell, NumberOfMaterial
NT = NumberOfTotalNode
allocate(CoordOfGauss(NumberOfCell*4,2), InitCoordOfNode(NT,2), &
     CoordOfNode(NT,2), DisplOfNode(NT,2), VelOfNode(NT,2),PreVelOfNode(NT,2),&
     AccelOfNode(NT,2), MassOfNode(NT), Phi4Node(N20,NT), dxPhi4Node(N20,NT),&
     dyPhi4Node(N20,NT), Phi4Gauss(NT,N20), dxPhi4Gauss(NT,N20), &
     dyPhi4Gauss(NT,N20), NumberOfDOI(NT), Jacobian(NumberOfCell), &
     NodeOfDOI(N20,NT), NumberOfGDOI(NT), GaussOfDOI(N20,NT), &
     NodeOfCell(NumberOfCell,4), dxPhi(N20,NumberOfCell*4), &
     dyPhi(N20,NumberOfCell*4), NodeOfDOI2(N20,NumberOfCell*4), &
     NumberOfDOI2(NumberOfCell*4), &
     stat = AllocateStatus)
if(AllocateStatus /= 0) stop "==== not enough memory ===="
!
! material property input
read(11,*) dummy
do i=1, NumberOfMaterial
   read(11,*) MaterialProperty(1,i), MaterialProperty(2,i), &
        MaterialProperty(3,i)
end do
!
! read final time
read(11,*) dummy
read(11,*) TimeToStop, TimeInterval
!
! Coordinate of Node input
read(11,*) dummy
do i=1,NumberOfTotalNode
   read(11,*) j, InitCoordOfNode(j,1), InitCoordOfNode(j,2)
end do
!
! Initially, Coordinates are same to initial value.
CoordOfNode = InitCoordOfNode
!
! Cell constitution
read(11,*) dummy
do i=1,NumberOfCell
   read(11,*) j, NodeOfCell(j,1), NodeOfCell(j,2), NodeOfCell(j,3), &
        NodeOfCell(j,4)
end do
!
! Fixed boundary condition
read(11,*) dummy
read(11,*) NumberOfFixedBC
allocate(NodeFixed(NumberOfFixedBC), stat=AllocateStatus)
if(AllocateStatus /= 0) stop "====== not enough memory - BC  ======="
do i=1, NumberOfFixedBC 
   read(11,*) NodeFixed(i)
end do
!
! External force per node
read(11,*) dummy
read(11,*) NumberOfExtForce
allocate(NodeExtF(NumberOfExtForce), ExtForce(NumberOfExtForce,2), &
     stat=AllocateStatus)
if(AllocateStatus /= 0) stop "====== not enough memory - 2nd ======="
do i=1, NumberOfExtForce
   read(11,*) NodeExtF(i), ExtForce(i,1), ExtForce(i,2)
end do
!
! Initial velocity condition
read(11,*) dummy
read(11,*) NumberOfVelCond
allocate(NodeInitVel(NumberOfVelCond), InitVelocity(NumberOfVelCond,2), &
     stat=AllocateStatus)
if(AllocateStatus /= 0) stop "====== not enough memory - 3rd ======="
do i=1, NumberOfVelCond
   read(11,*) NodeInitVel(i), InitVelocity(i,1), InitVelocity(i,2)
end do
!
! Contact data input
! master node set
read(11,*) dummy
read(11,*) NumberOfMasterNode
allocate(NodeOfMaster(NumberOfMasterNode), stat=AllocateStatus)
if(AllocateStatus /= 0) stop "====== not enough memory - 4th ======="
do i=1, NumberOfMasterNode
   read(11,*) NodeOfMaster(i)
end do
!
! slave node set
read(11,*) dummy
read(11,*) NumberOfSlaveNode
allocate(NodeOfSlave(NumberOfSlaveNode), stat=AllocateStatus)
if(AllocateStatus /= 0) stop "====== not enough memory - 5th ======="
do i=1, NumberOfSlaveNode
   read(11,*) NodeOfSlave(i)
end do

!==============================================================================
!
! Integration point estimation
call GaussPoint(NumberOfTotalNode, InitCoordOfNode, NumberOfCell, &
     NodeOfCell, CoordOfGauss, Jacobian)
!
! Shape function 계산	
!
! For node	 
call ShapeFtnNode(N20, NumberOfTotalNode, InitCoordOfNode, &
     Phi4Node, dxPhi4Node, dyPhi4Node, NumberOfDOI, NodeOfDOI)
!
! For Gauss point
call ShapeFtnGauss(N20, NT, NumberOfCell, InitCoordOfNode, CoordOfGauss, &
     Phi4Gauss, dxPhi4Gauss, dyPhi4Gauss, NumberOfGDOI, GaussOfDOI, Jacobian, &
     dxPhi, dyPhi, NodeOfDOI2, NumberOfDOI2)
print *, "shape routine end"
!
!
call MassDistribution(N20, NumberOfTotalNode, NumberOfCell, MassOfNode, &
     NumberOfDOI, NodeOfDOI, NumberOfGDOI, GaussOfDOI, Jacobian, Phi4Gauss)
DisplOfNode = 0.0
VelOfNode = 0.0
PreVelOfNode = 0.0
AccelOfNode = 0.0
TimeIncrement = 0.0
TimeAccumulated = 0.0
!
! initial condition 
do i=1, NumberOfVelCond
	PreVelOfNode(NodeInitVel(i),1) = InitVelocity(i,1)
	PreVelOfNode(NodeInitVel(i),2) = InitVelocity(i,2)
end do
!
! Setting for interval output
j =1 
!
!100 call IncrementOfTime(N20, NumberOfTotalNode, CoordOfNode, TimeIncrement, &
!         TimeIncrBack, TimeAccumulated, NumberOfDOI, NodeOfDOI)
100 TimeIncrement = 1.e-6
TimeIncrBack = TimeIncrement
TimeAccumulated = TimeAccumulated + 1.e-6
!
call Predict(NumberOfTotalNode, DisplOfNode,VelOfNode, PreVelOfNode, &
     AccelOfNode, TimeIncrement,InitCoordOfNode, CoordOfNode, &
     NodeInitVel,InitVelocity,NumberOfVelCond)
!
! position update
do i=1, NumberOfTotalNode
	CoordOfNode(i,1) = InitCoordOfNode(i,1) + DisplOfNode(i,1)
	CoordOfNode(i,2) = InitCoordOfNode(i,2) + DisplOfNode(i,2)
end do
!
! 
call Acceleration(N20, NumberOfTotalNode, NumberOfCell, DisplOfNode, &
     AccelOfNode, &
     MassOfNode, dxPhi4Node, dyPhi4Node, dxPhi4Gauss, dyPhi4Gauss, &
     NumberofDOI, NodeOfDOI, NumberOfGDOI, GaussOfDOI, Jacobian, &
     NumberOfExtForce, NodeExtF, ExtForce, CoordOfNode, &
     NumberOfMasterNode, NumberOfSlaveNode, NodeOfMaster, NodeOfSlave, &
     TimeIncrement, dxPhi, dyPhi, NodeOfDOI2, NumberOfDOI2 )
!
! Fixed B/C
do i=1, NumberOfFixedBC 
   AccelOfNode(NodeFixed(i),1) = 0.0
   AccelOfNode(NodeFixed(i),2) = 0.0
end do
!
! Constant B/C
do i=1, NumberOfVelCond
   AccelOfNode(NodeInitVel(i),1) = 0.0
   AccelOfNode(NodeInitVel(i),2) = 0.0
   PreVelOfNode(NodeInitVel(i),1) = InitVelocity(i,1)
   PreVelOfNode(NodeInitVel(i),2) = InitVelocity(i,2)
end do
!
call Correct(NumberOfTotalNode, VelOfNode, AccelOfNode, TimeIncrement)
!
! write interval result
if(TimeAccumulated > TimeInterval*j) then
   if( int(j/10) < 1) then
      filename(7:7) = number(1:1)
      filename(8:8) = number(j+1:j+1)
   else
      filename(7:7) = number(int(j/10)+1:int(j/10)+1)
      filename(8:8) = number(j - int(j/10)*10 +1:j - int(j/10)*10 +1)
   end if
   open(unit=16, file=filename)

   call Postprocess(16, NumberOfTotalNode, DisplOfNode, TimeAccumulated, &
        VelOfNode)
close(16)
j = j + 1
end if
!
if(TimeAccumulated < TimeToStop) then
	goto 100
end if
call Postprocess(15, NumberOfTotalNode, CoordOfNode, TimeAccumulated, &
     VelOfNode)
!
t2 = secnds(t1)
print *, 'used cpu time is', t2 , 'seconds'
close(11)
close(15)
stop 
end program Julia2B
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine Postprocess(N15, NumberOfTotalNode, DisplOfNode, TimeAccumulated, &
     VelOfNode)
!==============================================================================
implicit none
!==============================================================================
integer::NumberOfTotalNode, i, N15
real*8::DisplOfNode(NumberOfTotalNode,2), TimeAccumulated, &
     VelOfNode(NumberOfTotalNode,2)
INTEGER::NDMAX, NWidth
REAL*8::DefMax, xxx, yyy
CHARACTER(LEN=80)::title, subtitle, subtitle2
print *, "output is processed at time =", TimeAccumulated
!do i=1, NumberOfTotalNode
!	write(N15,*) i , CoordOfNode(i,1), CoordOfNode(i,2), VelOfNode(i,2)
!end do
NWidth = 2     ! DOF of data
NDMAX = NumberOfTotalNode     ! Maximum number of nodes
DefMax = 1.    ! Maximum data value
title = "title"
subtitle = "subtitle"
subtitle2 = "subtitle2"
write(16,2001) title
write(16,2002) NumberOfTotalNode, NumberOfTotalNode, DefMax, NDMAX, NWidth
write(16,2001) subtitle
write(16,2001) subtitle2
do i=1, NumberOfTotalNode
   if (abs(DisplOfNode(i,1)) < 1.e-9 ) then
      xxx = 0.0
   else  
      xxx = DisplOfNode(i,1)
   end if
   if (abs(DisplOfNode(i,2)) < 1.e-9 ) then
      yyy = 0.0
   else 
      yyy = DisplOfNode(i,2)
   end if
   write(16,2003) i, xxx, yyy
end do
2001 format(80A1)
2002 format(2I9,E15.6E1,2I9)
2003 format(I8,  (5E13.7E1))
return
end subroutine Postprocess
!

