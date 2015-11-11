!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
subroutine MassDistribution(N20, NT, NumberOfCell, MassOfNode, &
     NumberOfDOI, NodeOfDOI, NumberOfGDOI, GaussOfDOI, Jacobian, Phi4Gauss)
!==============================================================================
implicit none
!==============================================================================
real*8::MaterialProperty(3,5)
common/material/MaterialProperty
integer::N20, NT, NumberOfTotalNode, NumberOfCell, i, j, k, NumberofDOI(NT), &
     NodeOfDOI(N20, NT), NumberOfGDOI(NT), GaussOfDOI(N20,NT)
real*8::MassOfNode(NT), Jacobian(NumberOfCell), Phi4Gauss(NT,N20),&
     density
!==============================================================================
NumberOfTotalNode = NT
MassOfNode = 0.0
density = 1.4089e-4/32.17
do i=1, NumberOfTotalNode
   do j=1, NumberOfGDOI(i)
      !
      ! k is cell of jth gauss point
      k = int((GaussOfDOI(j,i) - 1)/4) + 1
      MassOfNode(i) = MassOfNode(i) + Phi4Gauss(i,j)*Jacobian(k)*density
   end do
end do
!
return
end subroutine MassDistribution
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
subroutine IncrementOfTime(NumberOfTotalNode, CoordOfNode, TimeIncrement, &
     TimeIncrBack, TimeAccumulated, NumberOfDOI, NodeOfDOI)
!==============================================================================
implicit none
!==============================================================================
real*8::MaterialProperty(3,5)
common/material/MaterialProperty
integer::NumberOfTotalNode, i,j, &
     NumberofDOI(NumberOfTotalNode), NodeOfDOI(20,NumberOfTotalNode)
real*8::CoordOfNode(NumberOfTotalNode,2), length, TimeIncrement, &
     TimeIncrBack, TimeAccumulated, length2(20), short
!==============================================================================
TimeIncrBack = TimeIncrement
length = 1.0
do i=1,NumberOfTotalNode
   length2 = 1.0
   do j=1,NumberOfDOI(i)
      if(i /= NodeOfDOI(j,i)) then
         length2(j) = ((CoordOfNode(i,1) - &
              CoordOfNode(NodeOfDOI(j,i),1))**2 + &
              (CoordOfNode(i,2) - CoordOfNode(NodeOfDOI(j,i),2))**2 )
      end if
   end do
   short = minval(length2)
   !아래 부분 고칠것. 매 루프마다 계산해야할 양이 많음.
   length = min(short*MaterialProperty(1,1)/ &
        MaterialProperty(2,1), length)
end do
TimeIncrement = sqrt(length)
TimeAccumulated = TimeAccumulated + TimeIncrement
return
end subroutine
!
! Gauss-point estimation loop
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine GaussPoint(NT, InitCoordOfNode, NumberOfCell, NodeOfCell, &
     CoordOfGauss, Jacobian)
!==============================================================================
implicit none
!==============================================================================
!
INTEGER::i,j, NT, NumberOfTotalNode, NumberOfCell, NodeOfCell(NumberOfCell,4)
REAL*8::InitCoordOfNode(NT,2), CoordOfGauss(NumberOfCell*4,2), psi1, psi2, &
     psi3, psi4, psi(4,4), t, xi(4), eta(4), x1, x2, x3, x4, y1, y2, y3, y4, &
     Jacobian(NumberOfCell), For4Pjacobian
!
! Gauss point position in local coordinate
t = 1./3.
xi(1) =  - sqrt(t)
xi(2) =  sqrt(t)
xi(3) = xi(2)
xi(4) = xi(1)
eta(1) = xi(1)
eta(2) = xi(1)
eta(3) = xi(2)
eta(4) = xi(2)
!
! factor for gauss-point estimation
do i=1, 4
   psi(1,i) = psi1(xi(i),eta(i))
   psi(2,i) = psi2(xi(i),eta(i))
   psi(3,i) = psi3(xi(i),eta(i))
   psi(4,i) = psi4(xi(i),eta(i))
end do
!
! Gauss-point estimation loop
do i=1, NumberOfCell
   x1 = InitCoordOfNode(NodeOfCell(i,1),1)
   y1 = InitCoordOfNode(NodeOfCell(i,1),2)
   x2 = InitCoordOfNode(NodeOfCell(i,2),1)
   y2 = InitCoordOfNode(NodeOfCell(i,2),2)
   x3 = InitCoordOfNode(NodeOfCell(i,3),1)
   y3 = InitCoordOfNode(NodeOfCell(i,3),2)
   x4 = InitCoordOfNode(NodeOfCell(i,4),1)
   y4 = InitCoordOfNode(NodeOfCell(i,4),2)
   CoordOfGauss((i-1)*4+1,1) = x1*psi(1,1)+x2*psi(2,1)+x3*psi(3,1)+x4*psi(4,1)
   CoordOfGauss((i-1)*4+2,1) = x1*psi(1,2)+x2*psi(2,2)+x3*psi(3,2)+x4*psi(4,2)
   CoordOfGauss((i-1)*4+3,1) = x1*psi(1,3)+x2*psi(2,3)+x3*psi(3,3)+x4*psi(4,3)
   CoordOfGauss((i-1)*4+4,1) = x1*psi(1,4)+x2*psi(2,4)+x3*psi(3,4)+x4*psi(4,4)
   CoordOfGauss((i-1)*4+1,2) = y1*psi(1,1)+y2*psi(2,1)+y3*psi(3,1)+y4*psi(4,1)
   CoordOfGauss((i-1)*4+2,2) = y1*psi(1,2)+y2*psi(2,2)+y3*psi(3,2)+y4*psi(4,2)
   CoordOfGauss((i-1)*4+3,2) = y1*psi(1,3)+y2*psi(2,3)+y3*psi(3,3)+y4*psi(4,3)
   CoordOfGauss((i-1)*4+4,2) = y1*psi(1,4)+y2*psi(2,4)+y3*psi(3,4)+y4*psi(4,4)
   Jacobian(i) = For4Pjacobian(x1,x2,x3,x4,y1,y2,y3,y4)
end do
return
end subroutine GaussPoint
!
! shape function for convenctional FEM in 4-point integration
function psi1(x,y)
  real*8::x, y, psi1
  psi1 = 0.25*(1. - x) * (1 - y)
end function psi1
function psi2(x,y)
  real*8::x, y, psi2
  psi2 = 0.25*(1. + x) * (1 - y)
end function psi2
function psi3(x,y)
  real*8::x, y, psi3
  psi3 = 0.25*(1. + x) * (1 + y)
end function psi3
function psi4(x,y)
  real*8::x, y, psi4
  psi4 = 0.25*(1. - x) * (1 + y)
end function psi4
!
! Jacobian estimation for each bucket
function For4Pjacobian(x1,x2,x3,x4,y1,y2,y3,y4)
  real*8::x1, x2, x3, x4, y1, y2, y3, y4, For4Pjacobian
  For4Pjacobian = ((-x1 + x2 + x3 - x4)*(-y1 -y2 + y3 + y4) - &
       (-x1 - x2 + x3 + x4)*(-y1 + y2 + y3 - y4))/16.0
end function For4Pjacobian


