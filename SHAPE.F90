subroutine ShapeFtnNode(N20, NT, InitCoordOfNode, &
     Phi4Node, dxPhi4Node, dyPhi4Node, NumberOfDOI, NodeOfDOI)
!==============================================================================
implicit none
!==============================================================================
!
integer::i, j, NT, N20, NumberOfTotalNode, NumberOfDOI(NT), NodeOfDOI(N20,NT),&
     k, l
real*8::InitCoordOfNode(NT,2), Phi4Node(N20,NT), dxPhi4Node(N20,NT), &
     dyPhi4Node(N20,NT), rd, w(N20), dw(N20,2), dwdr(N20), SumOfWeight, &
     dSumOfWeightdx, dSumOfWeightdy, distance, summation, Determinant,&
     DilationParameterA, DifferenceOfX, DifferenceOfY
!==============================================================================
NumberOfTotalNode = NT
NumberOfDOI = 0
NodeOfDOI = 0
DilationParameterA = 1.5
rd = DilationParameterA*abs(InitCoordOfNode(2,1)-InitCoordOfNode(1,1))
Phi4Node = 0.0
dxPhi4Node = 0.0
dyPhi4Node = 0.0
!
FirstLoop: do i=1, NumberOfTotalNode
!
! i mean 임의의 x변수
SumOfWeight = 0.0
dSumOfWeightdx = 0.0
dSumOfWeightdy = 0.0
w = 0.0
dwdr = 0.0
dw = 0.0
k = 0
do j=1,NumberOfTotalNode
!
! Material/Domain check
!   if(NFlagOfMatrl(i) == NFlagOfMatrl(j) ) then
!
! loop starts
      DifferenceOfX = InitCoordOfNode(i,1) - InitCoordOfNode(j,1)
      DifferenceOfY = InitCoordOfNode(i,2) - InitCoordOfNode(j,2)
      distance = sqrt(DifferenceOfX**2 + DifferenceOfY**2)
      if (i == j) then
         distance = 1.0e-6
         k = k + 1
         NodeOfDOI(k,i) = j
         w(k) = (rd/distance)**2 - (2.0*rd/distance) + 1.0
         dwdr(k) = ((-2.0*(rd**2)/(distance**4)) + &
              (2.0*rd/(distance**3))) 
         dw(k,1) = dwdr(k)*abs(DifferenceOfX)/distance
         dw(k,2) = dwdr(k)*abs(DifferenceOfY)/distance
         SumOfWeight = SumOfWeight + w(k)
         dSumOfWeightdx = dSumOfWeightdx + dw(k,1)
         dSumOfWeightdy = dSumOfWeightdy + dw(k,2)
      else if (distance <= rd) then
         k = k + 1
         NodeOfDOI(k,i) = j
         w(k) = (rd/distance)**2 - (2.0*rd/distance) + 1.0
         dwdr(k) = ((-2.0*(rd**2)/(distance**4)) + &
              (2.0*rd/(distance**3))) 
         dw(k,1) = dwdr(k)*abs(DifferenceOfX)/distance
         dw(k,2) = dwdr(k)*abs(DifferenceOfY)/distance
         SumOfWeight = SumOfWeight + w(NumberOfDOI(i))
         dSumOfWeightdx = dSumOfWeightdx + dw(k,1)
         dSumOfWeightdy = dSumOfWeightdy + dw(k,2)
      else
!         w(j) = 0.0
!         dwdr(j) = 0.0
      end if
!   else 
!      w(j) = 0.0
!      dwdr(j) = 0.0
!   end if
end do
NumberOfDOI(i) = k
call ShapeFtnSub(N20, NumberOfDOI(i), NodeOfDOI, NT, SumOfWeight, &
     dSumOfWeightdx, dSumOfWeightdy, w, dw, InitCoordOfNode, i, &
     Phi4Node, dxPhi4Node, dyPhi4Node)
end do Firstloop
!loop end for arbitrary nodal point x
!
return
end subroutine ShapeFtnNode
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Shape function routine for Gauss-point
subroutine ShapeFtnGauss(N20, NT, NumberOfCell, InitCoordOfNode, &
     CoordOfGauss, Phi4Gauss, dxPhi4Gauss, dyPhi4Gauss, NumberOfGDOI, &
     GaussOfDOI, Jacobian, dxPhi, dyPhi, NodeOfDOI2, NumberOfDOI2)
!==============================================================================
implicit none
!==============================================================================
!
integer::i, j, NT, N20, NumberOfTotalNode, NumberOfCell, Node, &
     NumberOfDOI2(NumberOfCell*4), NodeOfDOI2(N20,NumberOfCell*4),&
     k, l, NumberOfGDOI(NT), GaussOfDOI(N20,NT) 
real*8::InitCoordOfNode(NT,2), Phi4Gauss(NT,N20), dxPhi4Gauss(NT,N20), &
     dyPhi4Gauss(NT,N20), rd, w(N20), dw(N20,2), dwdr(N20), SumOfWeight, &
     dSumOfWeightdx, dSumOfWeightdy, distance, summation, Determinant,&
     DilationParameterA, DifferenceOfX, DifferenceOfY, Jacobian(NumberOfCell),&
     Phi(N20,NumberOfCell*4), dxPhi(N20,NumberOfCell*4), &
     dyPhi(N20,NumberOfCell*4), CoordOfGauss(NumberOfCell*4,2)
!==============================================================================
NumberOfTotalNode = NT
NumberOfDOI2 = 0
NodeOfDOI2 = 0
NumberOfGDOI = 0
GaussOfDOI = 0
DilationParameterA = 1.0
rd = DilationParameterA*abs(InitCoordOfNode(2,1)-InitCoordOfNode(1,1))
Phi4Gauss = 0.0
dxPhi4Gauss = 0.0
dyPhi4Gauss = 0.0
!
FirstLoop: do i=1, NumberOfCell*4
!
! i mean 임의의 x변수
SumOfWeight = 0.0
dSumOfWeightdx = 0.0
dSumOfWeightdy = 0.0
w = 0.0
dwdr = 0.0
dw = 0.0
k = 0
do j=1,NumberOfTotalNode
!
! Material/Domain check
!   if(NFlagOfMatrl(i) == NFlagOfMatrl(j) ) then
!
! loop starts
      DifferenceOfX = CoordOfGauss(i,1) - InitCoordOfNode(j,1)
      DifferenceOfY = CoordOfGauss(i,2) - InitCoordOfNode(j,2)
      distance = sqrt(DifferenceOfX**2 + DifferenceOfY**2)
      if (distance <= rd) then
         k = k + 1
         NodeOfDOI2(k,i) = j
         w(k) = (rd/distance)**2 - (2.0*rd/distance) + 1.0
         dwdr(k) = ((-2.0*(rd**2)/(distance**4)) + &
              (2.0*rd/(distance**3))) 
         dw(k,1) = dwdr(k)*abs(DifferenceOfX)/distance
         dw(k,2) = dwdr(k)*abs(DifferenceOfY)/distance
         SumOfWeight = SumOfWeight + w(k)
         dSumOfWeightdx = dSumOfWeightdx + dw(k,1)
         dSumOfWeightdy = dSumOfWeightdy + dw(k,2)
      else
!         w(j) = 0.0
!         dwdr(j) = 0.0
      end if
!   else 
!      w(j) = 0.0
!      dwdr(j) = 0.0
!   end if
end do
NumberOfDOI2(i) = k
call ShapeFtnSub2(N20, NumberOfDOI2(i), NodeOfDOI2, NT, NumberOfCell, &
     SumOfWeight, dSumOfWeightdx, dSumOfWeightdy, w, dw, InitCoordOfNode, i, &
     Phi, dxPhi, dyPhi, CoordOfGauss(i,1), CoordOfGauss(i,2))
end do Firstloop
!
! Redistribution
do i = 1, NumberOfCell*4
   do j = 1, NumberOfDOI2(i)
      Node = NodeOfDOI2(j,i)
      NumberOfGDOI(Node) = NumberOfGDOI(Node) + 1
      GaussOfDOI(NumberOfGDOI(Node), Node) = i
      Phi4Gauss(Node, NumberOfGDOI(Node)) = Phi(j,i)
      dxPhi4Gauss(Node, NumberOfGDOI(Node)) = dxPhi(j,i)
      dyPhi4Gauss(Node, NumberOfGDOI(Node)) = dyPhi(j,i)
   end do
end do
!
return
end subroutine ShapeFtnGauss
