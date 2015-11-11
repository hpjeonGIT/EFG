subroutine ShapeFtnSub(N20, NN, NodeOfDOI, NT, SumOfWeight, &
     dSumOfWeightdx, dSumOfWeightdy, w, dw, InitCoordOfNode, NodeCurrent, &
     Phi4Node, dxPhi4Node, dyPhi4Node)     
!==============================================================================
implicit none
!==============================================================================
integer::i, j, k, N20, NumberOfDOI, NT, NN, Node, NodeOfDOI(N20,NT), &
     NodeCurrent
real*8::InitCoordOfNode(NT,2), Phi4Node(N20,NT), dxPhi4Node(N20,NT), &
     dyPhi4Node(N20,NT), DilationParameterA, rd, w(N20), dw(N20,2), &
     SumOfWeight, dSumOfWeightdx, dSumOfWeightdy, v(NN), dvdx(NN), dvdy(NN), &
     u1, u2rear, u2(NN), du2dx, du2dy, u3rear, u3(NN), du3dx, du3dy, &
     matrixI(NN,NN), matrixV(NN,NN), matrixdVdx(NN,NN), matrixdVdy(NN,NN), &
     B(2,NN), dBdx(2,NN), dBdy(2,NN), temp1(NN,1), temp2(NN,1), temp3(NN,1), &
     g(2,1), dgdx(2,1), dgdy(2,1), A(2,2), Ainv(2,2), distance, dAdx(2,2), &
     dAdy(2,2), summation, Determinant, DifferenceOfX, DifferenceOfY
!==============================================================================
!sum = 0.0
NumberOfDOI = NN
u1 = 1.0 / sqrt(SumOfWeight)
!
! Initialization ==========================================
u2rear = 0.0
u3rear = 0.0
du2dx = 0.0
du2dy = 0.0
du3dx = 0.0
du3dy = 0.0
v = 0.0
matrixV = 0.0
matrixdVdx = 0.0
matrixdVdy = 0.0
matrixI = 0.0
do k=1,NN
   matrixI(k,k) = 1.0
end do
! End of Initialization ========================================
!
do j=1, NumberOfDOI
   Node = NodeOfDOI(j, NodeCurrent) 
   v(j) = w(j) / SumOfWeight
   dvdx(j) = (dw(j,1)*SumOfWeight - w(j)*dSumOfWeightdx)/(SumOfWeight**2)
   dvdy(j) = (dw(j,2)*SumOfWeight - w(j)*dSumOfWeightdy)/(SumOfWeight**2)
   u2rear = u2rear + (InitCoordOfNode(Node,1)*v(j))
   u3rear = u3rear + (InitCoordOfNode(Node,2)*v(j))
   du2dx = du2dx - (InitCoordOfNode(Node,1)*dvdx(j))
   du2dy = du2dy - (InitCoordOfNode(Node,1)*dvdy(j))
   du3dx = du3dx - (InitCoordOfNode(Node,2)*dvdx(j))
   du3dy = du3dy - (InitCoordOfNode(Node,2)*dvdy(j))
   do k=1, NumberOfDOI
      matrixV(k,j) = v(j)
      matrixdVdx(k,j) = dvdx(j)
      matrixdVdy(k,j) = dvdy(j)
   end do
end do
!
B = 0.0
dBdx = 0.0
dBdy = 0.0
u2 = 0.0
u3 = 0.0
do k=1, NumberOfDOI
   Node = NodeOfDOI(k, NodeCurrent)
   u2(k) = InitCoordOfNode(Node,1) - u2rear
   u3(k) = InitCoordOfNode(Node,2) - u3rear
   B(1,k) = w(k)*u2(k)
   B(2,k) = w(k)*u3(k)
   dBdx(1,k) = dw(k,1)*u2(k) + w(k)*du2dx
   dBdx(2,k) = dw(k,1)*u3(k) + w(k)*du3dx
   dBdy(1,k) = dw(k,2)*u2(k) + w(k)*du2dy
   dBdy(2,k) = dw(k,2)*u3(k) + w(k)*du3dy
end do
!
g(1,1) = InitCoordOfNode(NodeCurrent,1) - u2rear
g(2,1) = InitCoordOfNode(NodeCurrent,2) - u3rear
dgdx(1,1) = 1.0 + du2dx
dgdx(2,1) = du3dx
dgdy(1,1) = du2dy
dgdy(2,1) = 1.0 + du3dy
!
A = 0.0
dAdx = 0.0
dAdy = 0.0
do k=1,NumberOfDOI
   Node = NodeOfDOI(k,NodeCurrent)
   if (Node == NodeCurrent) then
      A = A + 0.0
      dAdx = dAdx + 0.0
      dAdy = dAdy + 0.0 
   else
      A(1,1) = A(1,1) + (u2(k)*w(k)*u2(k))
      A(1,2) = A(1,2) + (u2(k)*w(k)*u3(k))
      A(2,2) = A(2,2) + (u3(k)*w(k)*u3(k))
      dAdx(1,1) = dAdx(1,1) + (2.0*du2dx*w(k)*u2(k)) + (u2(k)*dw(k,1)*u2(k))
      dAdx(1,2) = dAdx(1,2) + (du2dx*w(k)*u3(k)) + (u2(k)*dw(k,1)*u3(k)) + &
           (u2(k)*w(k)*du3dx)
      dAdx(2,2) = dAdx(2,2) + (2.0*du3dx*w(k)*u3(k)) + (u3(k)*dw(k,1)*u3(k))
      dAdy(1,1) = dAdy(1,1) + (2.0*du2dy*w(k)*u2(k)) + (u2(k)*dw(k,2)*u2(k))
      dAdy(1,2) = dAdy(1,2) + (du2dy*w(k)*u3(k)) + (u2(k)*dw(k,2)*u3(k)) + &
           (u2(k)*w(k)*du3dy)
      dAdy(2,2) = dAdy(2,2) + (2.0*du3dy*w(k)*u3(k)) + (u3(k)*dw(k,2)*u3(k))
     end if
  end do
  A(2,1) = A(1,2)
  dAdx(2,1) = dAdx(1,2)
  dAdy(2,1) = dAdy(1,2)
  !call DLINRG(2, A, 2, Ainv, 2)
  Determinant = (A(1,1)*A(2,2)) - (A(1,2)*A(2,1))
  Ainv(1,1) = A(2,2)/Determinant
  Ainv(2,2) = A(1,1)/Determinant
  Ainv(1,2) = -1.0*A(1,2)/Determinant
  Ainv(2,1) = -1.0*A(2,1)/Determinant
  !
  !temp=>temp(m,1)
  temp1 = matmul(matmul(matmul( transpose(matrixI - matrixV), transpose(B)), &
       Ainv),g)
  temp2 = matmul(transpose(matrixI - matrixV), (matmul(matmul(transpose(dBdx),&
       Ainv), g) - &
       matmul(matmul(transpose(B),Ainv), matmul(matmul(dAdx, Ainv),g)) + &
       matmul(matmul(transpose(B),Ainv),dgdx) ) ) - &
       matmul(matmul(matmul( transpose(matrixdVdx), transpose(B)), Ainv),g)
  temp3 = matmul(transpose(matrixI - matrixV), (matmul(matmul(transpose(dBdy),&
       Ainv), g) - &
       matmul(matmul(transpose(B),Ainv), matmul(matmul(dAdy, Ainv),g)) + &
       matmul(matmul(transpose(B),Ainv),dgdy) ) ) - &
       matmul(matmul(matmul( transpose(matrixdVdy), transpose(B)), Ainv),g)
do j=1,NumberOfDOI
!   Node = NodeOfDOI(j,NodeCurrent)
!
!  여기서 ( ,x)는 임의의 x, (j, )는 루프상의 지점을 의미한다.
!
   Phi4Node(j,NodeCurrent) = v(j) + temp1(j,1)
   dxPhi4Node(j,NodeCurrent) = dvdx(j) + temp2(j,1)
   dyPhi4Node(j,NodeCurrent) = dvdy(j) + temp3(j,1)
!   sum = sum + Phi4Node(j,NodeCurrent)
end do
!print *, NodeCurrent,"=", sum
!sum = 0.0
return
end subroutine ShapeFtnSub
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
! sub-shape function routine for Gauss-point
subroutine ShapeFtnSub2(N20, NN, NodeOfDOI2, NT, NumberOfCell, SumOfWeight, &
     dSumOfWeightdx, dSumOfWeightdy, w, dw, InitCoordOfNode, GaussCurrent, &
     Phi, dxPhi, dyPhi, Coord1, Coord2)
!==============================================================================
implicit none
!==============================================================================
integer::i, j, k, N20, NT, NN, NumberOfCell, Node, NumberOfDOI2, &
     NodeOfDOI2(N20,NumberOfCell*4), GaussCurrent
real*8::InitCoordOfNode(NT,2), DilationParameterA, rd, w(N20), dw(N20,2), &
     SumOfWeight, dSumOfWeightdx, dSumOfWeightdy, v(NN), dvdx(NN), dvdy(NN), &
     u1, u2rear, u2(NN), du2dx, du2dy, u3rear, u3(NN), du3dx, du3dy, &
     matrixI(NN,NN), matrixV(NN,NN), matrixdVdx(NN,NN), matrixdVdy(NN,NN), &
     B(2,NN), dBdx(2,NN), dBdy(2,NN), temp1(NN,1), temp2(NN,1), temp3(NN,1), &
     g(2,1), dgdx(2,1), dgdy(2,1), A(2,2), Ainv(2,2), distance, dAdx(2,2), &
     dAdy(2,2), summation, Determinant, DifferenceOfX, DifferenceOfY, &
     Phi(N20,NumberOfCell*4), dxPhi(N20,NumberOfCell*4), &
     dyPhi(N20,NumberOfCell*4), Coord1, Coord2
!==============================================================================
!sum = 0.0
NumberOfDOI2 = NN
u1 = 1.0 / sqrt(SumOfWeight)
!
! Initialization ==========================================
u2rear = 0.0
u3rear = 0.0
du2dx = 0.0
du2dy = 0.0
du3dx = 0.0
du3dy = 0.0
v = 0.0
matrixV = 0.0
matrixdVdx = 0.0
matrixdVdy = 0.0
matrixI = 0.0
do k=1,NN
   matrixI(k,k) = 1.0
end do
! End of Initialization ========================================
!
do j=1, NumberOfDOI2
   Node = NodeOfDOI2(j, GaussCurrent) 
   v(j) = w(j) / SumOfWeight
   dvdx(j) = (dw(j,1)*SumOfWeight - w(j)*dSumOfWeightdx)/(SumOfWeight**2)
   dvdy(j) = (dw(j,2)*SumOfWeight - w(j)*dSumOfWeightdy)/(SumOfWeight**2)
   u2rear = u2rear + (InitCoordOfNode(Node,1)*v(j))
   u3rear = u3rear + (InitCoordOfNode(Node,2)*v(j))
   du2dx = du2dx - (InitCoordOfNode(Node,1)*dvdx(j))
   du2dy = du2dy - (InitCoordOfNode(Node,1)*dvdy(j))
   du3dx = du3dx - (InitCoordOfNode(Node,2)*dvdx(j))
   du3dy = du3dy - (InitCoordOfNode(Node,2)*dvdy(j))
   do k=1, NumberOfDOI2
      matrixV(k,j) = v(j)
      matrixdVdx(k,j) = dvdx(j)
      matrixdVdy(k,j) = dvdy(j)
   end do
end do
!
B = 0.0
dBdx = 0.0
dBdy = 0.0
u2 = 0.0
u3 = 0.0
do k=1, NumberOfDOI2
   Node = NodeOfDOI2(k, GaussCurrent)
   u2(k) = InitCoordOfNode(Node,1) - u2rear
   u3(k) = InitCoordOfNode(Node,2) - u3rear
   B(1,k) = w(k)*u2(k)
   B(2,k) = w(k)*u3(k)
   dBdx(1,k) = dw(k,1)*u2(k) + w(k)*du2dx
   dBdx(2,k) = dw(k,1)*u3(k) + w(k)*du3dx
   dBdy(1,k) = dw(k,2)*u2(k) + w(k)*du2dy
   dBdy(2,k) = dw(k,2)*u3(k) + w(k)*du3dy
end do
!
g(1,1) = Coord1 - u2rear
g(2,1) = Coord2 - u3rear
dgdx(1,1) = 1.0 + du2dx
dgdx(2,1) = du3dx
dgdy(1,1) = du2dy
dgdy(2,1) = 1.0 + du3dy
!
A = 0.0
dAdx = 0.0
dAdy = 0.0
do k=1,NumberOfDOI2
   Node = NodeOfDOI2(k,GaussCurrent)
   A(1,1) = A(1,1) + (u2(k)*w(k)*u2(k))
   A(1,2) = A(1,2) + (u2(k)*w(k)*u3(k))
   A(2,2) = A(2,2) + (u3(k)*w(k)*u3(k))
   dAdx(1,1) = dAdx(1,1) + (2.0*du2dx*w(k)*u2(k)) + (u2(k)*dw(k,1)*u2(k))
   dAdx(1,2) = dAdx(1,2) + (du2dx*w(k)*u3(k)) + (u2(k)*dw(k,1)*u3(k)) + &
        (u2(k)*w(k)*du3dx)
   dAdx(2,2) = dAdx(2,2) + (2.0*du3dx*w(k)*u3(k)) + (u3(k)*dw(k,1)*u3(k))
   dAdy(1,1) = dAdy(1,1) + (2.0*du2dy*w(k)*u2(k)) + (u2(k)*dw(k,2)*u2(k))
   dAdy(1,2) = dAdy(1,2) + (du2dy*w(k)*u3(k)) + (u2(k)*dw(k,2)*u3(k)) + &
        (u2(k)*w(k)*du3dy)
   dAdy(2,2) = dAdy(2,2) + (2.0*du3dy*w(k)*u3(k)) + (u3(k)*dw(k,2)*u3(k))
end do
A(2,1) = A(1,2)
dAdx(2,1) = dAdx(1,2)
dAdy(2,1) = dAdy(1,2)
Determinant = (A(1,1)*A(2,2)) - (A(1,2)*A(2,1))
Ainv(1,1) = A(2,2)/Determinant
Ainv(2,2) = A(1,1)/Determinant
Ainv(1,2) = -1.0*A(1,2)/Determinant
Ainv(2,1) = -1.0*A(2,1)/Determinant
!
!temp=>temp(m,1)
temp1 = matmul(matmul(matmul( transpose(matrixI - matrixV), transpose(B)), &
     Ainv),g)
temp2 = matmul(transpose(matrixI - matrixV), (matmul(matmul(transpose(dBdx),&
     Ainv), g) - &
     matmul(matmul(transpose(B),Ainv), matmul(matmul(dAdx, Ainv),g)) + &
     matmul(matmul(transpose(B),Ainv),dgdx) ) ) - &
     matmul(matmul(matmul( transpose(matrixdVdx), transpose(B)), Ainv),g)
temp3 = matmul(transpose(matrixI - matrixV), (matmul(matmul(transpose(dBdy),&
     Ainv), g) - &
     matmul(matmul(transpose(B),Ainv), matmul(matmul(dAdy, Ainv),g)) + &
     matmul(matmul(transpose(B),Ainv),dgdy) ) ) - &
     matmul(matmul(matmul( transpose(matrixdVdy), transpose(B)), Ainv),g)
do j=1,NumberOfDOI2
!
!  여기서 ( ,x)는 임의의 x, (j, )는 루프상의 지점을 의미한다.
!
   Phi(j,GaussCurrent) = v(j) + temp1(j,1)
   dxPhi(j,GaussCurrent) = dvdx(j) + temp2(j,1)
   dyPhi(j,GaussCurrent) = dvdy(j) + temp3(j,1)
!   sum = sum + phi(j,GaussCurrent)
end do
!print *, GaussCurrent,"=", sum
!sum = 0.0
return
end subroutine ShapeFtnSub2


