Subroutine Predict(NumberOfTotalNode, DisplOfNode,VelOfNode, PreVelOfNode, &
           AccelOfNoDe, TimeIncrement,InitCoordOfNode, CoordOfNode)
!=============================================================================
implicit none
!=============================================================================
INTEGER::NumberOfTotalNode, i
REAL*8::DisplOfNode(NumberOfTotalNode,2), VelOfNOde(NumberOfTotalNode,2), &
	PreVelOfNode(NumberOfTotalNode,2), AccelOfNode(NumberOfTotalNode,2), &
	TimeIncrement, InitCoordOfNode(NumberOfTotalNode,2), &
        CoordOfNode(NumberOfTotalNode,2)
!==============================================================================
do i=1, NumberOfTotalNode
   DisplOfNode(i,1) = DisplOfNode(i,1) + TimeIncrement*PreVelOfNode(i,1) + &
        0.5*(TimeIncrement**2)*AccelOfNode(i,1)
   DisplOfNode(i,2) = DisplOfNode(i,2) + TimeIncrement*PreVelOfNode(i,2) + &
        0.5*(TimeIncrement**2)*AccelOfNode(i,2)
   VelOfNode(i,1) = PreVelOfNode(i,1) + 0.5*TimeIncrement*AccelOfNode(i,1)
   VelOfNode(i,2) = PreVelOfNode(i,2) + 0.5*TimeIncrement*AccelOfNode(i,2)
end do
	PreVelOfNode = VelOfNode
!
return
end subroutine Predict
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine Correct(NumberOfTotalNode, VelOfNode, AccelOfNode, TimeIncrement)
!==============================================================================
implicit none
!==============================================================================
INTEGER::i, NumberOfTotalNode
REAL*8::VelOfNOde(NumberOfTotalNode,2), AccelOfNode(NumberOfTotalNode,2), &
     TimeIncrement
!==============================================================================
do i=1, NumberOfTotalNode
   VelOfNode(i,1) = VelOfNode(i,1) + (0.5*TimeIncrement*AccelOfNode(i,1))
   VelOfNode(i,2) = VelOfNode(i,2) + (0.5*TimeIncrement*AccelOfNode(i,2))
end do
return
end subroutine Correct
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine Acceleration(N20, NT, NumberOfCell, DisplOfNode, AccelOfNode, &
      MassOfNode, dxPhi4Node, dyPhi4Node, dxPhi4Gauss, dyPhi4Gauss, &
      NumberofDOI, NodeOfDOI, NumberOfGDOI, GaussOfDOI, Jacobian, &
      NumberOfExtForce, NodeExtF, ExtForce, CoordOfNode, &
      NumberOfMasterNode, NumberOfSlaveNode, NodeOfMaster, NodeOfSlave, &
      TimeIncrement, dxPhi, dyPhi, NodeOfDOI2, NumberOfDOI2)
!==============================================================================
implicit none
!==============================================================================
INTEGER::N20, NT, NumberOfTotalNode, NumberOfMasterNode, NumberOfSlaveNode, &
         NumberOfExtForce, i, j, k, l, NumberOfDOI(NT), &
         NodeOfDOI(N20,NT), NodeExtF(NumberOfExtForce), &
         NumberOfGDOI(NT), GaussOfDOI(N20, NT), NumberOfCell, &
         NodeOfMaster(NumberOfMasterNode), NodeOfSlave(NumberOfSlaveNode), &
         NumberOfDOI2(NumberOfCell*4), NodeOfDOI2(N20,NumberOfCell*4)
REAL*8::DisplOfNode(NT,2), AccelOfNode(NT,2), MassOfNode(NT), &
        dxPhi4Node(N20,NT), dyPhi4Node(N20,NT), dxPhi4Gauss(NT,N20), &
        dyPhi4Gauss(NT,N20), IntForceOfNode(NT,2), Jacobian(NumberOfCell), &
        ExtForceOfNode(NT,2), ExtForce(NumberOfExtForce,2), &
        ContactForce(NT,2), CoordOfNode(NT,2), TimeIncrement, stress(2,2),&
        Stress1stPK(NumberOfCell*4,2,2), dxPhi(N20,NumberOfCell*4), &
        dyPhi(N20,NumberOfCell*4)
!==============================================================================
NumberOfTotalNode = NT
IntForceOfNode = 0.0
ExtForceOfNode = 0.0
ContactForce = 0.0
!
! 1st Piola-Kirhoff stress estimation
!
do i=1, NumberOfCell*4
   Call StressPK(N20, NumberOfTotalNode, NumberOfCell, i, DisplOfNode, dxPhi, &
        dyPhi, NumberofDOI2, NodeOfDOI2, stress)
   Stress1stPK(i,1,1) = stress(1,1)
   Stress1stPK(i,1,2) = stress(1,2)
   Stress1stPK(i,2,1) = stress(2,1)
   Stress1stPK(i,2,2) = stress(2,2)
end do
!
! Internal Force estimation
!
do i = 1, NumberOfTotalNode
   do j = 1, NumberOfGDOI(i)
      !
      ! k is cell of jth gauss point
      l = GaussOfDOI(j,i)
      k = int((l - 1)/4) + 1
      IntForceOfNode(i,1) = IntForceOfNode(i,1) + &
           ((Stress1stPK(l,1,1)*dxPhi4Gauss(i,j)) + &
           (Stress1stPK(l,2,1)*dyPhi4Gauss(i,j)))*Jacobian(k)
      IntForceOfNode(i,2) = IntForceOfNode(i,2) + &
           ((Stress1stPK(l,1,2)*dxPhi4Gauss(i,j)) + &
           (Stress1stPK(l,2,2)*dyPhi4Gauss(i,j)))*Jacobian(k)
   end do
end do

!==============================================================================
!
call ExternalForce(NumberOfTotalNode, ExtForceOfNode, NumberOfExtForce, &
     NodeExtF, ExtForce)
!
! Contact treatment
!Do i=1,NumberOfContactPair
call Contact(NumberOfTotalNode, CoordOfNode, ContactForce, &
     NumberOfMasterNode, NumberOfSlaveNode, NodeOfMaster, NodeOfSlave,&
     MassOfNode, TimeIncrement)
!End do
!============================================================================
do i=1, NumberOfTotalNode
   AccelOfNode(i,1) = (ExtForceOfNode(i,1) + ContactForce(i,1) - &
        IntForceOfNode(i,1))/ MassOfNode(i)
   AccelOfNode(i,2) = (ExtForceOfNode(i,2) + ContactForce(i,2) - &
        IntForceOfNode(i,2))/ MassOfNode(i)
end do
return
end subroutine Acceleration
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine ExternalForce(NT, ExtForceOfNode, NumberOfExtForce, &
     NodeExtF, ExtForce)
!============================================================================
implicit none
!==============================================================================
INTEGER::NT, NumberOfExtForce, NodeExtF(NumberOfExtForce), i
REAL*8::ExtForceOfNode(NT,2), ExtForce(NumberOfExtForce,2)
do i=1, NumberOfExtForce
   ExtForceOfNode(NodeExtF(i),1) = ExtForce(i,1)
   ExtForceOfNode(NodeExtF(i),2) = ExtForce(i,2)
end do
return
end subroutine ExternalForce
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
Subroutine StressPK(N20, NT, NumberOfCell, GaussCurrent, &
     DisplOfNode, dxPhi, dyPhi, NumberofDOI2, NodeOfDOI2, stress1stPK)
!==============================================================================
implicit none
!
! Implementation for hyperelasticity 2000. 12. 30
! Update for 4-point integration 2001. 2. 4
!==============================================================================
REAL*8::MaterialProperty(3,5)
common/material/MaterialProperty
INTEGER::N20, NT, NumberOfTotalNode, NumberOfCell, i, j, k, NumberOfGDOI(NT), &
     GaussOfDOI(N20,NT), GaussCurrent, NodeOfDOI2(N20,NumberOfCell*4), &
     NumberOfDOI2(NumberOfCell*4)
REAL*8::stress2ndPK(2,2), &
     Gradient(2,2), DisplOfNode(NT,2), JacobiGrad, GreenStrain(2,2), &
     stress1stPK(2,2), dxPhi(N20,NumberOfCell*4), &
     dyPhi(N20,NumberOfCell*4) 
!
! hyperelastic variable
REAL*8::GradC(2,2), DetGradC, InvGradC(2,2), Inv1, Inv2, Inv3,  matrixI(2,2), &
     C1, C2, lambda
!=============================================================================
NumberOfTotalNode = NT
!
! Initialization
Gradient = 0.0
!
! Deformation gradient and right Cauch-Green deformation tensor
!
do j=1, NumberOfDOI2(GaussCurrent)
   i = NodeOfDOI2(j, GaussCurrent)
   Gradient(1,1) = Gradient(1,1) + (dxPhi(j,GaussCurrent)*DisplOfNode(i,1))
   Gradient(1,2) = Gradient(1,2) + (dyPhi(j,GaussCurrent)*DisplOfNode(i,1))
   Gradient(2,1) = Gradient(2,1) + (dxPhi(j,GaussCurrent)*DisplOfNode(i,2))
   Gradient(2,2) = Gradient(2,2) + (dyPhi(j,GaussCurrent)*DisplOfNode(i,2))
end do
Gradient(1,1) = Gradient(1,1) + 1.0
Gradient(2,2) = Gradient(2,2) + 1.0
JacobiGrad = (Gradient(1,1)*Gradient(2,2)) - (Gradient(1,2)*Gradient(2,1))
GradC = matmul(transpose(Gradient), Gradient)
!
DetGradC = (GradC(1,1)*GradC(2,2)) - (GradC(1,2)*GradC(2,1))
!
! Invariant estimation
!
Inv1 = GradC(1,1) + GradC(2,2)
Inv2 = (GradC(1,1)*GradC(2,2)) - (GradC(1,2)*GradC(2,1))
Inv3 = DetGradC
!
! Green strain 
! (1,2) is not e(1,2) but gamma(1,2), this means 2 times of e(1,2)
GreenStrain = 0.0
do i=1,2
   do j=1,2
      do k=1,2
         GreenStrain(i,j) = GreenStrain(i,j) + (Gradient(k,i)*Gradient(k,j))
      end do
   end do
end do
GreenStrain(1,1) = (GreenStrain(1,1) - 1.0)/2.0
GreenStrain(2,2) = (GreenStrain(2,2) - 1.0)/2.0
!
!
matrixI = 0.0
matrixI(1,1) = 1.0
matrixI(2,2) = 1.0
InvGradC(1,1) = GradC(2,2)/DetGradC
InvGradC(1,2) = - GradC(1,2)/DetGradC
InvGradC(2,1) = - GradC(2,1)/DetGradC
InvGradC(2,2) = GradC(1,1)/DetGradC
C1 = MaterialProperty(1,1)
C2 = MaterialProperty(2,1)
lambda = MaterialProperty(3,1)
stress2ndPK = 2.0*((C1 + Inv1*C2)*matrixI- C2*GradC - &
     ((C1*Inv3**(1./3.) + 2.0*C2*Inv3**(2./3.) - lambda*log(Inv3))*InvGradC))
!
stress1stPK = matmul(stress2ndPK, transpose(Gradient))
!
return
end subroutine StressPK




