!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
subroutine Contact(NumberOfTotalNode, CoordOfNode, ContactForce, &
     NumberOfMasterNode, NumberOfSlaveNode, NodeOfMaster, NodeOfSlave, &
     MassOfNode, TimeIncrement)
!==============================================================================
implicit none
!==============================================================================
integer::NumberOfMasterNode, NumberOfSlaveNode, &
     NodeOfMaster(NumberOfMasterNode), NodeOfSlave(NumberOfSlaveNode), &
     NFlag(NumberOfSlaveNode), NumberOfTotalNode, &
     NumberOfPair(NumberOfSlaveNode),&
     NPair(NumberOfSlaveNode, NumberOfMasterNode-1),&
     i,j,k,l,N1,N2,N3, N4(1)
real*8::CoordOfNode(NumberOfTotalNode,2), ContactForce(NumberOfTotalNode,2), &
	LengthOfSurface(NumberOfMasterNode-1), Length4,&
	Area(NumberOfMasterNode-1,NumberOfSlaveNode), &
	xmin(NumberOfMasterNode-1), xmax(NumberOfMasterNode-1), &
        ymin(NumberOfMasterNode-1), ymax(NumberOfMasterNode-1), &
        a(NumberOfMasterNode-1), b(NumberOfMasterNode-1), c, &
        Delta(NumberOfMasterNode-1), theta1, theta2, &
        MassOfNode(NumberOfTotalNode), TimeIncrement, DOP, normalx, normaly
!==============================================================================
!  Caution! : argument of Length~( ), Area( ), Delta is not node number.
!             it is only numbering of surface series
!==============================================================================
!
! Check loop - step #1
do i=1, NumberOfSlaveNode
   NumberOfPair(i) = 0
   NFlag(i) = 0
   N3 = NodeOfSlave(i)       ! N3 is node of slave contact
   do j=1,NumberOfMasterNode-1
      N1 = NodeOfMaster(j)     ! N2 is first node of master contact
      N2 = NodeOfMaster(j+1)    ! N3 is second node of master contact
      if ( i == 1) then   ! estimate length of master surface at only initial stage
         LengthOfSurface(j) = sqrt(( CoordOfNode(N1,1) - &
              CoordOfNode(N2,1) )**2 + &
              ( CoordOfNOde(N1,2) - CoordOfNode(N2,2) )**2 )
         Length4 = LengthOfSurface(i)/4.0       ! bias for bucket size
         !
         ! max/min for bucket range about master surface j
         xmin(j) = min(CoordOfNode(N1,1) - Length4, CoordOfNode(N2,1) - &
              Length4 )
         xmax(j) = max(CoordOfNode(N1,1) + Length4, CoordOfNode(N2,1) + &
              Length4 )
         ymin(j) = min(CoordOfNode(N1,2) - Length4, CoordOfNode(N2,2) - &
              Length4 )
         ymax(j) = max(CoordOfNode(N1,2) + Length4, CoordOfNode(N2,2) + &
              Length4 )
      end if
      if (((CoordOfNode(N3,1)>xmin(j)).AND.(CoordOfNode(N3,1)<xmax(j))) .AND. &
           ((CoordOfNode(N3,2)>ymin(j)).AND.(CoordOfNode(N3,2)<ymax(j)))) then
!
! If slave node exists in bucket range, estimate triangle
         Area(j,i) = (CoordOfNode(N2,1)*CoordOfNode(N3,2) + &
              CoordOfNode(N1,1)*CoordOfNode(N2,2) + &
              CoordOfNode(N3,1)*CoordOfNode(N1,2) - &
              CoordOfNode(N1,2)*CoordOfNode(N2,1) - &
              CoordOfNode(N2,2)*CoordOfNode(N3,1) - &
              CoordOfNode(N3,2)*CoordOfNode(N1,1) ) * 0.5
         if (Area(j,i) > 0.0 ) then   ! watch out that "Area(Master, Slave)"
! 
! slave is seized and penetrated
! if NFlag = 1, slave node attacks master surface, 
!               contact force will be estimated
! if NFlag = 0/-1, slave node avoids master surface, 
!                  contact routine will be neglected
            NFlag(i) = 1
!
! NumberOfPair(i) is the number of possible contact pairs for slave i
            NumberOfPair(i) = NumberOfPair(i) + 1
!
! NPair(M,N) is series of master node, M is series of slave node, 
! N is number of the possible contact pair
            NPair(i, NumberOfPair(i)) = j
            k = k + 1
         else
! 
! slave is seized but not penetrated
!NFlag(i) = -1
!go to 1001
         end if
      end if
   end do
1001 end do
!
! estimate loop
Do i = 1, NumberOfSlaveNode
Delta = 10.0   ! initialize as just big value.
if(NFlag(i) == 1) then
   N3 = NodeOfSlave(i)  ! N3 is slave node
   Do j=1, NumberOfPair(i)
      k = NPair(i,j) 
      N1 = NodeOfMaster(k)
      N2 = NodeOfMaster(k+1)
      ! is'nt c is length of surface?
      !c = sqrt((CoordOfNode(N2,1) - CoordOfNode(N1,1))**2 + (CoordOfNode(N2,2) - &
      !		  CoordOfNode(N1,2))**2)
      c = LengthOfSurface(k)
      b(k) = sqrt((CoordOfNode(N3,1) - CoordOfNode(N1,1))**2 + &
           (CoordOfNode(N3,2) - CoordOfNode(N1,2))**2)
      a(k) = sqrt((CoordOfNode(N2,1) - CoordOfNode(N3,1))**2 + &
           (CoordOfNode(N2,2) - CoordOfNode(N3,2))**2)
      ! 
      ! estimate difference b/w slave-first master and slave-second master
      ! pair which consists of minimum difference will be used as contact pair
      Delta(k) = abs(acos((a(k)**2-c**2-b(k)**2)/(2*c*b(k))) - &
           acos((b(k)**2-c**2- a(k)**2)/(2*c*a(k))))
   end do
   ! 
   ! l is master node of true contact pair for slave node i
   N4 = minloc(Delta)
   l = N4(1)
   !
   DOP = 2.*area(l,i)/LengthOfSurface(l)
   N2 = NodeOfMaster(l)    ! N2 is first node of master surface
   N1 = NodeOfMaster(l+1)     ! N1 is second node of master surface
   ! 
   ! Treatment of slave node
   normalx = (CoordOfNode(N2,2) - CoordOfNode(N1,2))/LengthOfSurface(l)
   normaly = (CoordOfNode(N1,1) - CoordOfNode(N2,1))/LengthOfSurface(l)
   ContactForce(N3,1) = -1.0*normalx*DOP*MassOfNode(N3)/(TimeIncrement**2)
   CoordOfNode(N3,1) = CoordOfNode(N3,1) - DOP*normalx
   ContactForce(N3,2) = -1.0*normaly*DOP*MassOfNode(N3)/(TimeIncrement**2)
   CoordOfNode(N3,2) = CoordOfNode(N3,2) - DOP*normaly
   !
   ! Treatment of master nodes
   if(abs(DOP-b(l)) < 0.0001) then 
      theta1 = 1.0e9
   else
      theta1 = tan(asin(DOP/b(l)))
   end if
   if(abs(DOP-a(l)) < 0.0001) then 
      theta2 = 1.0e9
   else
      theta2 = tan(asin(DOP/a(l)))
   end if
   !
   ! ¼öÁ¤ ¿ä¸Á 11.19
   ContactForce(N2,1) = ContactForce(N2,1) - ContacTForce(N3,1) * theta1 / &
        (theta1 + theta2)
   ContactForce(N2,2) = ContactForce(N2,2) - ContacTForce(N3,2) * theta1 / &
        (theta1 + theta2)
   ContactForce(N1,1) = ContactForce(N1,1) - ContacTForce(N3,1) * theta2 / &
        (theta1 + theta2)
   ContactForce(N1,2) = ContactForce(N1,2) - ContacTForce(N3,2) * theta2 / &
        (theta1 + theta2)
end if
end do
return
end subroutine Contact



