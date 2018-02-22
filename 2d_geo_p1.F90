MODULE l_inter_p
USE generalclass

implicit none


contains


!-----------------------------------------------------
!         <SUBROUTINE> CUTPOLYGON_VERTICAL
!-------------------------------------------------------
SUBROUTINE CUTPOLYGON_V(PLG,CL,PLGL,PLGR)      
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
REAL(KIND=8),INTENT(IN)   :: CL
TYPE(Node),POINTER        :: PLHEAD,PL,PR,PRHEAD
INTEGER                   :: i,j,NUML,NUMR,N
REAL(KIND=8)              :: P(2),Q(2)

ALLOCATE(PRHEAD)
ALLOCATE(PLHEAD)

PLhead%val(1) = 1.0e+8
PLhead%val(2) = 1.0e+8
PRhead%val(1) = 1.0e+8
PRhead%val(2) = 1.0e+8

PL => PLhead
PR => PRhead

DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL DETERMINESIDE_V(PLG%Nodes(i),PLG%Nodes(i+1),CL,PL,PR)
ENDDO
   CALL DETERMINESIDE_V(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),CL,PL,PR)

nullify(pl%next,pr%next)


PR=>PRHEAD
PR=>PR%NEXT
P=PR%VAL
N=1
do while(ASSOCIATED(PR%NEXT))
   PR=>PR%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P - PR%VAL)) < eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF



PL=>PLHEAD
PL=>PL%NEXT
P=PL%VAL
N=1

do while(ASSOCIATED(PL%NEXT))

   PL=>PL%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P-PL%VAL)) < eps  ) then 

  PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

ALLOCATE(PLGL%NODES(PLGL%NODESNUM),PLGR%NODES(PLGR%NODESNUM))

PL=>PLHEAD
do i=1,PLGL%NODESNUM
   PL=>PL%NEXT
   PLGL%NODES(i)%VAL=PL%VAL
ENDDO


PR=>PRHEAD
do i=1,PLGR%NODESNUM
   PR=>PR%NEXT
   PLGR%NODES(i)%VAL=PR%VAL
ENDDO

NULLIFY(PL%NEXT,PR%NEXT)
NULLIFY(PLHEAD,PRHEAD)
CALL PLGDEL(PLG)
END SUBROUTINE CUTPOLYGON_V

!---------------------------------------------------
!       <SUBROUTINE> CUTPOLYGON_HORIZONTAL
!---------------------------------------------------
SUBROUTINE CUTPOLYGON_H(PLG,CL,PLGL,PLGR) !SUBROUTINE CUTPOLYGON_H
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
REAL(KIND=8),INTENT(IN)   :: CL
TYPE(Node),POINTER        :: PLHEAD1,PL1,PR1,PRHEAD1
INTEGER                   :: i,j,NUML,NUMR,N
REAL(KIND=8)              :: P(2),Q(2)

ALLOCATE(PRHEAD1)
ALLOCATE(PLHEAD1)

PLhead1%val(1) = 1.0e+8
PLhead1%val(2) = 1.0e+8
PRhead1%val(1) = 1.0e+8
PRhead1%val(2) = 1.0e+8

PL1 => PLhead1
PR1 => PRhead1

DO i=1,PLG%NodesNum-1                
   CALL DETERMINESIDE_H(PLG%Nodes(i),PLG%Nodes(i+1),CL,PL1,PR1)
ENDDO
   CALL DETERMINESIDE_H(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),CL,PL1,PR1)

nullify(pl1%next,pr1%next)


PR1=>PRHEAD1
PR1=>PR1%NEXT
P=PR1%VAL
N=1
do while(ASSOCIATED(PR1%NEXT))
   PR1=>PR1%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P - PR1%VAL)) < eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF


PL1=>PLHEAD1
PL1=>PL1%NEXT
P=PL1%VAL
N=1
do while(ASSOCIATED(PL1%NEXT))
   PL1=>PL1%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P - PL1%VAL)) < eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

ALLOCATE(PLGL%NODES(PLGL%NODESNUM),PLGR%NODES(PLGR%NODESNUM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PL1=>PLHEAD1
do i=1,PLGL%NODESNUM
   PL1=>PL1%NEXT
   PLGL%NODES(i)%VAL=PL1%VAL
ENDDO


PR1=>PRHEAD1
do i=1,PLGR%NODESNUM
   PR1=>PR1%NEXT
   PLGR%NODES(i)%VAL=PR1%VAL
ENDDO
                                
NULLIFY(PL1%NEXT,PR1%NEXT)
NULLIFY(PLHEAD1,PRHEAD1)

CALL PLGDEL(PLG)
END SUBROUTINE CUTPOLYGON_H
!--------------------------------------------------------

!--------------------------------------------------------
SUBROUTINE DETERMINESIDE(a,b,mx,my,alpha,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
real(kind=8),intent(in)  :: mx,my,alpha
real(kind=8)             :: nx,ny,nalpha


TYPE(POINTS)             :: Z
TYPE(Node),POINTER       :: PL,PR
real(kind=8)             :: aval,bval


nx = mx
ny = my
nalpha=alpha


if (nx .lt. 0.0d0) then
   nx = -1.0d0*nx
   ny = -1.0d0*ny
   nalpha = -1.0d0*nalpha
endif

if(nx .eq. 0.0d0  .and. ny .lt. 0.0d0) then
   ny = -1.0d0*ny
   nalpha = -1.0d0*nalpha
endif



call PinL(nx,ny,nalpha,a,aval)
call PinL(nx,ny,nalpha,b,bval)


IF(((aval.GE. 0.0d0) .and. (bval.GT. 0.0d0)).OR.&
     &((aval .GT. 0.0d0).and. (bval .GE. 0.0d0))) THEN
   IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then

     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((aval .LE. 0.0d0) .and. (bval.LT. 0.0d0)).OR.&
     &((aval .LT. 0.0d0).and. (bval .LE. 0.0d0))) THEN
    IF( maxval(abs( a%VAL - PL%VAL)) < eps ) Then

     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
    ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
    ENDIF 

ELSEIF((aval.GT. 0.0d0) .and. (bval .LT. 0.0d0)) THEN
  CALL LxL(a,b,mx,my,alpha,z)
    IF( maxval(abs( a%VAL - PR%VAL)) < eps ) Then

    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ELSE
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=a%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((aval.LT. 0.0d0 ).and. (bval.GT. 0.0d0)) THEN
  CALL LxL(a,b,mx,my,alpha,z)
    IF( maxval(abs( a%VAL - PL%VAL)) < eps ) Then

    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ELSE
    ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=Z%VAL
    ENDIF
 
   ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=Z%VAL
    ALLOCATE(PR%NEXT)
    PR=>PR%NEXT
    PR%VAL=B%VAL
!------------------------------------------------------------
ELSEIF( aval .eq. 0.0d0  .and. bval .eq. 0.0d0 )  THEN
  if (ny .eq. 0.0d0 .and. nx .ne. 0.0d0) then
    
   IF(a%val(2) .GT. b%val(2))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
   ELSEIF(a%val(2) .LT. b%val(2))THEN
     IF( maxval ( abs ( a%VAL - PL%VAL) ) < eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
     else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
     endif  
   else
     write(*,*) 'cutpolygon:ERROR1'
   endif


 elseif(nx .eq. 0.0d0 .and. ny .ne. 0.0d0 ) then
     IF(a%val(1) .lT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .gT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
 elseif(nx .ne. 0.0d0 .and. ny .ne. 0.0d0)then
   if(nx*ny .gt. 0.0d0) then
     IF(a%val(1) .lT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .gT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
   ELSE
     IF(a%val(1) .GT. b%val(1))THEN
       IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then
                ALLOCATE(PR%NEXT)
                PR=>PR%NEXT
                PR%VAL=b%VAL
       else
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=a%VAL
            ALLOCATE(PR%NEXT)
            PR=>PR%NEXT
            PR%VAL=b%VAL    
       endif
     ELSEIF(a%val(1) .LT. b%val(1))THEN
      IF( maxval ( abs ( a%VAL - PL%VAL) ) < eps ) Then
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL 
      else 
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL
      endif  
     else
       write(*,*) 'cutpolygon: ERROR2'
     endif
   endif
 else
   print *, "cutpolygon: nx = 0, ny = 0"
   
ENDIF
!------------------------------------------------------------------
ELSE
   print *, 'ERROR3'
ENDIF

END SUBROUTINE DetermineSide


!--------------------------------------------------------
subroutine cutpolygon(a,b,alpha,PLG,PLGL,PLGR)
!----------------------------------------------------------
!   ax+by = alpha
!--------------------------------------------------------
implicit none

REAL(KIND=8),intent(in)    :: a,b,alpha
type(polygon)              :: plg
type(polygon),intent(out)  :: plgl,plgr

TYPE(Node),POINTER         :: PLHEAD1,PL1,PR1,PRHEAD1
INTEGER                    :: i,j,NUML,NUMR,N
REAL(KIND=8)               :: P(2),Q(2)


ALLOCATE(PRHEAD1)
ALLOCATE(PLHEAD1)

PLhead1%val(1) = 1.0e+8
PLhead1%val(2) = 1.0e+8
PRhead1%val(1) = 1.0e+8
PRhead1%val(2) = 1.0e+8

PL1 => PLhead1
PR1 => PRhead1


if (a .eq. 0.0d0 .and. b .eq. 0.0d0) then
WRITE(*,*) " a and b both 0"
endif

DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL DETERMINESIDE(PLG%Nodes(i),PLG%Nodes(i+1),a,b,alpha,PL1,PR1)
ENDDO
   CALL DETERMINESIDE(PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),a,b,alpha,PL1,PR1)


nullify(pl1%next,pr1%next)


PR1=>PRHEAD1
if (.not. associated(pr1%next)) then
  PLGR%NODESNUM = 0
else
PR1=>PR1%NEXT
P=PR1%VAL
N=1
do while(ASSOCIATED(PR1%NEXT))
   PR1=>PR1%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P - PR1%VAL)) < eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF

endif


PL1=>PLHEAD1
if(.not. associated(PL1%next)) then
   PLGL%NODESNUM = 0 
ELSE
PL1=>PL1%NEXT
P=PL1%VAL
N=1
do while(ASSOCIATED(PL1%NEXT))
   PL1=>PL1%NEXT
   N=N+1
ENDDO

IF( maxval(abs(P - PL1%VAL)) < eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

endif

!-------------------------------------------------
IF(PLGL%NODESNUM .ne. 0) then
ALLOCATE(PLGL%NODES(PLGL%NODESNUM))
PL1=>PLHEAD1
do i=1,PLGL%NODESNUM
   PL1=>PL1%NEXT
   PLGL%NODES(i)%VAL=PL1%VAL
ENDDO

endif

if(PLGR%NODESNUM .ne. 0) then
ALLOCATE(PLGR%NODES(PLGR%NODESNUM))
PR1=>PRHEAD1
do i=1,PLGR%NODESNUM
   PR1=>PR1%NEXT
   PLGR%NODES(i)%VAL=PR1%VAL
ENDDO

endif 


NULLIFY(PL1%NEXT,PR1%NEXT)
NULLIFY(PLHEAD1,PRHEAD1)


end subroutine cutpolygon

!//////////////////////////////////////////////////////////////
SUBROUTINE l_DetermineSide(alpha,beta,zeta,a,b,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
TYPE(POINTS)             :: Z
TYPE(NODE),POINTER       :: PL,PR
real(kind=8)             :: vala,valb
real(kind=8),intent(in)  :: alpha, beta, zeta

call cut_DIFF(alpha,beta,zeta,a,vala)
call cut_DIFF(alpha,beta,zeta,b,valb)


IF(((vala .GE. 0) .and. (valb.GT. 0)).OR.&
     &((vala.GT. 0).and. (valb .GE. 0))) THEN
   IF( maxval ( abs ( a%VAL - PR%VAL) ) < eps ) Then
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ElSE
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=> PR%NEXT
     PR%VAL=b%VAL
   ENDIF
ELSEIF (((vala.LE. 0) .and. (valb.LT. 0)).OR.&
     &((vala.LT. 0).and. (valb .LE. 0))) THEN
    IF( maxval(abs( a%VAL - PL%VAL)) < eps ) Then
      ALLOCATE(PL%NEXT)
      PL=> PL%NEXT
      PL%VAL=b%VAL
    ElSE
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=a%VAL
     ALLOCATE(PL%NEXT)
     PL=> PL%NEXT
     PL%VAL=b%VAL
   ENDIF 

ELSEIF((vala.GT. 0) .and. (valb.LT. 0)) THEN
  CALL l_FINDINTERSECTION(alpha,beta,zeta,a,b,z)                           !!!!!!!!!!!!!!!!
    IF( maxval(abs( a%VAL - PR%VAL)) < eps ) Then
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=z%VAL
    ELSE
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=a%VAL
     ALLOCATE(PR%NEXT)
     PR=>PR%NEXT
     PR%VAL=z%VAL
    ENDIF
     
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=z%VAL
    ALLOCATE(PL%NEXT)
    PL=>PL%NEXT
    PL%VAL=b%VAL
 

ELSEIF((vala .LT. 0 ).and. (valb .GT. 0)) THEN
  CALL l_FINDINTERSECTION(alpha,beta,zeta,a,b,z)
   IF( maxval(abs( a%VAL - PL%VAL)) < eps ) Then
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=Z%VAL
   ELSE
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=A%VAL
      ALLOCATE(PL%NEXT)
      PL=>PL%NEXT
      PL%VAL=Z%VAL
   ENDIF
 
   ALLOCATE(PR%NEXT)
   PR=>PR%NEXT
   PR%VAL=Z%VAL
   ALLOCATE(PR%NEXT)
   PR=>PR%NEXT
   PR%VAL=B%VAL

ELSEIF(vala .EQ. 0 .and. valb .EQ. 0) THEN
    IF(a%VAL(2).GT.b%VAL(2))THEN
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=a%VAL
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=b%VAL    
    ELSEIF(a%VAL(2).LT. b%VAL(2))THEN
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL    
    ENDIF

ELSE
   !print *, 'ERROR4'

ENDIF
END SUBROUTINE l_DetermineSide

!----------------------------------------------------------------
SUBROUTINE l_CUTPOLYGON(PLG,a,b,d,PLGL,PLGR)             
IMPLICIT NONE

TYPE(POLYGON)             :: PLG
TYPE(POLYGON),INTENT(OUT) :: PLGL,PLGR
TYPE(Node),POINTER        :: PLHEAD,PL,PR,PRHEAD
INTEGER                   :: i,j,NUML,NUMR,N
REAL(KIND=8)              :: P(2),Q(2)
REAL(KIND=8),INTENT(IN)   :: a,b,d

plg%nodesnum = 4

ALLOCATE(PRHEAD)
ALLOCATE(PLHEAD)

PLhead%val(1) = 1.0e+8
PLhead%val(2) = 1.0e+8
PRhead%val(1) = 1.0e+8
PRhead%val(2) = 1.0e+8

PL => PLhead
PR => PRhead



DO i=1,PLG%NodesNum-1                  !not allocate memory PLG!
   CALL L_DETERMINESIDE(a,b,d,PLG%Nodes(i),PLG%Nodes(i+1),PL,PR)   
ENDDO
   CALL L_DETERMINESIDE(a,b,d,PLG%Nodes(PLG%NODESNUM),PLG%Nodes(1),PL,PR)
NULLIFY(PL%NEXT,PR%NEXT)

!--------count--------------------

PR=>PRHEAD
if (.not. associated(pr%next)) then
  PLGR%NODESNUM = 0
else
PR=>PR%NEXT
P=PR%VAL
N=1
do while(ASSOCIATED(PR%NEXT))
   PR=>PR%NEXT
   N=N+1
ENDDO


IF( maxval(abs(P - PR%VAL)) < eps ) THEN
   PLGR%NODESNUM=N-1
ELSE 
   PLGR%NODESNUM=N
ENDIF

endif
!-------------------------------
PL=>PLHEAD
if(.not. associated(PL%next)) then
   PLGL%NODESNUM = 0 
ELSE
PL=>PL%NEXT
P=PL%VAL
N=1
do while(ASSOCIATED(PL%NEXT))
   PL=>PL%NEXT
   N=N+1
ENDDO


IF( maxval(abs(P - PL%VAL)) < eps ) THEN
   PLGL%NODESNUM=N-1
ELSE 
   PLGL%NODESNUM=N
ENDIF

endif
!--------------------------------
IF(PLGL%NODESNUM .NE. 0) THEN
  ALLOCATE(PLGL%NODES(PLGL%NODESNUM))
  PL=>PLHEAD
  do i=1,PLGL%NODESNUM
    PL=>PL%NEXT
    PLGL%NODES(i)%VAL=PL%VAL
  ENDDO
ENDIF

IF(PLGR%NODESNUM .NE. 0) THEN
  ALLOCATE(PLGR%NODES(PLGR%NODESNUM))
  PR=>PRHEAD
  do i=1,PLGR%NODESNUM
    PR=>PR%NEXT
    PLGR%NODES(i)%VAL=PR%VAL
  ENDDO
ENDIF

NULLIFY(PL%NEXT,PR%NEXT)
NULLIFY(PLHEAD,PRHEAD)

END SUBROUTINE l_CUTPOLYGON

Subroutine L_FindIntersection(a,b,d,x1,x2,zout)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x1,x2
TYPE(POINTS),INTENT(OUT) :: zout
REAL(KIND=8)             :: m
REAL(KIND=8),INTENT(IN)  :: a,b,d

if (x2%val(1) .ne. x1%val(1)) THEN
   m= (x2%val(2)-x1%val(2))/(x2%val(1)-x1%val(1))

   IF((b .NE. 0.0d0) .and. (a .NE. 0.0d0)) then
     zout%val(1)= (m*x1%val(1)-x1%val(2)-d/b)/(m+a/b)
     zout%val(2)= -a/b*zout%val(1)-d/b

   ELSEIF((b .EQ. 0.0d0) .and. (a .NE. 0.0d0)) then
    zout%val(1)= -d/a
    zout%val(2)= x1%val(2)+m*(-d/a-x1%val(1)) 
   ELSE 
    zout%val(2)= -d/b
    zout%val(1)= (-d/b-x1%val(2))/m+x1%val(1) 
   ENDIF
ELSE
   zout%val(1)= x1%val(1)
   zout%val(2)= (-d-a*x1%val(1))/b
   
ENDIF

END Subroutine l_FindIntersection


















end module
