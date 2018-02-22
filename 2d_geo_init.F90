MODULE GeneralCLASS

IMPLICIT NONE


real(kind=8),parameter :: eps = 1e-8
real(kind=8),parameter :: tol = 1e-6
real(kind=8),parameter :: pi  = 4.0d0 * atan(1.0d0)

!//////////////////////////////////////
!    TYPE DEFINE
!/////////////////////////////////////
TYPE::Node
  REAL(KIND=8)       :: VAL(2)
  TYPE(Node),POINTER :: Next
END TYPE Node

TYPE::POINTS
    REAL(KIND=8)       :: VAL(2)
END TYPE

TYPE:: SIDE
    TYPE(POINTS)       :: PT(2)
END TYPE 

type:: IFSEG
    integer             :: flag
    TYPE(POINTS)        :: PT(2)
END TYPE

TYPE::polygon
    INTEGER                  :: NodesNum
!    INTEGER                  :: TAG
    TYPE(POINTS)             :: CENTER
    TYPE(POINTS),ALLOCATABLE :: Nodes(:)
!    INTEGER                  :: nsig(4)         ! (1 for pos   -1 for neg)
    TYPE(SIDE)               :: SIDES(4)
    TYPE(POINTS)             :: side_cen(4)
!    INTEGER                  :: csig(4)   
    TYPE(points)             :: centroid      
END TYPE

TYPE:: LINE
    INTEGER                :: TAG
    REAL(KIND=8)           :: VAL(3)
END TYPE

TYPE:: PLG_LIST
    INTEGER                       :: TAG
    INTEGER                       :: NUM
    TYPE(POLYGON),ALLOCATABLE     :: PLG(:) 
END TYPE


CONTAINS
!-----------------------------------------------
!   <SUBROUTINE> polygoncopy
!----------------------------------------------
subroutine polygoncopy(p_id,plga,plgb)             
implicit none
integer, intent(in) :: p_id
TYPE(POLYGON),INTENT(IN)   ::PLGA
TYPE(POLYGON),INTENT(OUT)  ::PLGB
INTEGER                    :: i

if(plga%nodesnum .le. 2) then
  print *,"invalid polygon in polygon copy nodesum=",plga%nodesnum
  print *,"p_id= ",p_id
  stop
endif 

PLGB%NODESNUM=PLGA%NODESNUM

if(plgb%nodesnum .ne. 0) then
allocate(PLGB%nodes(PLGB%nodesnum))

 do i=1,plgb%nodesnum
   plgb%nodes(i)%VAL=plga%nodes(i)%VAL
 enddo

endif

end subroutine polygoncopy


!------------------------------------------------
!     <subroutine> POLYGON DELETE
!------------------------------------------------
SUBROUTINE PLGDEL(plga)
IMPLICIT NONE

TYPE(POLYGON)    :: PLGA

if (plga%nodesnum .ne. 0) then
   PLGA%NODESNUM=0
   DEALLOCATE (PLGA%NODES)
endif

plga%nodesnum = 0

END SUBROUTINE PLGDEL
!---------------------------------------------------
!   <subroutine> output polygon
!-------------------------------------------------- 
subroutine outputplg(plg)
implicit none

type(polygon)        :: plg
integer              :: i

if (plg%nodesnum .eq. 0) then
   print *, "0 polygon"
else
    do i=1,plg%nodesnum
       print *, i, plg%nodes(i)%val
    enddo
endif

end subroutine

subroutine PinL(mx,my,alpha,X,VOUT)
implicit none

REAL(KIND=8),INTENT(IN)   :: mx,my,alpha
TYPE(POINTS)              :: x
REAL(KIND=8)              :: VOUT

  VOUT= mx*x%val(1)+my*x%val(2)-alpha

end subroutine PinL
!----------------------------------------------------------
subroutine LxL(a,b,mx,my,alpha,z)
implicit none

type(points),intent(in)  :: a,b
real(kind=8),intent(in)  :: mx,my,alpha
real(kind=8)             :: slope

type(points),intent(out) :: z

slope = 0.0d0

if(abs(a%val(1) - b%val(1)) .lt. eps) then
   z%val(1) = a%val(1)
   z%val(2) = -(mx*z%val(1) - alpha)/my
elseif(abs(a%val(2) - b%val(2)) .lt. eps) then
   z%val(2) = a%val(2)
   z%val(1) = -(my*z%val(2)-alpha)/mx
else
   slope = (b%val(2)-a%val(2))/(b%val(1)-a%val(1))
   z%val(1) = (my*(slope*a%val(1)-a%val(2))+alpha)/(mx+my*slope)
   z%val(2) = slope*(z%val(1)-a%val(1))+a%val(2)
endif

end subroutine LxL

!------------------------------------------------
!    <SUBROUTINE> FINDINTERSECTION_VERTICAL
!------------------------------------------------
Subroutine FindIntersection(x,y,clc,zv)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x,y
TYPE(POINTS),INTENT(OUT) :: zv
REAL(KIND=8),INTENT(IN):: clc

 zv%VAL(1)=clc
 zv%VAL(2)=(clc-y%VAL(1))/(x%VAL(1)-y%VAL(1))*x%VAL(2)+&
          &(clc-x%VAL(1))/(y%VAL(1)-x%VAL(1))*y%VAL(2)

return

END Subroutine FindIntersection
!------------------------------------------------
!     <SUBROUTINE> FINDINTERSECTION_HORIZONTAL
!------------------------------------------------
Subroutine FindIntersection_H(x,y,clc,zh)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: x,y
TYPE(POINTS),INTENT(OUT) :: zh
REAL(KIND=8),INTENT(IN):: clc

 zh%VAL(2)=clc
 zh%VAL(1)=(clc-y%VAL(2))/(x%VAL(2)-y%VAL(2))*x%VAL(1)+&
          &(clc-x%VAL(2))/(y%VAL(2)-x%VAL(2))*y%VAL(1)

return

END Subroutine FindIntersection_H
!----------------------------------------------------
!     <SUBROUTINE> DETERMINESIDE_VERTICAL
!----------------------------------------------------

SUBROUTINE DETERMINESIDE_V(a,b,clc,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)  :: a,b
TYPE(POINTS)             :: Z
Real(kind=8),INTENT(IN)  :: clc
TYPE(Node),POINTER       :: PL,PR
!TYPE(Node),POINTER       :: PLHEAD,PRHEAD

IF(((a%VAL(1).GE. clc) .and. (b%VAL(1).GT. clc)).OR.&
     &((a%VAL(1).GT.clc).and. (b%VAL(1) .GE. clc))) THEN
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
ELSEIF (((a%VAL(1).LE. clc) .and. (b%VAL(1).LT. clc)).OR.&
     &((a%VAL(1).LT.clc).and. (b%VAL(1) .LE. clc))) THEN
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

ELSEIF((a%VAL(1).GT. clc) .and. (b%VAL(1).LT. clc)) THEN
  CALL FINDINTERSECTION(a,b,clc,z)
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
 

ELSEIF((a%VAL(1).LT.clc ).and. (b%VAL(1).GT.clc)) THEN
  CALL FINDINTERSECTION(a,b,clc,z)
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

ELSEIF( abs(a%VAL(1) - clc) < eps  .and. abs( b%VAL(1) - clc) < eps )  THEN
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
   !print *, 'ERROR5'

ENDIF
END SUBROUTINE DetermineSide_V

!###################################################
!#######SUBROUTINE HORIZONTAL_DETERMINESIDE
SUBROUTINE DetermineSide_H(a,b,CLC,PL,PR)
IMPLICIT NONE

TYPE(POINTS),INTENT(IN)   :: a,b
TYPE(POINTS)              :: Z
Real(kind=8),INTENT(IN)   :: clc
TYPE(Node),POINTER        :: PL,PR

IF(((a%VAL(2).GE. clc) .and. (b%VAL(2).GT. clc)).OR.&
     &((a%VAL(2).GT.clc).and. (b%VAL(2) .GE. clc))) THEN
   IF( maxval(abs(a%VAL - PR%VAL)) < eps  ) Then

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
ELSEIF (((a%VAL(2).LE. clc) .and. (b%VAL(2).LT. clc)).OR.&
     &((a%VAL(2).LT.clc).and. (b%VAL(2) .LE. clc))) THEN

   IF( maxval(abs(a%VAL - PL%VAL)) < eps  ) Then
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

ELSEIF((a%VAL(2).GT. clc) .and. (b%VAL(2).LT. clc)) THEN
  CALL FINDINTERSECTION_H(a,b,clc,z)
   IF( maxval(abs(a%VAL - PR%VAL)) < eps  ) Then

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
 

ELSEIF((a%VAL(2).LT.clc .and. b%VAL(2).GT.clc)) THEN
  CALL FINDINTERSECTION_H(a,b,clc,z)
   IF( maxval(abs(a%VAL - PL%VAL)) < eps  ) Then

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

ELSEIF( abs(a%VAL(2) - clc) < eps .and. abs(b%VAL(2) - clc) < eps ) THEN
    IF(a%VAL(1).GT.b%VAL(1))THEN
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=a%VAL
       ALLOCATE(PL%NEXT)
       PL=>PL%NEXT
       PL%VAL=b%VAL    

    ELSEIF(a%VAL(1).LT. b%VAL(1))THEN
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=a%VAL
       ALLOCATE(PR%NEXT)
       PR=>PR%NEXT
       PR%VAL=b%VAL    
    ENDIF

ELSE
   !print *, 'ERROR2'
ENDIF

END SUBROUTINE DetermineSide_H



END MODULE
