module mesh_cut
use generalclass
implicit none


contains

!----------------------------------------------------------------
!            <SUBROUTINE> DETERMINE_RANGE
!----------------------------------------------------------------
SUBROUTINE DetermineRange(POLYGONIN,X,Y,N,M,il,ir,id,iu)
IMPLICIT NONE
TYPE(POLYGON),INTENT(IN)              :: POLYGONIN
REAL(KIND=8),INTENT(IN)               :: X(N+1),Y(M+1)
REAL(KIND=8)                          :: CLD,CLU,CLL,CLR
REAL(KIND=8),ALLOCATABLE              :: TRANS(:)
INTEGER,INTENT(OUT)                   :: il,ir,iu,id         
INTEGER                               :: i,N,M

print *,"i=0,N-1 and j=0,N-1"
stop

ALLOCATE(TRANS(POLYGONIN%NODESNUM))
do i=1,polygonin%nodesnum
   trans(i)=polygonin%nodes(i)%val(2)
ENDDO  

 CLD=MINVAL(trans,POLYGONIN%NODESNUM)
 CLU=MAXVAL(trans,POLYGONIN%NODESNUM)

do i=1,polygonin%nodesnum
   trans(i)=polygonin%nodes(i)%val(1)
ENDDO  
 CLL=MINVAL(trans,POLYGONIN%NODESNUM)
 CLR=MAXVAL(trans,POLYGONIN%NODESNUM)

DEALLOCATE(TRANS)

print *,"i=0,M-1 now ?"
stop

do i=1,M
   if ((CLD .ge. (Y(i)-eps)) .and. (CLD .LT. Y(i+1))) then
      id=i 
      exit  
   ENDIF
enddo

print *,"i=0,M-1 now ?"
stop

do i=1,M
   if (CLU.GT. Y(i) .and. CLU.LE.(Y(i+1)+eps)) then
      iu=i+1
      exit
   ENDIF
!print*,iu
enddo

print *,"i=0,N-1 now ?"
stop

do i=1,N
   if (CLL.ge. (x(i)-eps) .and. CLL .LT. X(i+1) )  then
      il=i
      exit
   ENDIF
enddo

print *,"i=0,N-1 now ?"
stop

do i=1,N
   if (CLR.GT. X(i) .and.  CLR .LE.(x(i+1)+eps)) then
      ir=i+1
      exit
   ENDIF
enddo


END SUBROUTINE DetermineRange

!-------------------------------------------------------------------
!                    MESH_CUT   
!-------------------------------------------------------------------
SUBROUTINE HorizontalMeshCut(polygonIn,X,Y,NX,MY,polygontemp,id,iu)
IMPLICIT NONE

TYPE(POLYGON)                         :: POLYGONIN
TYPE(POLYGON),ALLOCATABLE,Intent(out) :: POLYGONTEMP(:)
TYPE(POLYGON)                         :: PLGL,PLGR,PLGRTEMP
INTEGER                  ,INTENT(IN)  :: NX,MY
REAL(KIND=8)             ,INTENT(IN)  :: X(NX+1),Y(MY+1)
INTEGER                               :: RY
INTEGER                               :: i,j,k,S
INTEGER                               :: il,ir,iu,id


CALL DetermineRange(POLYGONIN,X,Y,NX,MY,il,ir,id,iu) 


RY= iu-id

ALLOCATE(POLYGONTEMP(RY))

IF(iu.eq.id+1) THEN
  CALL polygoncopy(10,POLYGONIN,POLYGONTEMP(RY))
  CALL PLGDEL(POLYGONIN)
ELSE
     CALL polygoncopy(11,PolygonIn,PLGRTEMP)
     CALL PLGDEL(POLYGONIN)
     DO i=id+1,iu-1
        CALL CUTPOLYGON_H(PLGRTEMP,Y(i),PLGL,PLGR)
        CALL polygoncopy(12,PLGL,POLYGONTEMP(i-id))
        CALL PLGDEL(PLGL)
        CALL polygoncopy(13,PLGR,PLGRTEMP)
        CALL PLGDEL(PLGR)
     ENDDO
     CALL polygoncopy(14,PLGRTEMP,POLYGONTEMP(RY))
     CALL PLGDEL(PLGRTEMP)
ENDIF
 

END SUBROUTINE HorizontalMeshCut


!----------------------------------------------------------
SUBROUTINE VerticalMeshCut(PolygonIn, X,Y,NX,MY,il,ir,PolygonOut)

TYPE(POLYGON)                         :: PolygonIn

TYPE(POLYGON),INTENT(OUT),ALLOCATABLE :: PolygonOut(:)
TYPE(POLYGON)                         :: POLYGONTEMP,PLGL,PLGR,PLGRTEMP

INTEGER                  ,INTENT(IN)  :: NX,MY
REAL(KIND=8)             ,INTENT(IN)  :: X(NX+1),Y(MY+1)
INTEGER                               :: il,ir,id,iu
INTEGER                               :: RX
INTEGER                               :: i,j,K


  CALL DETERMINERANGE(PolygonIn,X,Y,NX,MY,il,ir,id,iu)
  RX= ir-il
  ALLOCATE(PolygonOut(RX))
  IF(ir.eq.il+1) THEN
    CALL polygoncopy(15,PolygonIn,polygonout(1)) 
    CALL PLGDEL(POLYGONIN)
  ELSE      
    CALL polygoncopy(16,POLYGONIN,PLGRTEMP)
    CALL PLGDEL(POLYGONIN)
    do j=il+1,ir-1
       CALL CUTPOLYGON_V(PLGRTEMP,X(j),PLGL,PLGR)

       CALL polygoncopy(17,PLGL,PolygonOut(j-il))
       CALL PLGDEL(PLGL)
       CALL polygoncopy(18,PLGR,PLGRTEMP)
       CALL PLGDEL(PLGR)
    enddo
    CALL polygoncopy(19,PLGRTEMP,PolygonOut(RX))
    CALL PLGDEL(PLGRTEMP)
  ENDIF

END SUBROUTINE VerticalMeshCut

!////////////////////////////////////////////////////////////////


