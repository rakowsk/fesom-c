!=================================================================================|
! Set up sigma and more, to be continued....                                                       !
!								                                                  !
! case(1): SIGMA LEVELS                                                           !
! Sigma levels are determined by a formula of                                     !
!                       sigma(k)=-[(k-1)/(kb-1)]^pow_sig                          !
!    pow_sig=1: uniform sigma layers                                              !
!    pow_sig=2: layers satisfying a parabolic function with high                  !
!               vertical resolution near the surface and bottom.                  !
!    pow_sig can be used any real number                                          !
!									                                              !
! case(2): GENERAL VERTICAL LEVEL                                                 !
! Vertical levels are determined by the formula                                   !     
!                                                                                 !
!          tanh[(lev_bot+lev_sur)*(nsigma-k)/(nsigma-1)-lev_bot]+tanh(lev_bot)    !                   
! sigma(k)= ------------------------------------------------------------------ - 1!                 
!                      tanh(lev_bot) + tanh(lev_sur)                              !      
!                                                                                 !
!                                                                                 !
! case(3): CONSTANT LAYER TRANSFORMATION                                          !
! It is not a real sigma, this is useful approach for the sediment solutions      !  
! Near bottom or/and near surface you have layers with constant depth.            !
! If the depth is too large, than you have case of ind_sigma=2. This model        !
! is useful for some task in frame sediment movement simulation.                  !
! For this case you need prescribe the number of these constant layers            !
! and its thickness. Also you prescribe the critical depth, from which you will   !
! shift to usual sigma distribution (ind_sigma=2).If the depth is less than sum   !
! of the depths of 'constant layers', than the uniform sigma layers               !
! will be prescribed. Useful refernce is:                                         !
! Pietrzak, Jakobson, Burchard, Vested, Petersen, 2002. A three-dimensional       !
! hydrostatic model for coastal and ocean modelling using a generalised topography!  
! following co-ordinate system. Ocean Modelling 4, 173-205                        !
!=================================================================================|
       
SUBROUTINE SET_SIGMA         

!=================================================================================|
USE o_MESH
USE o_ARRAYS
USE o_PARAM

IMPLICIT NONE
   INTEGER :: K,KK,lb,ls
   INTEGER :: I
   REAL(WP):: X1,X2,X3                 
   REAL(WP):: DR,CL,CU, CB, CS, D_ref_min 
   Real(kind=WP), allocatable   :: ZKBw(:), ZKSw(:)

!==============================================================================|

   SELECT CASE(ind_sigma)

   CASE(1)
allocate(sigma(nsigma))
   IF(pow_sig > 1 .AND. MOD(nsigma,2) == 0)THEN
     WRITE(*,*) 'nsigma shoulde be an odd number,stop ....'
     STOP
   END IF

   IF(pow_sig == 1)THEN
     DO K=1,nsigma
       sigma(K) = -((K-1)/DFLOAT(nsigma-1))**pow_sig + 1.0_WP
     END DO
   ELSE
     DO K=1,(nsigma+1)/2
       sigma(K) = -((K-1)/DFLOAT((nsigma+1)/2-1))**pow_sig*0.5_WP +1.0_WP
     END DO

     DO K=(nsigma+1)/2+1,nsigma
       sigma(K) = ((nsigma-K)/DFLOAT((nsigma+1)/2-1))**pow_sig*0.5_WP 
     END DO
   END IF

!==============================================================================|   
   CASE(2)
   allocate(sigma(nsigma))
   ! Sigma 
   DO K=1,nsigma
     X1=TANH((lev_bot+lev_sur)*(nsigma-K)/(nsigma-1) - lev_bot)
     X2=TANH(lev_bot)
     X3=X2+TANH(lev_sur)
     sigma(K)=(X1+X2)/X3
    END DO

!==============================================================================| 
   CASE(3)
   allocate(Zl(nsigma,nod2D))
   
 !  IF((sum(KSw)+sum(KBw))>D_ref_min)THEN
 !     WRITE(*,*) '(sum(KSw)+sum(KBw)) should be less then or equal to D_ref_min....'
 !    STOP
 !  END IF

   D_ref_min=sum(KSw)+sum(KBw)
  
   DO I=1,nod2D
    IF(depth(I) < D_ref_min)THEN	 
	 DO K=1,nsigma
	   Zl(K,I)=-((K-1)/DFLOAT(nsigma-1)) + 1.0_WP
     END DO
  ELSE
	 IF(depth(I) > D_ref_max)THEN
	 lev_bot=1.0_WP+size(KBw)/nsigma
	 lev_sur=1.0_WP+size(KSw)/nsigma
	 
       DO K=1,nsigma
     X1=TANH((lev_bot+lev_sur)*(nsigma-K)/(nsigma-1) - lev_bot)
     X2=TANH(lev_bot)
     X3=X2+TANH(lev_sur)
     Zl(K,I)=(X1+X2)/X3
    END DO
	
	ELSE 
    lb=size(KBw)
	ls=size(KSw)
	allocate(ZKSw(ls),ZKBw(lb)) 
	   CS=-sum(KSw)/depth(I)
       CB=sum(KBw)/depth(I)-1.0_WP
       DR=(CB-CS)/(nsigma-ls-lb-1)

       DO K=1,ls
         ZKSw(K)=KSw(K)/depth(I)
       END DO
       DO K=1,lb
         ZKBw(K)=KBw(K)/depth(I)
       END DO

       Zl(1,I)=1.0_WP

       DO K=2,ls+1
         Zl(K,I)=Zl(K-1,I)-ZKSw(K-1)
       END DO

       DO K=ls+2,nsigma-lb
         Zl(K,I)=Zl(K-1,I)+DR
       END DO

       KK=0
       DO K=nsigma-lb+1,nsigma
         KK=KK+1
         Zl(K,I)=Zl(K-1,I)-ZKBw(KK)
       END DO
	   Zl(nsigma,I)=0.0_WP !can be very small, but non-zero value, after the loop above  
	   deallocate(ZKSw,ZKBw)
     END IF
	 END IF
   END DO   
   END SELECT   

   RETURN
   END SUBROUTINE SET_SIGMA
!==============================================================================|
