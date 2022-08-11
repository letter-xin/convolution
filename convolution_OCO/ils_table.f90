MODULE ILS_table
USE commonData

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE ils_conv(obs_spec,plist_int,ils_delta,ils_res,hres,lres)
IMPLICIT NONE

! Passed variables
REAL(8),INTENT(IN) :: obs_spec(1016)
INTEGER,INTENT(IN) :: plist_int
REAL(8),INTENT(IN) :: ils_delta(200,1016)
REAL(8),INTENT(IN) :: ils_res(200,1016)
REAL(8),INTENT(IN) :: hres(2,12442)
REAL(8),INTENT(OUT) :: lres

! FFTW parameters
REAL(8) :: lambda,delta

REAL(8) :: ils_delta_hres(1,2)
REAL(8) :: ils_res_hres(1,2)
REAL(8) :: bound_low,bound_high
REAL(8) :: delta_x,y,sumsum
REAL(8) :: ils_table,sum_part
! Other parameters
INTEGER	:: i,k,j,i1,i2,j1,j2,k1,k2,k3

!!!!!!!!!!!!!!!!!!!!!!!!!! End variable declaration !!!!!!!!!!!!!!!!!!!!!!!!!!!



lambda=obs_spec(plist_int) 

if ((lambda+ils_delta(1,plist_int))>=hres(1,1))THEN
	bound_low=lambda+ils_delta(1,plist_int)
else
	bound_low=hres(1,1) 
end if
if (lambda+ils_delta(200,plist_int)<=hres(1,12442)) THEN
	bound_high=lambda+ils_delta(200,plist_int) 
else
	bound_high=hres(1,12442) 
end if
k1=1
k2=1
DO i=1,12442
	iF (hres(1,i)-bound_low<0.d0) Then
		k1=k1+1
	END IF
	iF (hres(1,i)-bound_high<0.d0) Then
		k2=k2+1
	END IF
END DO

i1=k1
i2=k2-1 

sumsum=0.0d0
lres=0.0d0


DO i=i1,i2
	delta=hres(1,i)-lambda 
	k3=1
	if (ils_delta(199,plist_int)-delta<=0) THEN
		j1=INT(199.d0+1.E-6) 
		j2=INT(200.d0+1.E-6) 
	else
		DO j=1,200
			IF (ils_delta(j,plist_int)-delta<0) Then
				k3=k3+1
			END IF
		END DO
	END if
	j1=k3
	j2=j1+1
	ils_delta_hres(1,1)=ils_delta(j1,plist_int) 
	ils_delta_hres(1,2)=ils_delta(j2,plist_int) 
	ils_res_hres(1,1)=ils_res(j1,plist_int) 
	ils_res_hres(1,2)=ils_res(j2,plist_int) 
	
	CALL LinTerpDouble(ils_res_hres(1,1),ils_res_hres(1,2),ils_table,ils_delta_hres(1,1),ils_delta_hres(1,2),delta)

	sum_part=hres(2,i)*ils_table*(hres(1,i+1)-hres(1,i)) 
	lres=lres+sum_part 
	delta_x=hres(1,i+1)-hres(1,i)
	y=ils_table
	sumsum=sumsum+delta_x*y
end DO
lres=lres/sumsum



END SUBROUTINE ils_conv


END MODULE ILS_table
