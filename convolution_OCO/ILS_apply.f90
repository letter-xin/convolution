PROGRAM ILS_apply
USE ILS_table
USE commonData
IMPLICIT NONE


REAL(8)::hres(2,12442)
REAL(8)::obs_spec(2,1016)
REAL(8)::plist(652)
REAL(8)::ils_delta(200,1016)
REAL(8)::ils_res(200,1016)
REAL(8), ALLOCATABLE:: lres(:,:)
REAL(4)::t1,t2
INTEGER:: i, j, k
INTEGER::plist_int(652)
INTEGER:: ierr

OPEN(UNIT = 21, FILE= 'input/hres.dat', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF
DO i=1,2
    READ(21,*)(hres(i,j),j=1,12442)
END DO
OPEN(UNIT = 22, FILE= 'input/plist.dat', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF
READ(22,*)(plist(i),i=1,652)

OPEN(UNIT = 23, FILE= 'input/obs_spec.dat', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF
DO i=1,2
    READ(23,*)(obs_spec(i,j),j=1,1016)
END DO





OPEN(UNIT = 31, FILE= 'output/lres.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF

OPEN(UNIT = 41, FILE= 'input/ils_delta.dat', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF
DO i=1,200
    READ(41,*)(ils_delta(i,j),j=1,1016)
END DO
OPEN(UNIT = 42, FILE= 'input/ils_res.dat', STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierr)
IF(ierr /= 0) THEN
	WRITE(*,*) 'Error! Could not open output file!'
ENDIF
DO i=1,200
    READ(42,*)(ils_res(i,j),j=1,1016)
END DO



DO i=1,652
    plist_int(i) = INT(plist(i)+1.E-6)
END DO
call cpu_time(t1)
ALLOCATE(lres(2,652))
DO i = 1,652
write(*,*) "Now calculating n=", i
    lres(1,i)=obs_spec(1,plist_int(i))
    CALL ils_conv(obs_spec(1,:),plist_int(i),ils_delta,ils_res,hres,lres(2,i))
END DO
call cpu_time(t2)
WRITE(*,*) "Time=",t2-t1
Do j=1,652
    WRITE(31,"(*(E12.6,:,','))") (lres(i,j) , i= 1 , 2)
END DO


close(21)
close(22)
close(23)
close(31)
close(41)
close(42)


END PROGRAM ILS_apply