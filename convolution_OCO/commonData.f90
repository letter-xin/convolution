MODULE commonData

PUBLIC :: LinTerp, LinTerpDouble


CONTAINS


   SUBROUTINE LinTerp(Y1, Y2, Yout, XX, x)
   ! This routine performs a linear interpolation
   ! Input: 
   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
   !        XX     !! array of length 2 containing concentrations to interpolate over.
   !        x      !! Concentration to calc absco at.

   ! Output:
   !        Yout   !! The interpolated absco.

      IMPLICIT NONE
      REAL(4), INTENT(IN) :: Y1, Y2
      REAL(4), INTENT(OUT) :: Yout
      REAL(8), INTENT(IN) :: XX(2),x      ! X(1) = lower value; X(2) = upper value
      REAL(8) :: diff1, diff2

      diff1 = x - XX(1)
      diff2 = XX(2) - XX(1)
      Yout = Y1 + diff1/diff2*(Y2 - Y1)
   END SUBROUTINE

   SUBROUTINE LinTerpDouble(Y1, Y2, Yout, X1, X2, x)
   ! This routine performs a linear interpolation
   ! Input: 
   !        Y1     !! wvnmst must be the same for Y1 and Y2 so they are of the same length.
   !        Y2     !! Otherwise, you will be interpolating between abscos of different wavenumbers.
   !        XX     !! array of length 2 containing concentrations to interpolate over.
   !        x      !! Concentration to calc absco at.

   ! Output:
   !        Yout   !! The interpolated absco.

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: Y1, Y2
      REAL(8), INTENT(OUT) :: Yout
      REAL(8), INTENT(IN) :: X1,X2,x      ! X(1) = lower value; X(2) = upper value
      REAL(8) :: diff1, diff2

      diff1 = x - X1
      diff2 = X2 - X1
      Yout = Y1 + diff1/diff2*(Y2 - Y1)
   END SUBROUTINE


END MODULE commonData

