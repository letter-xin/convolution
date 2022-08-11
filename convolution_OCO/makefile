FILES = commonData.f90 ils_table.f90 ILS_apply.f90
EXE = ILS

F95 = ifort 
F95FLAGS = -FR -w95
F95FLAGS1 = -FR -w95 
F95FLAGS5 = -FR -w -cm -CB -traceback -g -assume byterecl


LDFLAGS =-r8 -assume byterecl

OBJS  = commonData.o ils_table.o ILS_apply.o

all: $(EXE)

$(EXE): $(OBJS)
	$(F95) $(LDFLAGS) -o $(EXE) $(OBJS) -L /public/home/lcc-rt/opt/fftw/.libs/ -lfftw3 -lm


%.o: %.f90
	$(F95) $(F95FLAGS5) -c $*.f90

clean:
	-rm -f $(OBJS) *.mod *.o
