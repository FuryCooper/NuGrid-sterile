EXEC   = Readxxx

OBJS   = vars.o param.o read_snapshot.o 

INCL   = Makefile

# module load mvapich2-1.2-sdr-gnu/4.1.2
SYSTYPE="zzc"

FFTW_INCL = -I/home/rui/library/fftw3/include
FFTW_LIBS = -L/home/rui/library/fftw3/lib

CC       =   mpicc        # sets the C-compiler (default)

FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(FFTW_LIBS) 
CFLAGS =   $(FFTW_INCL) 

$(EXEC): $(OBJS) 
	$(CC)  $(OBJS) -lm $(LIBS) -lfftw3 -g  $(CFLAGS)  -o  $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



