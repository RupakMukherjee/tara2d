F95 = gfortran
LOCAL=local
INIFLAG=-I$(LOCAL)/include -L$(LOCAL)/lib -lcfgio
C_DIR := $(shell pwd)

# Definition of the Flags
INC=-I$(LOCAL)/include
LIB=-L$(LOCAL)/lib
FFLAGS = -w -lfftw3_threads -lfftw3 -lm
OMPFLAGS = -O -fopenmp

EXEC=tara2d

SDIR	= src
ODIR  = src/obj
LDIR	= lib
OUTDIR = output

SRC_ 	= main.f95 hyd1main.f95 hyd2main.f95 hyd2mainOMP.f95 mhd2mainOMP.f95 mhd3mainOMP.f95 derive.f95 ab.f95
OBJ_	= $(SRC_:.f95=.o)

SRC = $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

all: $(EXEC)

$(EXEC) : $(ODIR)/main.o $(OBJ)
	@echo "Linking TARA2D"
	@$(F95) $^ -o $@ $(INIFLAG) $(FFLAGS) $(INC) $(LIB) $(OMPFLAGS)
	@echo "TARA2D is built"
	@mkdir -p $(OUTDIR)

$(ODIR)/%.o: $(SDIR)/%.f95
	@echo "Compiling $<"
	@mkdir -p $(ODIR)
	@$(F95) -c $< -o $@ $(FFLAGS) $(INC) $(LIB) $(OMPFLAGS)

subsystems:
	@cd $(LDIR)/iniparser && $(MAKE)
	@cd $(LDIR)/fftw && ./configure --enable-threads --with-openmp --prefix=$(C_DIR)/$(LOCAL)
	@cd $(LDIR)/fftw && $(MAKE)
	@cd $(LDIR)/fftw && $(MAKE) install
	export PATH=$(C_DIR)/$(LOCAL)/bin
	export C_INCLUDE_PATH=$(C_DIR)/$(LOCAL)/include
	export LIBRARY_PATH=$(C_DIR)/$(LOCAL)/lib

clean:
	@echo "Cleaning compiled files"
	@echo "run make veryclean to remove executables"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~
	@rm -rf $(OUTDIR)

veryclean: clean
	@echo "Cleaning executables and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) clean > /dev/null 2>&1
	@cd $(LDIR)/fftw && $(MAKE) clean > /dev/null 2>&1

run:
	@echo "Running TARA2D"
	./tara2d
