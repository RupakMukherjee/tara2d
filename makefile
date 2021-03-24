F95 = gfortran
ILOCAL=local/include
LLOCAL=local/lib
INIFLAG=-I$(ILOCAL) -L$(LLOCAL) -lcfgio

# Definition of the Flags
INC=-I/usr/local/include
LIB=-L/usr/local/lib
FFLAGS = -w -lfftw3 -lm

EXEC=tara2d

SDIR	= src
ODIR  = src/obj
LDIR	= lib
OUTDIR = output

SRC_ 	= main.f95 derive.f95 ab.f95
OBJ_	= $(SRC_:.f95=.o)

SRC = $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

all: $(EXEC)

$(EXEC) : $(ODIR)/main.o $(OBJ)
	@echo "Linking TARA2D"
	@$(F95) $^ -o $@ $(INIFLAG) $(FFLAGS) $(INC) $(LIB)
	@echo "TARA2D is built"
	@mkdir -p $(OUTDIR)

$(ODIR)/%.o: $(SDIR)/%.f95
	@echo "Compiling $<"
	@mkdir -p $(ODIR)
	@$(F95) -c $< -o $@ $(FFLAGS) $(INC) $(LIB)

subsystem:
	@cd $(LDIR)/iniparser && $(MAKE)

clean:
	@echo "Cleaning compiled files"
	@echo "run make veryclean to remove executables"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~
	@rm -rf $(OUTDIR)

veryclean: clean
	@echo "Cleaning executables and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) clean > /dev/null 2>&1

run:
	@echo "Running TARA2D"
	./tara2d
