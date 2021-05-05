.SUFFIXES:


SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
EXE     = $(BIN_DIR)/jdm_rmc
FC      = gfortran
FCLAGS  = -std=f2018
LDFLAGS = -L/home/joe/QUIP/build/linux_x86_64_gfortran          # LD directories
LDLIBS  = -lblas -latoms -llapack #-lgsl    # library names
INCLUDES = -I. -I/home/joe/QUIP/build/linux_x86_64_gfortran -I/p_2/Documents/MoS/RMC/mt_stream_f90-1.11 -I/p_2/Documents/MoS/RMC/mt_stream_f90-1.11/f_jump_ahead_coeff
COMPILE.f08 = $(FC) $(INCLUDES) $(FCFLAGS) $(TARGET_ARCH) -c

SRC = $(wildcard $(SRC_DIR)/*.f95)
OBJ = $(patsubst $(SRC_DIR)/%.f95, $(OBJ_DIR)/%.o,  $(SRC))
MT = /p_2/Documents/MoS/RMC/mt_stream_f90-1.11
EX_OBJ = $(MT)/mt_stream.o $(MT)/f_jump_ahead_coeff/gf2xe.o $(MT)/f_jump_ahead_coeff/f_get_coeff.o

all: $(EXE)

.PHONY: all clean

$(EXE): $(OBJ) $(EX_OBJ) | $(BIN_DIR)
	$(FC) $(FCLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@  

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f95 | $(OBJ_DIR)
	$(COMPILE.f08) $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@	

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR) # The @ disables the echoing of the command

## compilation rules
$(OBJ_DIR)/%.o $(OBJ_DIR)/%.mod $(OBJ_DIR)/%.smod: $(SRC_DIR)/%.f95
	$(COMPILE.f08) -o $(OBJ_DIR)/$*.o $<
	@touch $@

## .o -> .mod of the modules it uses
$(OBJ_DIR)/mc.o: $(OBJ_DIR)/numeric_kinds.mod
