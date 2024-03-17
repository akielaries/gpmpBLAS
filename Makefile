# Compiler
FC 			= gfortran

# Flags
FFLAGS 		= -Wall -Wextra
BUILDDIR 	= build

# Files
SRCS 		= dgemm.f xerbla.f lsame.f ex_dgemm0.f90
OBJS 		= $(patsubst %.f,$(BUILDDIR)/%.o,$(patsubst %.f90,$(BUILDDIR)/%.o,$(SRCS)))

# Executable
EXEC 		= ex_dgemm

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(BUILDDIR)/%.o: %.f90
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.f
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) $(EXEC)

.PHONY: all clean

