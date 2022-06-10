# Makefile for test of LBFGSB

SHELL = cmd
FC = gfortran
FLINKER = $(FC)
ARCH = ar
ARCHFLAGS = -rsc
OPTS = -s

.PHONY: all
all: LBFGSB_wrapper.mod libLBFGSB_wrapper.a

.PHONY: test
test: test.exe
	.\$<

test.exe: test.obj all
	@echo Linking program $@ against $^ ...
	$(FLINKER) -o $@ $< -L . -l LBFGSB_wrapper -static $(OPTS)

test.obj: test.F90 LBFGSB_wrapper.mod
	@echo Compiling $@ from $< ...
	$(FC) -o $@ -c $< $(OPTS) -cpp -I .

libLBFGSB_wrapper.a: LBFGSB_wrapper.obj LBFGSB.obj linpack.obj blas.obj
	@echo Generating archive $@ from $^ ...
	$(ARCH) $(ARCHFLAGS) $@ $^

LBFGSB_wrapper.mod: LBFGSB_wrapper.obj

LBFGSB_wrapper.obj: LBFGSB_wrapper.F90
	@echo Compiling $@ from $< ...
	$(FC) -o $@ -c $< $(OPTS) -cpp

%.obj: %.f
	@echo Compiling $@ from $< ...
	$(FC) -o $@ -c $< $(OPTS)

.PHONY: clean
clean: clean_test
	-del /q *.a *.mod 2> NUL

.PHONY: clean_test
clean_test:
	-del /q *.exe *.obj iterate.dat 2> NUL

