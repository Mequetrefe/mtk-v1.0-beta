# Makefile for the SRC folder.

include ../makefile_inc

all: datastructs exectools meshesandgrids mimopers nummethods \
  input solvers outputandvis archive

DATASTRUCTS  = \
	mtk_dense_matrix.o \
	mtk_1d_array_dense_stencil_matrix.o \
	mtk_2d_array_dense_matrix.o \
	mtk_crs_sparse_matrix.o \
	mtk_ccs_sparse_matrix.o \
	mtk_dok_sparse_matrix.o

EXECTOOLS = mtk_tool_manager.o \
	mtk_scalar_vector_manager.o

MESHESANDGRIDS = \
	mtk_grid_manager.o \
	mtk_grid_1D.o \
	mtk_lr_1d_uniform_nodal_grid.o \
	mtk_lr_1d_uniform_staggered_grid.o \
	mtk_lr_1d_uniform_source_grid.o \
	mtk_lr_1d_uniform_known_solution_holder_grid.o

MIMOPERS = \
	mtk_1d_cgm_dirichlet_operator.o \
	mtk_1d_cgm_neumann_operator.o \
	mtk_1d_cgm_gradient.o \
	mtk_1d_cgm_laplacian.o \
	mtk_2d_cgm_gradient.o \
	mtk_2d_cgm_divergence.o

NUMMETHODS = mtk_system.o

INPUT = mtk_steady_1d_problem.o

SOLVERS =

OUTPUTANDVIS = mtk_1d_plotter.o

datastructs: $(DATASTRUCTS)
	@echo ' '
	@echo '----- Compiling 3rd concern (Data structures): Done!'
	@echo ' '

exectools: $(EXECTOOLS)
	@echo ' '
	@echo '----- Compiling 4th concern (Execution tools): Done!'
	@echo ' '

meshesandgrids: $(MESHESANDGRIDS)
	@echo ' '
	@echo '----- Compiling 5th concern (Meshes and grids): Done!'
	@echo ' '

mimopers: $(MIMOPERS)
	@echo ' '
	@echo '----- Compiling 6th concern (Mimetic operators): Done!'
	@echo ' '

nummethods: $(NUMMETHODS)
	@echo ' '
	@echo '----- Compiling 7th concern (Numerical methods): Done!'
	@echo ' '

input: $(INPUT)
	@echo ' '
	@echo '----- Compiling 8th concern (Input): Done!'
	@echo ' '

solvers: $(SOLVERS)
	@echo ' '
	@echo '----- Compiling 9th concern (Solvers): Done!'
	@echo ' '

outputandvis: $(OUTPUTANDVIS)
	@echo ' '
	@echo '----- Compiling 10th concern (Output and visualization): Done!'
	@echo ' '

# Compilation statements:

%.o: %.cc
	$(CC) $(CCFLAGS) $<

# Notice we create the building directory with "-p" so no prob if it exists:

archive: $(DATASTRUCTS) $(EXECTOOLS) $(MESHESANDGRIDS) $(MIMOPERS) \
  $(NUMMETHODS) $(INPUT) $(SOLVERS) $(OUTPUTANDVIS)
	$(ARCH) $(ARCHFLAGS) $(MTKLIB) \
	  $(DATASTRUCTS) \
	  $(EXECTOOLS) \
	  $(MESHESANDGRIDS) \
	  $(MIMOPERS) \
	  $(NUMMETHODS) \
	  $(INPUT) \
	  $(SOLVERS) \
	  $(OUTPUTANDVIS)
	$(RANLIB) $(MTKLIB)
	@echo ' '
	@echo '----- Library created! Check:' $(LIB)
	@echo ' '

clean:
	rm -f *~ *.o
