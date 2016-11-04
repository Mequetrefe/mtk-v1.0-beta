# Makefile for the EXAMPLES folder.

include ../makefile_inc

all: ex01 ex02 ex03 ex04 ex05 ex06 ex07 ex08 ex09 ex10 ex11 ex12 ex13 ex14

EXE = \
	ex01-1d_poisson_driver plot.gpl solution.dat \
	ex02-crs_sparse_matrices_driver \
	ex03-ccs_sparse_matrices_driver \
	ex04-solution_holder_grid_heun_method \
	ex05-1d_diff_react \
	ex06-solution_holder_grid_access \
	ex07-stencil_matrices_mutators_accesors \
	ex08-scalar_vector_product \
	ex09-curvilinear_hyperbolic_time_dependent \
	ex10-accuracy_1d_poisson_driver \
	ex11-1d_poisson_driver_batista_castillo \
	ex12-accuracy_1d_poisson_driver_batista_castillo \
	ex13-dok_sparse_matrices_driver \
	ex14-2d_cgm_operators_driver

# Compiling statements:

%.o: %.cc $(INCLUDE)/mtk.h
	$(CC) $(CCFLAGS) $<

# Linking statements:
ex01: ex01-1d_poisson_driver.o ex01-1d_poisson_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex01-1d_poisson_driver.o $(LIBS) -o ex01-1d_poisson_driver
	@echo ' '

ex02: ex02-crs_sparse_matrices_driver.o ex02-crs_sparse_matrices_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex02-crs_sparse_matrices_driver.o $(LIBS) -o ex02-crs_sparse_matrices_driver
	@echo ' '

ex03: ex03-ccs_sparse_matrices_driver.o ex03-ccs_sparse_matrices_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex03-ccs_sparse_matrices_driver.o $(LIBS) -o ex03-ccs_sparse_matrices_driver
	@echo ' '

ex04: ex04-solution_holder_grid_heun_method.o ex04-solution_holder_grid_heun_method.cc $(INCLUDE)/mtk.h
	$(LINKER) ex04-solution_holder_grid_heun_method.o $(LIBS) -o ex04-solution_holder_grid_heun_method
	@echo ' '

ex05: ex05-1d_diff_react.o ex05-1d_diff_react.cc $(INCLUDE)/mtk.h
	$(LINKER) ex05-1d_diff_react.o $(LIBS) -o ex05-1d_diff_react
	@echo ' '

ex06: ex06-solution_holder_grid_access.o ex06-solution_holder_grid_access.cc $(INCLUDE)/mtk.h
	$(LINKER) ex06-solution_holder_grid_access.o $(LIBS) -o ex06-solution_holder_grid_access
	@echo ' '

ex07: ex07-stencil_matrices_mutators_accesors.o ex07-stencil_matrices_mutators_accesors.cc $(INCLUDE)/mtk.h
	$(LINKER) ex07-stencil_matrices_mutators_accesors.o $(LIBS) -o ex07-stencil_matrices_mutators_accesors
	@echo ' '

ex08: ex08-scalar_vector_product.o ex08-scalar_vector_product.cc $(INCLUDE)/mtk.h
	$(LINKER) ex08-scalar_vector_product.o $(LIBS) -o ex08-scalar_vector_product
	@echo ' '

ex09: ex09-curvilinear_hyperbolic_time_dependent.o ex09-curvilinear_hyperbolic_time_dependent.cc $(INCLUDE)/mtk.h
	$(LINKER) ex09-curvilinear_hyperbolic_time_dependent.o $(LIBS) -o ex09-curvilinear_hyperbolic_time_dependent
	@echo ' '

ex10: ex10-accuracy_1d_poisson_driver.o ex10-accuracy_1d_poisson_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex10-accuracy_1d_poisson_driver.o $(LIBS) -o ex10-accuracy_1d_poisson_driver
	@echo ' '

ex11: ex11-1d_poisson_driver_batista_castillo.o ex11-1d_poisson_driver_batista_castillo.cc $(INCLUDE)/mtk.h
	$(LINKER) ex11-1d_poisson_driver_batista_castillo.o $(LIBS) -o ex11-1d_poisson_driver_batista_castillo
	@echo ' '

ex12: ex12-accuracy_1d_poisson_driver_batista_castillo.o ex12-accuracy_1d_poisson_driver_batista_castillo.cc $(INCLUDE)/mtk.h
	$(LINKER) ex12-accuracy_1d_poisson_driver_batista_castillo.o $(LIBS) -o ex12-accuracy_1d_poisson_driver_batista_castillo
	@echo ' '

ex13: ex13-dok_sparse_matrices_driver.o ex13-dok_sparse_matrices_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex13-dok_sparse_matrices_driver.o $(LIBS) -o ex13-dok_sparse_matrices_driver
	@echo ' '

ex14: ex14-2d_cgm_operators_driver.o ex14-2d_cgm_operators_driver.cc $(INCLUDE)/mtk.h
	$(LINKER) ex14-2d_cgm_operators_driver.o $(LIBS) -o ex14-2d_cgm_operators_driver
	@echo ' '

clean:
	rm -f *~ *.o $(EXE) *.gif *.png
