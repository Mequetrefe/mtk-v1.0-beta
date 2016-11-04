# Makefile for the INCLUDE folder.

include ../makefile_inc

all: check

CHECK = mtk.gch \
\
	mtk_roots.gch \
	mtk_constants.gch \
	mtk_enums.gch \
\
	mtk_tool_manager.gch \
	mtk_scalar_vector_manager.gch \
\
	mtk_dense_matrix.gch \
	mtk_1d_array_dense_stencil_matrix.gch \
	mtk_2d_array_dense_matrix.gch \
	mtk_crs_sparse_matrix.gch \
	mtk_ccs_sparse_matrix.gch \
	mtk_dok_sparse_matrix.gch \
\
	mtk_grid_manager.gch \
	mtk_1d_node.gch \
	mtk_grid_1D.gch \
	mtk_lr_1d_uniform_nodal_grid.gch \
	mtk_lr_1d_uniform_staggered_grid.gch \
	mtk_lr_1d_uniform_source_grid.gch \
	mtk_lr_1d_uniform_known_solution_holder_grid.gch \
	mtk_lr_1d_nonuniform_curvilinear_staggered_grid.gch \
\
	mtk_1d_cgm_dirichlet_operator.gch \
	mtk_1d_cgm_neumann_operator.gch \
	mtk_1d_cgm_divergence.gch \
	mtk_1d_cgm_gradient.gch \
	mtk_1d_cgm_laplacian.gch \
	mtk_2d_cgm_gradient.gch \
	mtk_2d_cgm_divergence.gch \
\
	mtk_steady_1d_problem.gch \
	mtk_system.gch \
\
	mtk_1d_plotter.gch \
	mtk_1d_plot_properties.gch

check: $(CHECK)
	@echo 'Checking the header files... Done!'

%.gch: %.h
	$(CC) $(CCFLAGS) $<

clean:
	rm -f *~ *.gch
