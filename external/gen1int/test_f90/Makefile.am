AUTOMAKE_OPTIONS = foreign
INCLUDES = -I../include -I../src -I../tools
FCFLAGS = -O3 $(WITH_DEBUG) $(WITH_XTIME)
bin_PROGRAMS = test_gen1int
test_gen1int_SOURCES = check_difference.F90               \
                       xyz_function.F90                   \
                       recur/recur_hgto.F90               \
                       recur/recur_mag.F90                \
                       recur/recur_carmom.F90             \
                       recur/recur_delta.F90              \
                       recur/recur_nucpot.F90             \
                       tools/test_norm_contr_gto.F90      \
                       tools/test_reorder_ints.F90        \
                       tools/test_trace_ints.F90          \
                       basic/test_binom_coeff.F90         \
                       basic/test_const_contr_ints.F90    \
                       basic/test_hgto_to_cgto.F90        \
                       basic/test_next_permutation.F90    \
                       basic/test_shell_scatter.F90       \
                       basic/test_sort_cents.F90          \
                       auxfun/test_aux_boys_vec.F90       \
                       geom/test_geom_total.F90           \
                       geom/test_geom_part_zero.F90       \
                       geom/test_geom_part_one.F90        \
                       mag/test_hgto_to_lcgto.F90         \
                       mag/test_london_mom_hgto.F90       \
                       carmom/test_carmom_deriv.F90       \
                       carmom/test_carmom_hbra.F90        \
                       carmom/test_carmom_hrr_ket.F90     \
                       carmom/test_carmom_moment.F90      \
                       carmom/test_prim_hgto_carmom.F90   \
                       carmom/test_contr_cgto_carmom.F90  \
                       carmom/test_contr_sgto_carmom.F90  \
                       carmom/test_contr_csgto_carmom.F90 \
                       carmom/test_contr_sgto_kinene.F90  \
                       carmom/test_contr_csgto_kinene.F90 \
                       delta/test_delta_geom.F90          \
                       delta/test_delta_hket.F90          \
                       delta/test_delta_moment.F90        \
                       delta/test_prim_hgto_delta.F90     \
                       delta/test_contr_cgto_delta.F90    \
                       delta/test_contr_sgto_delta.F90    \
                       nucpot/test_nucpot_geom.F90        \
                       nucpot/test_nucpot_hket.F90        \
                       nucpot/test_nucpot_hbra.F90        \
                       nucpot/test_prim_hgto_nucpot.F90   \
                       nucpot/test_contr_cgto_nucpot.F90  \
                       nucpot/test_contr_sgto_nucpot.F90  \
                       gaupot/test_gaupot_geom.F90        \
                       odist/test_prim_hgto_odist.F90     \
                       value/test_prim_hgto_value.F90     \
                       test_gen1int.F90
if WITH_F_MOD
FCFLAGS += -DBUILD_F90_MODULE
test_gen1int_SOURCES += f90mod/test_f90mod_sgto_kinene.F90 \
                        f90mod/test_f90mod_sgto_carmom.F90 \
                        f90mod/test_f90mod_sgto_sphmom.F90
endif
test_gen1int_LDADD = ../src/libgen1int.a ../tools/libtools.a
clean-generic:
	rm -f *.o *.mod
