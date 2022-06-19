AUTOMAKE_OPTIONS = foreign
INCLUDES = -I../include
FCFLAGS = -O3 $(WITH_DEBUG) $(WITH_XTIME)
noinst_LIBRARIES = libgen1int.a
libgen1int_a_SOURCES = dump_info.F90                 \
                       error_stop.F90                \
                       xtimer.F90                    \
                       tools/norm_contr_cgto.F90     \
                       tools/norm_contr_sgto.F90     \
                       tools/reorder_ints.F90        \
                       tools/trace_ints.F90          \
                       tools/get_address_list.F90    \
                       basic/binom_coeff.F90         \
                       basic/const_contr_ints.F90    \
                       basic/hgto_to_cgto.F90        \
                       basic/hgto_to_sgto.F90        \
                       basic/next_permutation.F90    \
                       basic/shell_scatter.F90       \
                       basic/sort_cents.F90          \
                       auxfun/aux_boys_vec.F90       \
                       geom/geom_total.F90           \
                       geom/geom_part_zero.F90       \
                       geom/geom_part_one.F90        \
                       mag/hgto_to_lcgto.F90         \
                       mag/london_mom_hgto.F90       \
                       carmom/carmom_deriv.F90       \
                       carmom/carmom_hbra.F90        \
                       carmom/carmom_hrr_ket.F90     \
                       carmom/carmom_moment.F90      \
                       carmom/prim_hgto_carmom.F90   \
                       carmom/contr_cgto_carmom.F90  \
                       carmom/contr_sgto_carmom.F90  \
                       carmom/contr_csgto_carmom.F90 \
                       delta/delta_geom.F90          \
                       delta/delta_hket.F90          \
                       delta/delta_moment.F90        \
                       delta/prim_hgto_delta.F90     \
                       delta/contr_cgto_delta.F90    \
                       delta/contr_sgto_delta.F90    \
                       nucpot/nucpot_geom.F90        \
                       nucpot/nucpot_hket.F90        \
                       nucpot/nucpot_hbra.F90        \
                       nucpot/prim_hgto_nucpot.F90   \
                       nucpot/contr_cgto_nucpot.F90  \
                       nucpot/contr_sgto_nucpot.F90  \
                       gaupot/gaupot_geom.F90        \
                       gaupot/prim_hgto_gaupot.F90   \
                       gaupot/contr_cgto_gaupot.F90  \
                       gaupot/contr_sgto_gaupot.F90  \
                       odist/prim_hgto_odist.F90     \
                       odist/contr_cgto_odist.F90    \
                       odist/contr_sgto_odist.F90    \
                       value/prim_hgto_value.F90     \
                       value/const_contr_gto.F90     \
                       value/contr_cgto_value.F90    \
                       value/contr_sgto_value.F90
if WITH_F_MOD
libgen1int_a_SOURCES += gen1int_geom.F90     \
                        gen1int_carmom.F90   \
                        gen1int_nucpot.F90   \
                        gen1int_onehamil.F90 \
                        gen1int_gaupot.F90   \
                        gen1int.F90
endif
clean-generic:
	rm -f *.o *.mod
