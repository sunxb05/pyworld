# Host program could include this file by setting LIB_GEN1INT_PATH
IF("${LIB_GEN1INT_PATH}" STREQUAL "")
  SET(LIB_GEN1INT_PATH ${PROJECT_SOURCE_DIR})
ENDIF()

# Source codes of Gen1Int library
SET(Gen1Int_LIB_SRCS 
    ${LIB_GEN1INT_PATH}/src/xkind.F90
    ${LIB_GEN1INT_PATH}/src/dump_info.F90
    ${LIB_GEN1INT_PATH}/src/error_stop.F90
    ${LIB_GEN1INT_PATH}/src/xtimer.F90
    ${LIB_GEN1INT_PATH}/src/tools/norm_contr_cgto.F90
    ${LIB_GEN1INT_PATH}/src/tools/norm_contr_sgto.F90
    ${LIB_GEN1INT_PATH}/src/tools/reorder_ints.F90
    ${LIB_GEN1INT_PATH}/src/tools/trace_ints.F90
    ${LIB_GEN1INT_PATH}/src/tools/get_address_list.F90
    ${LIB_GEN1INT_PATH}/src/basic/binom_coeff.F90
    ${LIB_GEN1INT_PATH}/src/basic/const_contr_ints.F90
    ${LIB_GEN1INT_PATH}/src/basic/hgto_to_cgto.F90
    ${LIB_GEN1INT_PATH}/src/basic/hgto_to_sgto.F90
    ${LIB_GEN1INT_PATH}/src/basic/next_permutation.F90
    ${LIB_GEN1INT_PATH}/src/basic/shell_scatter.F90
    ${LIB_GEN1INT_PATH}/src/basic/sort_cents.F90
    ${LIB_GEN1INT_PATH}/src/auxfun/aux_boys_vec.F90
    ${LIB_GEN1INT_PATH}/src/geom/geom_total.F90
    ${LIB_GEN1INT_PATH}/src/geom/geom_part_zero.F90
    ${LIB_GEN1INT_PATH}/src/geom/geom_part_one.F90
    ${LIB_GEN1INT_PATH}/src/mag/hgto_to_lcgto.F90
    ${LIB_GEN1INT_PATH}/src/mag/london_mom_hgto.F90
    ${LIB_GEN1INT_PATH}/src/carmom/carmom_deriv.F90
    ${LIB_GEN1INT_PATH}/src/carmom/carmom_hbra.F90
    ${LIB_GEN1INT_PATH}/src/carmom/carmom_hrr_ket.F90
    ${LIB_GEN1INT_PATH}/src/carmom/carmom_moment.F90
    ${LIB_GEN1INT_PATH}/src/carmom/prim_hgto_carmom.F90
    ${LIB_GEN1INT_PATH}/src/carmom/contr_cgto_carmom.F90
    ${LIB_GEN1INT_PATH}/src/carmom/contr_sgto_carmom.F90
    ${LIB_GEN1INT_PATH}/src/carmom/contr_csgto_carmom.F90
    ${LIB_GEN1INT_PATH}/src/delta/delta_geom.F90
    ${LIB_GEN1INT_PATH}/src/delta/delta_hket.F90
    ${LIB_GEN1INT_PATH}/src/delta/delta_moment.F90
    ${LIB_GEN1INT_PATH}/src/delta/prim_hgto_delta.F90
    ${LIB_GEN1INT_PATH}/src/delta/contr_cgto_delta.F90
    ${LIB_GEN1INT_PATH}/src/delta/contr_sgto_delta.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/nucpot_geom.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/nucpot_hket.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/nucpot_hbra.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/prim_hgto_nucpot.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/contr_cgto_nucpot.F90
    ${LIB_GEN1INT_PATH}/src/nucpot/contr_sgto_nucpot.F90
    ${LIB_GEN1INT_PATH}/src/gaupot/gaupot_geom.F90
    ${LIB_GEN1INT_PATH}/src/gaupot/prim_hgto_gaupot.F90
    ${LIB_GEN1INT_PATH}/src/gaupot/contr_cgto_gaupot.F90
    ${LIB_GEN1INT_PATH}/src/gaupot/contr_sgto_gaupot.F90
    ${LIB_GEN1INT_PATH}/src/odist/prim_hgto_odist.F90
    ${LIB_GEN1INT_PATH}/src/odist/contr_cgto_odist.F90
    ${LIB_GEN1INT_PATH}/src/odist/contr_sgto_odist.F90
    ${LIB_GEN1INT_PATH}/src/value/prim_hgto_value.F90
    ${LIB_GEN1INT_PATH}/src/value/const_contr_gto.F90
    ${LIB_GEN1INT_PATH}/src/value/contr_cgto_value.F90
    ${LIB_GEN1INT_PATH}/src/value/contr_sgto_value.F90)

# We may not need F90 module sometimes, for instance, using Gen1Int in a C
# code program?
# F90 module (we need tools/str_decode.F90 in Dalton)
IF(NOT DISABLE_F90_MODULE)
  SET(Gen1Int_LIB_SRCS
      ${Gen1Int_LIB_SRCS}
      ${LIB_GEN1INT_PATH}/src/london_ao.F90
      ${LIB_GEN1INT_PATH}/src/gen1int_geom.F90
      ${LIB_GEN1INT_PATH}/src/gen1int_carmom.F90
      ${LIB_GEN1INT_PATH}/src/gen1int_nucpot.F90
      ${LIB_GEN1INT_PATH}/src/gen1int_onehamil.F90
      ${LIB_GEN1INT_PATH}/src/gen1int_gaupot.F90
      ${LIB_GEN1INT_PATH}/src/gen1int.F90
      ${LIB_GEN1INT_PATH}/tools/str_decode.F90)
  ADD_DEFINITIONS(-DBUILD_F90_MODULE)
ENDIF()
