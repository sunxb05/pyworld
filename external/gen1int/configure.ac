#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT([gen1int], [0.2.1], [bin.gao@uit.no])
AC_CONFIG_SRCDIR([src/geom/geom_total.F90])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([gen1int], [0.2.1])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC([pgf90 ifort gfortran])
AC_PROG_F77([pgf90 ifort gfortran])
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Fortran 90 module
AC_ARG_ENABLE(fmodule,
  [--enable-fmodule  enable Fortran 90 module],
  [case "${enableval}" in
     yes) fmodule=true ;;
     no) fmodule=false ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --enable-fmodule) ;;
   esac],
  [fmodule=false])
AM_CONDITIONAL(WITH_F_MOD, test x$fmodule = xtrue)

# Debug mode
AC_SUBST([WITH_DEBUG])
AC_ARG_WITH([debug],
  [AS_HELP_STRING([--with-debug],
    [enable debug mode (default: disabled)])],
  [WITH_DEBUG="-g -DDEBUG"],
  [WITH_DEBUG=""])

# Print time usage information
AC_SUBST([WITH_XTIME])
AC_ARG_WITH([time],
  [AS_HELP_STRING([--with-time],
    [print time usage information (default: disabled)])],
  [WITH_XTIME="-DXTIME"],
  [WITH_XTIME=""])

AC_CONFIG_FILES([Makefile src/Makefile tools/Makefile test_f90/Makefile])
AC_OUTPUT
