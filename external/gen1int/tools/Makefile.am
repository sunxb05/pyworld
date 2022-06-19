AUTOMAKE_OPTIONS = foreign
INCLUDES = -I../include
FCFLAGS = -O3 $(WITH_DEBUG) $(WITH_XTIME)
noinst_LIBRARIES = libtools.a
libtools_a_SOURCES = int_to_str.F90 \
                     html_log.F90   \
                     str_decode.F90
clean-generic:
	rm -f *.o *.mod
