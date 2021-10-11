if DOUBLE
PRG=doit
noinst_LIBRARIES=libpocketfft.a
endif

INCLUDES = $(INCLBENCH)

libpocketfft_a_SOURCES=pocketfft.c

doit_SOURCES=doit.c
doit_LDADD=libpocketfft.a $(LIBBENCH) @FLIBS@

include ../Makefile.common

