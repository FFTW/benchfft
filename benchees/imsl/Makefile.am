if HAVE_IMSL
PRG = doit
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES = doit.c
doit_LDADD = @LIBIMSL@ $(LIBBENCH) @FLIBS@

include ../Makefile.common
