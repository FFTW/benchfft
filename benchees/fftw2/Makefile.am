if HAVE_FFTW2
PRG=doitm doite doitm_nd doite_nd
else
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitnd.c

doitm_SOURCES=doitm.c
doitm_LDADD=$(LIBBENCH) @FLIBS@

doite_SOURCES=doite.c
doite_LDADD=$(LIBBENCH) @FLIBS@

doitm_nd_SOURCES=doitm_nd.c
doitm_nd_LDADD=$(LIBBENCH) @FLIBS@

doite_nd_SOURCES=doite_nd.c
doite_nd_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
