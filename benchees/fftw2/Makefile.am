if SINGLE
if HAVE_SFFTW2
PRG=doitm doite doitm_nd doite_nd
FFTW_LIBS=@LSFFTW@
endif
else
if HAVE_FFTW2
PRG=doitm doite doitm_nd doite_nd
FFTW_LIBS=@LFFTW@
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitnd.c doit-doc.h

doitm_SOURCES=doitm.c
doitm_LDADD=$(FFTW_LIBS) $(LIBBENCH) @FLIBS@

doite_SOURCES=doite.c
doite_LDADD=$(FFTW_LIBS) $(LIBBENCH) @FLIBS@

doitm_nd_SOURCES=doitm_nd.c
doitm_nd_LDADD=$(FFTW_LIBS) $(LIBBENCH) @FLIBS@

doite_nd_SOURCES=doite_nd.c
doite_nd_LDADD=$(FFTW_LIBS) $(LIBBENCH) @FLIBS@

include ../Makefile.common
