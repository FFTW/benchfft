if HAVE_MKL
PRG_MKL=doit doitf
endif
if HAVE_MKL_DFTI
PRG_DFTI=doit_dfti
endif


PRG=$(PRG_MKL) $(PRG_DFTI)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitf.c

doit_SOURCES=doit.c
doit_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doitf_SOURCES=doitf.c
doitf_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doit_dfti_SOURCES=doit_dfti.c
doit_dfti_LDADD=-lmkl $(LIBBENCH) @FLIBS@

include ../Makefile.common
