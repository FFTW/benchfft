if HAVE_MKL
PRG32=doit32 doit32f
endif
if HAVE_MKL_ITP
PRG64=doit64 doit64f
endif
if HAVE_MKL_DFTI
PRG_DFTI=doit_dfti
endif


PRG=$(PRG32) $(PRG64) $(PRG_DFTI)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitf.c

doit32_SOURCES=doit32.c
doit32_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doit32f_SOURCES=doit32f.c
doit32f_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doit64_SOURCES=doit64.c
doit64_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

doit64f_SOURCES=doit64f.c
doit64f_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

doit_dfti_SOURCES=doit_dfti.c
doit_dfti_LDADD=-lmkl_ia32 $(LIBBENCH) @FLIBS@

include ../Makefile.common
