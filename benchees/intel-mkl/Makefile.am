if HAVE_MKL_DEF
PRG_DEF=doit32_def doit32f_def
endif
if HAVE_MKL_P3
PRG_P3=doit32_p3 doit32f_p3
endif
if HAVE_MKL_P4
PRG_P4=doit32_p4 doit32f_p4
endif
if HAVE_MKL_ITP
PRG_ITP=doit64_itp doit64f_itp
endif

PRG=$(PRG_DEF) $(PRG_P3) $(PRG_P4) $(PRG_ITP)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitf.c

doit32_def_SOURCES=doit32_def.c
doit32_def_LDADD=-lmkl_def $(LIBBENCH) @FLIBS@

doit32f_def_SOURCES=doit32f_def.c
doit32f_def_LDADD=-lmkl_def $(LIBBENCH) @FLIBS@

doit32_p3_SOURCES=doit32_p3.c have-sse.c
doit32_p3_LDADD=-lmkl_p3 $(LIBBENCH) @FLIBS@

doit32f_p3_SOURCES=doit32f_p3.c have-sse.c
doit32f_p3_LDADD=-lmkl_p3 $(LIBBENCH) @FLIBS@

doit32_p4_SOURCES=doit32_p4.c have-sse2.c
doit32_p4_LDADD=-lmkl_p4 $(LIBBENCH) @FLIBS@

doit32f_p4_SOURCES=doit32f_p4.c have-sse2.c
doit32f_p4_LDADD=-lmkl_p4 $(LIBBENCH) @FLIBS@

doit64_itp_SOURCES=doit64_itp.c
doit64_itp_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

doit64f_itp_SOURCES=doit64f_itp.c
doit64f_itp_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

include ../Makefile.common
