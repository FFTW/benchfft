if HAVE_MKL
PRG32=doit32_ia32 doit32f_ia32
endif
if HAVE_MKL_ITP
PRG64=doit64_itp doit64f_itp
endif

PRG=$(PRG32) $(PRG64)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c doitf.c

doit32_ia32_SOURCES=doit32_ia32.c
doit32_ia32_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doit32f_ia32_SOURCES=doit32f_ia32.c
doit32f_ia32_LDADD=-lmkl $(LIBBENCH) @FLIBS@

doit64_itp_SOURCES=doit64_itp.c
doit64_itp_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

doit64f_itp_SOURCES=doit64f_itp.c
doit64f_itp_LDADD=-lmkl_itp $(LIBBENCH) @FLIBS@

include ../Makefile.common
