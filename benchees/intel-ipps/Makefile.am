if HAVE_IPPS
PRG32=doit
endif

if HAVE_IPPSEM64T
PRGEM64T=doitem64t
endif

if HAVE_IPPS64
PRG64=doit64
endif

PRG=$(PRG32) $(PRGEM64T) $(PRG64)
INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lipps $(LIBBENCH) @FLIBS@

doitem64t_SOURCES=doit.c
doitem64t_LDADD=-lippsem64t $(LIBBENCH) @FLIBS@

doit64_SOURCES=doit.c
doit64_LDADD=-lipps64 $(LIBBENCH) @FLIBS@

include ../Makefile.common
