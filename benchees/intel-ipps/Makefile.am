if HAVE_IPPS
PRG32=doit
endif

if HAVE_IPPS64
PRG64=doit64
endif

PRG=$(PRG32) $(PRG64)
INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lipps $(LIBBENCH) @FLIBS@

doit64_SOURCES=doit.c
doit64_LDADD=-lipps64 $(LIBBENCH) @FLIBS@

include ../Makefile.common
