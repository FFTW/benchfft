if SINGLE
PRG=doit
endif

PKG=vBigDSP
TARFILE=$(PKG)-20030322.tgz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

doit_SOURCES=doit.c
nodist_doit_SOURCES = $(PKG)/vBigDSP.c
doit_CFLAGS=$(AM_CFLAGS) @ALTIVEC_CFLAGS@ -Dvalloc=bench_malloc
doit_LDADD=$(LIBBENCH)

$(PKG)/vBigDSP.c: extract
doit.o: extract

include ../Makefile.common
