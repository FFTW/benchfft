if SINGLE
if HAVE_F77
PRG=doit doit_3d
noinst_LIBRARIES=libgpfa.a
endif
endif

PKG = gpfa
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libgpfa.a $(LIBBENCH) @FLIBS@ 

doit_3d_SOURCES=doit_3d.c 
doit_3d_LDADD=libgpfa.a $(LIBBENCH) @FLIBS@

libgpfa_a_SOURCES = $(PKG)/gpfa.f $(PKG)/gpfa2f.f $(PKG)/gpfa3f.f	\
$(PKG)/gpfa5f.f $(PKG)/setgpfa.f $(PKG)/gpfft3.f

include ../Makefile.common
