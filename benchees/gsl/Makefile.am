if HAVE_GSL
PRG=doit doit-dit doit-dif doit-real
else
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST=doit2.c

doit_SOURCES=doit.c
doit_LDADD=-lgsl -lgslcblas $(LIBBENCH)

doit_dit_SOURCES=doit-dit.c
doit_dit_LDADD=-lgsl -lgslcblas $(LIBBENCH)

doit_dif_SOURCES=doit-dif.c
doit_dif_LDADD=-lgsl -lgslcblas $(LIBBENCH)

doit_real_SOURCES=doit-real.c
doit_real_LDADD=-lgsl -lgslcblas $(LIBBENCH)

include ../Makefile.common
