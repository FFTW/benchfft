PRG=doit doit_fht

INCLUDES = $(INCLBENCH)

doit_SOURCES = doit.c quad2fft.c
doit_LDADD = $(LIBBENCH)

doit_fht_SOURCES = doit-fht.c b-fht.c
doit_fht_LDADD = $(LIBBENCH)
doit-fht.o: doit.c

include ../Makefile.common
