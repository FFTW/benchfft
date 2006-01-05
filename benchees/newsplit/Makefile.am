PRG=doit
EXTRA_PRG=doit_unsym doit_oldsplit

EXTRA_DIST = README doit.c

INCLUDES = $(INCLBENCH)

SPLIT = split.h split.c

doit_SOURCES = $(SPLIT) newsplit.c doit-newsplit.c
doit_LDADD=$(LIBBENCH)
doit-newsplit.c: doit.c split.h

doit_unsym_SOURCES = $(SPLIT) newsplit_unsym.c doit-unsym.c
doit_unsym_LDADD=$(LIBBENCH)
doit-unsym.c: doit.c split.h

doit_oldsplit_SOURCES = $(SPLIT) oldsplit.c doit-oldsplit.c
doit_oldsplit_LDADD=$(LIBBENCH)
doit-oldsplit.c: doit.c split.h

include ../Makefile.common
