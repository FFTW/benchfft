if HAVE_MKL32_DEF
PRG_DEF=doit32_def
endif
if HAVE_MKL32_P3
PRG_P3=doit32_p3
endif
if HAVE_MKL32_P4
PRG_P4=doit32_p4
endif

PRG=$(PRG_DEF) $(PRG_P3) $(PRG_P4)
INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c

doit32_def_SOURCES=doit32_def.c
doit32_def_LDADD=-lmkl32_def $(LIBBENCH) @FLIBS@

doit32_p3_SOURCES=doit32_p3.c
doit32_p3_LDADD=-lmkl32_p3 $(LIBBENCH) @FLIBS@

doit32_p4_SOURCES=doit32_p4.c
doit32_p4_LDADD=-lmkl32_p4 $(LIBBENCH) @FLIBS@

include ../Makefile.common
