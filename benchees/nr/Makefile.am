if SINGLE
if HAVE_F77
PRG=doit_c doit_f
noinst_LIBRARIES=libfft.a
else
PRG=doit_c
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST = doit.c
LOCAL_COPIES = four1f.f fournf.f realftf.f rlft3f.f

doit_c_SOURCES=doit-c.c my_nrerror.c
nodist_doit_c_SOURCES = four1.c fourn.c realft.c rlft3.c
doit_c_LDADD=$(LIBBENCH)
doit-c.o: doit.c

doit_f_SOURCES=doit-f.c
doit_f_LDADD=libfft.a $(LIBBENCH) @FLIBS@
doit-f.o: doit.c

nodist_libfft_a_SOURCES = four1f.f fournf.f realftf.f rlft3f.f

four1f.f: four1.f
	cp -f $< $@

fournf.f: fourn.f
	cp -f $< $@

realftf.f: realft.f
	cp -f $< $@

rlft3f.f: rlft3.f
	cp -f $< $@

include ../Makefile.common
