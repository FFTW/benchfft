if DOUBLE
if HAVE_CXX
PRG=doit_dif doit_dit doit_split doit_fht doit_4step doit_twodim	\
doit_ndim doitr_fht doitr_wrap doitr_split
endif
endif

INCLUDES = $(INCLBENCH)
PKG=fxt
TARFILE=$(PKG)-2003.03.26.tar.gz
EXTRA_DIST = $(TARFILE) doit.cc doitr.cc doit-doc.h

CPPFLAGS = @CPPFLAGS@ -I$(PKG)/include -I$(PKG)/auxil -I$(PKG)/permute -I$(PKG)/auxbit

LIBFXT = $(PKG)/libfxt.a

doit_dit_SOURCES = doit_dit.cc
doit_dit_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dit.o: doit.cc doit-doc.h $(LIBFXT)

doit_dif_SOURCES = doit_dif.cc
doit_dif_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dif.o: doit.cc doit-doc.h $(LIBFXT)

doit_split_SOURCES = doit_split.cc
doit_split_LDADD = $(LIBFXT) $(LIBBENCH)
doit_split.o: doit.cc doit-doc.h $(LIBFXT)

doit_fht_SOURCES = doit_fht.cc
doit_fht_LDADD = $(LIBFXT) $(LIBBENCH)
doit_fht.o: doit.cc doit-doc.h $(LIBFXT)

doit_4step_SOURCES = doit_4step.cc
doit_4step_LDADD = $(LIBFXT) $(LIBBENCH)
doit_4step.o: doit.cc doit-doc.h $(LIBFXT)

doit_twodim_SOURCES = doit_twodim.cc
doit_twodim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_twodim.o: doit.cc doit-doc.h $(LIBFXT)

doit_ndim_SOURCES = doit_ndim.cc
doit_ndim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_ndim.o: doit.cc doit-doc.h $(LIBFXT)

doitr_fht_SOURCES = doitr_fht.cc
doitr_fht_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_fht.o: doitr.cc doit-doc.h $(LIBFXT)

doitr_wrap_SOURCES = doitr_wrap.cc
doitr_wrap_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_wrap.o: doitr.cc doit-doc.h $(LIBFXT)

doitr_split_SOURCES = doitr_split.cc
doitr_split_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_split.o: doitr.cc doit-doc.h $(LIBFXT)

$(LIBFXT): extract
	cd $(PKG) && (unset CXXFLAGS; $(MAKE))

include ../Makefile.common
