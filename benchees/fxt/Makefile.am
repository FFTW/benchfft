if DOUBLE
if HAVE_CXX
PRG = doit_fht doit_cfht \
      doit_dif2 doit_depth_first_dif2 doit_dit2 doit_depth_first_dit2 \
      doit_dif4 doit_cdif4 doit_dit4 doit_cdit4 \
      doit_split doit_csplit \
      doit_twodim doit_ndim \
      doit_gray \
      doit_recursive_dit2 doit_recursive_dif2 \
      doit_matrixfft doit_cmatrixfft \
      doitr_fht doitr_wrap doitr_split
endif
endif

INCLUDES = $(INCLBENCH)
PKG=fxt
TARFILE=$(PKG)-2008.08.31.tar.gz
EXTRA_DIST = $(TARFILE) doit.cc doitr.cc doit-doc.h

AM_CPPFLAGS = -I$(PKG) -I$(PKG)/src

LIBFXT = $(PKG)/libfxt.a

doit_fht_SOURCES = doit_fht.cc
doit_fht_LDADD = $(LIBFXT) $(LIBBENCH)
doit_fht.o: doit.cc doit-doc.h $(LIBFXT)

doit_cfht_SOURCES = doit_cfht.cc
doit_cfht_LDADD = $(LIBFXT) $(LIBBENCH)
doit_cfht.o: doit.cc doit-doc.h $(LIBFXT)

doit_dif2_SOURCES = doit_dif2.cc
doit_dif2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dif2.o: doit.cc doit-doc.h $(LIBFXT)

doit_depth_first_dif2_SOURCES = doit_depth_first_dif2.cc
doit_depth_first_dif2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_depth_first_dif2.o: doit.cc doit-doc.h $(LIBFXT)

doit_dit2_SOURCES = doit_dit2.cc
doit_dit2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dit2.o: doit.cc doit-doc.h $(LIBFXT)

doit_depth_first_dit2_SOURCES = doit_depth_first_dit2.cc
doit_depth_first_dit2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_depth_first_dit2.o: doit.cc doit-doc.h $(LIBFXT)

doit_dif4_SOURCES = doit_dif4.cc
doit_dif4_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dif4.o: doit.cc doit-doc.h $(LIBFXT)

doit_cdif4_SOURCES = doit_cdif4.cc
doit_cdif4_LDADD = $(LIBFXT) $(LIBBENCH)
doit_cdif4.o: doit.cc doit-doc.h $(LIBFXT)

doit_dit4_SOURCES = doit_dit4.cc
doit_dit4_LDADD = $(LIBFXT) $(LIBBENCH)
doit_dit4.o: doit.cc doit-doc.h $(LIBFXT)

doit_cdit4_SOURCES = doit_cdit4.cc
doit_cdit4_LDADD = $(LIBFXT) $(LIBBENCH)
doit_cdit4.o: doit.cc doit-doc.h $(LIBFXT)

doit_split_SOURCES = doit_split.cc
doit_split_LDADD = $(LIBFXT) $(LIBBENCH)
doit_split.o: doit.cc doit-doc.h $(LIBFXT)

doit_csplit_SOURCES = doit_csplit.cc
doit_csplit_LDADD = $(LIBFXT) $(LIBBENCH)
doit_csplit.o: doit.cc doit-doc.h $(LIBFXT)

doit_twodim_SOURCES = doit_twodim.cc
doit_twodim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_twodim.o: doit.cc doit-doc.h $(LIBFXT)

doit_ndim_SOURCES = doit_ndim.cc
doit_ndim_LDADD = $(LIBFXT) $(LIBBENCH)
doit_ndim.o: doit.cc doit-doc.h $(LIBFXT)

doit_gray_SOURCES = doit_gray.cc
doit_gray_LDADD = $(LIBFXT) $(LIBBENCH)
doit_gray.o: doit.cc doit-doc.h $(LIBFXT)

doit_recursive_dif2_SOURCES = doit_recursive_dif2.cc
doit_recursive_dif2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_recursive_dif2.o: doit.cc doit-doc.h $(LIBFXT)

doit_recursive_dit2_SOURCES = doit_recursive_dit2.cc
doit_recursive_dit2_LDADD = $(LIBFXT) $(LIBBENCH)
doit_recursive_dit2.o: doit.cc doit-doc.h $(LIBFXT)

doit_matrixfft_SOURCES = doit_matrixfft.cc
doit_matrixfft_LDADD = $(LIBFXT) $(LIBBENCH)
doit_matrixfft.o: doit.cc doit-doc.h $(LIBFXT)

doit_cmatrixfft_SOURCES = doit_cmatrixfft.cc
doit_cmatrixfft_LDADD = $(LIBFXT) $(LIBBENCH)
doit_cmatrixfft.o: doit.cc doit-doc.h $(LIBFXT)

doitr_fht_SOURCES = doitr_fht.cc
doitr_fht_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_fht.o: doitr.cc doit-doc.h $(LIBFXT)

doitr_wrap_SOURCES = doitr_wrap.cc
doitr_wrap_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_wrap.o: doitr.cc doit-doc.h $(LIBFXT)

doitr_split_SOURCES = doitr_split.cc
doitr_split_LDADD = $(LIBFXT) $(LIBBENCH)
doitr_split.o: doitr.cc doit-doc.h $(LIBFXT)

if CXX_GNU
$(LIBFXT): extract
	cd $(PKG) && (unset CXXFLAGS; $(MAKE) CPUTYPE=$(GCC_ARCH))
else
$(LIBFXT): extract
	cd $(PKG) && $(MAKE)
endif

include ../Makefile.common
