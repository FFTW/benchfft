AUTOMAKE_OPTIONS = gnu

SUBDIRS = libbench benchees scripts

EXTRA_DIST = COPYRIGHT bootstrap.sh acx_f77_cmain_fflags.m4

OUTFILE=

collect:
	$(SHELL) ${top_srcdir}/scripts/collect ${top_builddir}/benchees $(OUTFILE)

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive accuracy-recursive fftinfo-recursive
benchmark: benchmark-recursive
benchmark-am:
accuracy: accuracy-recursive
accuracy-am:
fftinfo: fftinfo-recursive
fftinfo-am:

