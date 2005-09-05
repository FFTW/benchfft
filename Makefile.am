AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS=-I ./m4

SUBDIRS = libbench benchees scripts

EXTRA_DIST = COPYRIGHT bootstrap.sh m4

OUTFILE=

collect:
	$(SHELL) ${top_srcdir}/scripts/collect ${top_builddir}/benchees $(OUTFILE)

# handle recursive ``make benchmark'' etc.
BENCHFFT_RECURSIVE_TARGETS = benchmark-recursive accuracy-recursive	\
fftinfo-recursive

benchmark: benchmark-recursive
benchmark-am:
accuracy: accuracy-recursive
accuracy-am:
fftinfo: fftinfo-recursive
fftinfo-am:

$(BENCHFFT_RECURSIVE_TARGETS):
	@failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	  || eval $$failcom; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"
