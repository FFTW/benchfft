dnl @synopsis AX_F77_CMAIN_FFLAGS([ACTION-IF-SUCCEED], [ACTION-IF-FAIL])
dnl
dnl This macro figures out if extra Fortran compiler flags are required
dnl in order to use the Fortran linker to link programs where the main()
dnl function is defined via C (or other language).  On some systems,
dnl notably the Alpha with Compaq compilers, the Fortran libraries have
dnl their own main() function which must be disabled.
dnl
dnl Runs ACTION-IF-SUCCEED if successful, and ACTION-IF-FAIL if not.
dnl Defines the output variable F77_CMAIN_FFLAGS to any discovered flags.
dnl (If ACTION-IF-FAIL is not specified, defaults to halting with an error.)
dnl
dnl This macro is especially useful in conjunction with automake, since
dnl by default automake uses $F77 to link programs mixing C and Fortran,
dnl leading to a link error on some systems.  In this case, you should
dnl set the FFLAGS for that program to include F77_CMAIN_FFLAGS.
dnl
dnl @version $Id: ax_f77_cmain_fflags.m4,v 1.1 2005-02-12 23:31:11 stevenj Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>

AC_DEFUN([AX_F77_CMAIN_FFLAGS],
[AC_CACHE_CHECK([for f77 flags to use C main function], ax_cv_f77_cmain_fflags,
[ax_cv_f77_cmain_fflags="unknown"
AC_LANG_PUSH(C)
AC_COMPILE_IFELSE([[int main(void) { return 0; }]],
		  [mv conftest.$ac_objext conftest_cmain.$ac_objext],
		  [ax_cv_f77_cmain_fflags=error])
AC_LANG_POP(C)
if test "x$ax_cv_f77_cmain_fflags" != xerror; then
    AC_LANG_PUSH(Fortran 77)
    ax_save_LIBS=$LIBS
    LIBS="conftest_cmain.$ac_objext $LIBS"
    ax_save_FFLAGS=$FFLAGS
    for ax_flag in none -nofor_main; do
	case $ax_flag in
	    none) FFLAGS=$ax_save_FFLAGS ;;
	    *)    FFLAGS="$ax_save_FFLAGS $ax_flag" ;;
	esac
	AC_LINK_IFELSE([
      subroutine foobar()
      return
      end
], [ax_cv_f77_cmain_fflags=$ax_flag; break]);
    done
    FFLAGS=$ax_save_FFLAGS
    LIBS=$ax_save_LIBS
    AC_LANG_POP(Fortran 77)
fi])
    case $ax_cv_f77_cmain_fflags in
	error|unknown) 
	    F77_CMAIN_FFLAGS=""
	    ifelse([$2],,[AC_MSG_ERROR([cannot link C main with Fortran])],[$2])
	    ;;
	*) 
	    if test "x$ax_cv_f77_cmain_fflags" = xnone; then
		F77_CMAIN_FFLAGS=""
	    else
		F77_CMAIN_FFLAGS="$ax_cv_f77_cmain_fflags"
	    fi
	    $1
	    ;;
    esac
    AC_SUBST(F77_CMAIN_FFLAGS)
])
