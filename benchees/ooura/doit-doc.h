BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("package", "ooura")
BENCH_DOC("author", "Takuya Ooura")
BENCH_DOC("version", "2001/11/22")
BENCH_DOC("year", "2001")
BENCH_DOC("email", "ooura@mmm.t.u-tokyo.ac.jp")
BENCH_DOC("url", "http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html")
BENCH_DOC("url-was-valid-on", "Sat Feb 23 09:00:34 EST 2002")
BENCH_DOC("copyright", 
	  "Copyright(C) 1996-2001 Takuya OOURA\n"
	  "You may use, copy, modify this code for any purpose and\n"
	  "without fee. You may distribute this ORIGINAL package.")
#ifdef FORTRAN
  BENCH_DOC("language", "Fortran 77")
#else
  BENCH_DOC("language", "C")
#endif
BENCH_DOC("notes", "Inverse real transform is scaled by .5")
BENCH_DOC("notes", "Sign of forward real transform is +1")
END_BENCH_DOC
