#! /usr/bin/perl -w

# group data by problem

my $oprob="";
while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
  my $prob0 = $prob . $siz;
  if ($prob0 ne $oprob) { print "\n"; print "$prob $siz:\n"; }
  if ($mflops eq "FAILED") {
    printf("%30s %12s\n", $nam, $mflops);
  } else {
    printf("%30s %12.2f %14.5g %8.1g\n", $nam, $mflops, $tim, $setup_tim);
  }
  $oprob=$prob0;
}
