#! /usr/bin/perl -w

# group data by problem

my $oprob="";
while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim) = split / /;
  my $prob0 = $prob . $siz;
  if ($prob0 ne $oprob) { print "\n"; print "$prob $siz:\n"; }
  printf("%30s %12.2f %14.5g\n", $nam, $mflops, $tim);
  $oprob=$prob0;
}
