#! /usr/bin/perl -w

# insert an empty line when the problem changes

my $oprob="";
while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim) = split / /;
  $prob .= $siz;
  if ($prob ne $oprob) { print "\n"; }
  print;
  $oprob=$prob;
}
