#! /usr/bin/perl -w

# print only problems of rank argv[1]

@ARGV > 0 || die "usage: $0 RANK";
$r=shift;

while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
  my $rank = 1 + ($siz =~ s/x/*/g);
  my $tot = eval($siz);

  if ($rank == $r) {
    print;
  }
}

