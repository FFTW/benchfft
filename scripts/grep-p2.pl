#! /usr/bin/perl -w

# print only power-of-two problems

while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
  $siz =~ s/x/*/g;
  my $tot = eval($siz);

  if (($tot & ($tot - 1)) == 0) {
    print;
  }
}

