#! /usr/bin/perl -w

# print only power-of-two problems

while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim) = split / /;
  my $rank = ($siz =~ s/x/*/g);
  my $tot = eval($siz);

  if (($tot & ($tot - 1)) != 0) {
    print;
  }
}

