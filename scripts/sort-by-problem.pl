#! /usr/bin/perl -w

# inefficient
sub byproblem {
  my ($anam, $aprob, $asiz, $amflops, $atim) = split / /,$a;
  my ($bnam, $bprob, $bsiz, $bmflops, $btim) = split / /,$b;
  my $arank = ($asiz =~ s/x/*/g);
  my $brank = ($bsiz =~ s/x/*/g);
  my $atot = eval($asiz);
  my $btot = eval($bsiz);

  $aprob cmp $bprob ||    # sort problem kind
  $arank <=> $brank ||    # if tie, sort by rank
  $atot <=> $btot ||      # if tie, sort by total transform size
  $asiz cmp $bsiz ||      # if tie, sort alphabetically on problem description
  $bmflops <=> $amflops;  # if tie, reverse sort by speed
}

while (<>) {
  push(@records, $_);
}

@records = sort byproblem @records;

for (@records) {
  print;
}

