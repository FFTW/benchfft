#! /usr/bin/perl -w

while (<>) {
  push(@records, $_);
  my ($nam, $prob, $siz, $mflops, $tim) = split / /;
  my $rank = ($siz =~ s/x/*/g);
  my $tot = eval($siz);
  push(@probs, $prob);
  push(@ranks, $rank);
  push(@tots, $tot);
  push(@sizs, $siz);
  push(@mflopss, $mflops);
}

@records = @records[ sort {
  $probs[$a] cmp $probs[$b] ||    # sort by problem kind
  $ranks[$a] <=> $ranks[$b] ||    # if tie, sort by rank
  $tots[$a] <=> $tots[$b] ||      # if tie, sort by total transform size
  $sizs[$a] cmp $sizs[$b] ||      # if tie, sort alphabetically on problem description
  $mflopss[$b] <=> $mflopss[$a];  # if tie, reverse sort by speed
  } 0..$#records
];

for (@records) {
  print;
}

