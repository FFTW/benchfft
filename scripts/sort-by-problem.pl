#! /usr/bin/perl -w

my @records=();
my @probs=();
my @ranks=();
my @tots=();
my @sizs=();
my @mflops=();
my @nams=();
while (<>) {
  push(@records, $_);
  my ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
  my $rank = ($siz =~ s/x/*/g);
  my $tot = eval($siz);
  push(@nams, $nam);
  push(@probs, $prob);
  push(@ranks, $rank);
  push(@tots, $tot);
  push(@sizs, $siz);

  $mflops=0 if ($mflops eq "FAILED");
  push(@mflopss, $mflops);
}

@records = @records[ sort {
  $probs[$a] cmp $probs[$b] ||    # sort by problem kind
  $ranks[$a] <=> $ranks[$b] ||    # if tie, sort by rank
  $tots[$a] <=> $tots[$b] ||      # if tie, sort by total transform size
  $sizs[$a] cmp $sizs[$b] ||      # if tie, sort alphabetically on problem description
  $mflopss[$b] <=> $mflopss[$a] || # if tie, reverse sort by speed
  $nams[$a] cmp $nams[$b]         # if tie, sort  by name
  } 0..$#records
];

for (@records) {
  print;
}

