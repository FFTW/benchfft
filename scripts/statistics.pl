#! /usr/bin/perl -w

# Extract statistics: mean, median, min, and max mflops for each problem

%probdata = ();

while (<>) {
  my ($nam, $prob, $siz, $mflops, $tim, $setup_tim) = split / /;
  
  next if ($mflops == 0.0 || $mflops eq "");

  $problem = "$prob $siz";
  if (exists($probdata{$problem})) {
      $probdata{$problem} = "$probdata{$problem} $mflops";
  }
  else {
      $probdata{$problem} = $mflops;
  }
}

foreach $problem (keys %probdata) {
    @data = sort { 1000 * ($a - $b) } split(/ /, $probdata{$problem});
#    print "all $problem @data\n";
    $min = $data[0];
    $max = $data[$#data];
    if (($#data + 1) % 2 == 0) {
	$median = 0.5 * ($data[($#data + 1) / 2] + $data[($#data - 1) / 2]);
    }
    else {
	$median = $data[$#data / 2];
    }
    $median = $data[($#data + 1) / 2];
    $sum = 0;
    for ($i = 0; $i <= $#data; $i = $i + 1) {
	$sum = $sum + $data[$i];
    }
    $mean = $sum / ($#data + 1);
    print "min $problem $min\n";
    print "max $problem $max\n";
    print "median $problem $median\n";
    print "mean $problem $mean\n";
}
