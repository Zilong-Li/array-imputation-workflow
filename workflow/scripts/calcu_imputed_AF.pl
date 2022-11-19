#!/usr/bin/env perl -w

use strict;

my ($genfile, $infofile, $chr) = @ARGV;
$chr =~ s/^chr//;

my %h;
open(I, "gzip -dc $infofile |") or die $!;
while(<I>){
  chomp;
  next if $.==1;
  my @t = split;
  my $id = "$chr\t$t[2]\t$t[3]\t$t[4]";
  $h{$id} = $t[6];

}
close I;

print "CHR\tPOS\tREF\tALT\tImputedAF\tINFO\n";
open(I, "gzip -dc $genfile |") or die $!;
while(<I>){
  chomp;
  my @t = split;
  my $id = "$chr\t$t[2]\t$t[3]\t$t[4]";
  @t = @t[5..$#t];
  my $n = scalar(@t) / 3;
  my $dos = 0;
  for (my $i=0; $i<$n; $i++) {
    $dos += $t[$i * 3 + 1] + $t[$i * 3 + 2] * 2;
  }
  $dos /= 2 * $n;
  if (exists $h{$id}) {
    print "$id\t$dos\t$h{$id}\n";
  } else {
    die "$id doesn't exist in $infofile\n";
  }
}
close I;

