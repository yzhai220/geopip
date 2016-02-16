#!/usr/bin/perl -w

use strict;


my $merlin = "merlin-1.1.2/merlin";


if ($#ARGV < 3)
{
    print "USAGE: $0  [num sites] [in ped file] [in rec. rate file] [out score file]\n";
    exit(-34);
}

my $num_sites = shift(@ARGV);
my $in_ped = shift(@ARGV);
my $in_recrates = shift(@ARGV);
my $out_score = shift(@ARGV);

my $dat = "$in_ped.dat";
my $map = "$in_ped.map";


# make the dat file
open(OUT, ">$dat") || die("Cannot write to $dat");
print OUT " A DISEASE\n";
for (my $i=0; $i < $num_sites;  $i++)
{
  print OUT " M SNP-$i\n"
}
print OUT " E END-OF-DATA\n";
close(OUT);



#make_map_files.pl $num_sites $in_recrate $scratch/map.

# make the map file

open(MAP, ">$map") || die("Cannot write to $map");
print MAP "CHROMOSOME   MARKER          POSITION\n";
print MAP "1            SNP-0           0\n";
open(REC, $in_recrates) || die("Cannot read $in_recrates");
my $total = 0;
my $snp = 1;
while (my $line = <REC>)
{
    chomp($line);
    if ($line =~ /[.0-9]+/)
    {
	$total += $line * 100;
	print MAP "1            SNP-$snp           $total\n";
	$snp++;
    }
}
close(REC);
close(MAP);



`ulimit -v 8388608`;
my $output = `$merlin  --bits 1000 -d $dat  -p $in_ped -m $map --npl --markerNames --tabulate`;
`ulimit -v unlimited`;

print $output;

`mv merlin-nonparametric.tbl $out_score`

