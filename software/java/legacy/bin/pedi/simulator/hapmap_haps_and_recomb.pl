#!/usr/bin/perl

use strict;
use POSIX;

my $THETA = 9 * (10**-9);


my $data_path = "/ubc/cs/home/b/bbkirk/data/CEU-UNRELATED";

my $file_prefix = "hapmap3_r2_b36_fwd.consensus.qc.poly.chr";
my $file_suffix = "_ceu.unr.phased";


if ($#ARGV < 4)
{
    print "\nUSAGE: $0  [d=SNP density] [t=total num of SNPs] [output SNP list] [output hap prefix] [output rec. prefix]\n";
    print "\n";
    print "From among the segregating SNPs, choose every d-th SNP for a total of t SNPs.\n";
    print "\n";
    print "Print the possible founder haplotypes.\n";
    print "Generate the recomb. rates using the distances between loci.\n";
    print "\n";
    print "\n";
    exit(-1);
}

my $d = shift(@ARGV);
my $t = shift(@ARGV);
my $output_snps = shift(@ARGV);
my $output_haps = shift(@ARGV);
my $output_rrates = shift(@ARGV);



for (my $i = 10;  $i <= 10;  $i++)
{

    open(RATES, ">$output_rrates.$i") || die("Cannot write to $output_rrates.$i");
    open(SNPS, ">$output_snps.$i") || die("Cannot write to $output_snps.$i");


    open(FILE, "$data_path/$file_prefix$i$file_suffix") || die("Cannot read $data_path/$file_prefix$i$file_suffix");
    my @haps;
    my $count = 0;
    my $num_snps = 0;
    my $prev_pos = -1;
    my $line = <FILE>;  # column labels
    while ($line = <FILE>)
    {
	chomp($line);
	my @pieces = split(/\s+/, $line);

	

	my $isSegregating = 0;
	my $zero_allele = $pieces[2];
	for ($a = 2;  $a <= $#pieces;  $a++)
	{
	    if ($zero_allele ne $pieces[$a])
	    {
		$isSegregating = 1;
	    }
	}
	#print "$pieces[1] $isSegregating  mod=". ($count % $d) . "\n";

	if ($isSegregating)
	{
	    $count++;
	    if ($count % $d == 0)
	    {
		print "$pieces[0]\n";
		print SNPS "$pieces[0]\n";

		# save the haplotype
		for ($a = 2;  $a <= $#pieces;  $a++)
		{
		    if ($zero_allele eq $pieces[$a])
		    {
			$haps[$a-2] .= "0 ";
		    }
		    else
		    {
			$haps[$a-2] .= "1 ";
		    }
		}


		# compute recomb. rate
		my $this_pos = $pieces[1];
		my $diff = $this_pos - $prev_pos;
		my $pr_recomb = 1 - exp(-$THETA*$diff);

		if ($prev_pos != -1)
		{
		    printf RATES ("%f\n", $pr_recomb);
		}

		$num_snps++;
		$prev_pos = $this_pos;
	    }
	}

	if ($num_snps >= $t)
	{
	    last;
	}
    }
    close(FILE);

    close(SNPS);
    close(RATES);


    open(HAPS, ">$output_haps.$i") || die("Cannot write to $output_haps.$i");
    for (my $c = 0;  $c <= $#haps;  $c++)
    {
	print HAPS $haps[$c] . "\n";
    }
    close(HAPS);

}
