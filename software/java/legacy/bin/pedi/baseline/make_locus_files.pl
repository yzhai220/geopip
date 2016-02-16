#!/usr/bin/perl


####################
# Make locus files for
# Superlink.
####################



use strict;
use POSIX;


if ($#ARGV < 3)
{
    print "USAGE: $0 [block size] [rec rates file] [hold out file] [out prefix]\n";
    print "[block size] number of SNPs in a block\n";
    print "[rec rates file] input recombination rates\n";
    print "[hold out file] input SNP hold out list\n";
    print "[out prefix] prefix for output locus file names\n";
    exit(-1);
}



my $block_size = shift(@ARGV);
my $rec_rates_file = shift(@ARGV);
my $hold_out_file = shift(@ARGV);
my $out_prefix = shift(@ARGV);



# read the rec_rates file
my $block_snp_count = 1;
my $block_count = 0;
my @rec_rates;
open (HOLD, $hold_out_file) || die("Could not read $hold_out_file");
my $hline = <HOLD>;
my $line_count = 1;
open (REC, $rec_rates_file) || die("Could not read $rec_rates_file");
while (my $rline = <REC>)
{
    chomp($rline);
    if ($rline =~ /^([0-9.]+)$/){
	my $rate = $1;
	print "$rate ($block_snp_count, $block_count)\n";

	my ($blah1, $blah2, $unobserved) = split(/\s+/, $hline);
	if ($unobserved eq "UNOBSERVED")
	    {
		if ($#rec_rates >= 0){
		    $rec_rates[$#rec_rates] += $rate;
		}
		$hline = <HOLD>;
		next;
	    }

	if ($block_snp_count <= $block_size-1)
	{
	    push (@rec_rates, $rate);
	    $block_snp_count++;
	}
	else # if ($block_snp_count == $block_size-1)
	{
	    # this recomb. rate is lost, because it is at
	    # a block boundary

	    # output the rec rates for this block
	    print "&print_locus_file($out_prefix$block_count, $block_snp_count, @rec_rates);\n";
	    &print_locus_file("$out_prefix$block_count", $block_snp_count, @rec_rates);

	    splice(@rec_rates, 0, $#rec_rates+1);
	    $block_snp_count = 1;
	    $block_count++;
	}

	$hline = <HOLD>;
    }


}
close (REC);
close (HOLD);

print "&print_locus_file($out_prefix$block_count, $block_snp_count, @rec_rates);\n";
&print_locus_file("$out_prefix$block_count", $block_snp_count, @rec_rates);





sub print_locus_file
{
    my $file = shift(@_);
    my $snp_count = shift(@_);
    my @rec_rates = @_;

    open (OUT, ">$file") || die("Could not write to $file");
    print OUT "".($snp_count+1)." 0 0 5 0\n";
    print OUT "0 0.00 0.00 0\n";
    for (my $c = 1; $c <= $snp_count+1; $c++){
	print OUT "$c ";
    }
    print OUT "\n";
    # output affection locus
    print OUT "1 2 \# AFFECTATION\n";
    print OUT " 0.990000 0.010000\n";
    print OUT "1\n";
    print OUT " 0.0010 0.9000 0.9000\n";

    # output all other loci
    for (my $c = 1; $c <= $snp_count; $c++){
	print OUT "3 2\n 0.5 0.5\n";
    }
    print OUT "0 0\n";

    # output recomb. rates
    print OUT "0 ";
    foreach my $r (@rec_rates){
	printf(OUT "%f ", $r);
	#0.008787 0.000030 
    }
    print OUT "Haldane\n";
    print OUT "1 0.5 0.5\n";

    close(OUT);
}
