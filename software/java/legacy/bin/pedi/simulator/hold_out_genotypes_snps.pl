#!/usr/bin/perl


##########################
# Check the map from MS ped to SNP ped
# with the SNP pedigree.
# See if the gender, affection status, and 
# type state (typed or untyped) is the same.
##########################


use strict;


if ($#ARGV < 2)
{
    print "USAGE: $0 [hold-out individuals] [file to convert] [output file] \n";
    print "[hold-out individuals] file containing all individuals, hold out those classified affected\n";
    print "[file to convert] pedigree linkage file\n";
    print "[output file] output linkage file\n\n";
    exit(-1);
}


my $map_ped = shift(@ARGV);
my $file_to_convert = shift(@ARGV);
my $out_file = shift(@ARGV);




my $typed_count = 0;


# Read the hold-out list of individuals
my %isTyped;
open(IN, "$map_ped") || die ("Could not open $map_ped");
my $line;
my $line_count = 0;
while ($line = <IN>)
{
    $line_count++;
    chomp($line);
    $line =~ s/^\s+//;

    my @pieces = split(/\s+/, $line);
    if ($#pieces > 5)
    {
	die("Hold out list ($map_ped) has too many columns (should be exactly 6 columns)");
    }
    

    #if ($pieces[6] eq "TYPED")
    if ($pieces[5] eq "2")
    {
	$isTyped{"$pieces[0]:$pieces[1]"} = 1;
	$typed_count++;
    }
    #elsif ($pieces[6] eq "UNTYPED")
    elsif ($pieces[5] eq "1")
    {
	$isTyped{"$pieces[0]:$pieces[1]"} = 0;
    }
    else
    {
	die("ERROR: invalid type status at line $line_count: $pieces[0] $pieces[1] | $pieces[5] $pieces[6] $pieces[7]");
    }

    #print "$pieces[0]:$pieces[1] $pieces[5] ".$isTyped{"$pieces[0]:$pieces[1]"}."\n" ;


}
close(IN);



print "Typed count from hold-out file: $typed_count\n";
$typed_count=0;
my $non_zero=0;


# Read the file to convert
my $number_of_snps = -1;
$line_count = 0;
open(IN, $file_to_convert) || die("Could not open $file_to_convert");
open(OUT, ">$out_file") || die("Could not write to $out_file");
#open(CMP, ">$out_cmp") if ($out_cmp != -1);
while (my $line = <IN>)
{
    $line_count++;
    chomp($line);

    # 0) a pedigree number
    # 1) an individual identification number, or id
    # 2) father's id number
    # 3) mother's id number
    # perhaps (first child, next sibling by father, next sibling by mother)
    # 4) sex (male 1, female 2)  -- 72 male, 60 female
    # 5) affection status (affected 2, unaffected 1, unknown 0) 

    my @pieces = split(/\s+/, $line);


    my $isGenotyped = 0;
    if (defined($isTyped{"$pieces[0]:$pieces[1]"}))
    {
	if ($isTyped{"$pieces[0]:$pieces[1]"} == 1){
	    $isGenotyped = 1;
	    $typed_count++;
	} else {
	    $isGenotyped = 0;
	}
    }
    else
    {
	print "WARNING: could not find individual $pieces[0]:$pieces[1] (line $line_count) in the hold_out list";
	print "         assuming this person should be UNTYPED\n";
	$isGenotyped = 0;
    }

    #print "$pieces[0]:$pieces[1] $isGenotyped\n" ;


    for (my $i = 0;  $i < 6; $i++)
    {
	my $value = shift(@pieces);
	print OUT "$value ";
	#print CMP "$value " if ($out_cmp != -1);
    }
    #print OUT " ";

    my $locus_num = 0;
    while ($#pieces > -1)
    {
	my $value1 = shift(@pieces);
	my $value2 = shift(@pieces);

	if ($isGenotyped){
	    print OUT "  $value1 $value2";
	} else {
	    print OUT "  0 0";
	}

	$locus_num++;
    }

    print OUT "\n";
}
close(IN);
close(OUT);

print "Typed count in output ped: $typed_count\n";

