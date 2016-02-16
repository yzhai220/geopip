#!/usr/bin/perl

use strict;
use POSIX;

##########################
#
# Create a random pedigree with the structure:
#   l lineages
#   g generations
#   w total leaves
#   f fraction of overlap
#   c chaining generation
# where the lineages are chained by marriages at
# generation c with on average fw/l leaves decending
# from a pair of adjacent chained lineages.
#
#########################


my $STDDEV = 3;
my $DEBUG = 0;


if ($#ARGV < 5)
{
    print "$0 [family] [lineages] [generations] [avg total leaves] [fraction typed parents] [output file]\n";
    #print "$0 [family] [lineages] [generations] [avg total leaves] [fraction typed parents] [overlap fraction] [chaining generation] [output file]\n";
    print "\n";
    print "Create a random OUTBRED pedigree with the structure:\n";
    print "   f families\n";
    print "   l lineages\n";
    print "   g generations\n";
    print "   w total leaves\n";
    print "   o=0.2 fraction of overlap\n";
    print "   t fraction of leaf parents that are typed\n";
    #print "   c chaining generation\n";
    print "where the lineages are chained by marriages at\n";
    print "generation c with on average ow/l leaves decending\n";
    print "from a pair of adjacent chained lineages.\n";
    print "Note that all the leaves are typed\n";
    print "\n";
    exit(-1);
}

my $num_families = shift(@ARGV);
my $lineages = shift(@ARGV);
my $generations = shift(@ARGV);
my $avg_total_leaves = shift(@ARGV);
my $typed_parents = shift(@ARGV);
my $overlap = 0.2; #shift(@ARGV);
my $chain_generation = 2; #shift(@ARGV);
my $output = shift(@ARGV);



if ($chain_generation < $generations-1)
{
    print "At the moment, the chaining generation must be ".($generations-1);
    print "\n";
    $chain_generation = $generations-1;
}


print "PARAMS: families: $num_families, lineages: $lineages, generations: $generations, leaves: $avg_total_leaves, \n";
print "        overlap: $overlap, type parents; $typed_parents, chain generation: $chain_generation,\n";
print "        output: $output\n\n";


if ($generations <= 1){
    print "ERROR: need at least 3 generations for chaining\n";
    exit(-2);
}
if ($chain_generation < 2 || $chain_generation >= $generations)
{
    print "ERROR: chain generation should be in range [2,".($generations-1)."]\n";
    exit(-3);
}







######################
#
#  Box-Muller algorithm
#   for generating gaussian random variables 
#
# Global variable
my $SET = 0;
# contains a Gaussian r.v. when $SET == 1
my $Y = -1;
sub bmGaussianPolar
{
    if ($SET == 0)
    {
         my ($x1, $x2, $w, $y1, $y2);
 
	 $w = 2.0;
         while ( $w >= 1.0 )
	 {
                 $x1 = 2.0 * rand(1) - 1.0;
                 $x2 = 2.0 * rand(1) - 1.0;
                 $w = $x1 * $x1 + $x2 * $x2;
         }

         $w = sqrt( (-2.0 * log( $w ) ) / $w );
         $y1 = $x1 * $w;
         $y2 = $x2 * $w;

	 $Y = $y2;

	 return $y1;
     }

    $SET = 0;
    return $Y;

}



open (OUT, ">$output") || die("Could not open $output");
for (my $family = 1;  $family <= $num_families;  $family++)
{
    print "Family $family\n";
    &generate_pedigree($family);
}
close(OUT);





sub generate_pedigree
{
    my $family = shift(@_);

    my $individual_id = 1;
    #print "FAM: $family\n";

    # print the founders and
    # draw the number of leaves for each lineage
    my @leaves;
    my @branching;
    my $total_leaves = 0;
    my $avg_lineage_leaves = $avg_total_leaves/$lineages;
    my $parents = [];
    for (my $l = 0;  $l < $lineages;  $l++)
    {
	print "  Lineage $l\n";
	
	print OUT "$family $individual_id 0 0 1 1\n";
	print "$family $individual_id 0 0 1 1\n" if ($DEBUG);
	$individual_id++;
	print OUT "$family $individual_id 0 0 2 1\n";
	print "$family $individual_id 0 0 2 1\n" if ($DEBUG);
	$individual_id++;

	my $l_parents = [];
	push(@{$l_parents}, $individual_id-2, $individual_id-1);
	$parents->[$l] = $l_parents;


	my $r = $STDDEV * &bmGaussianPolar() + $avg_lineage_leaves;
	if ($r - floor($r) > 0.5) {
	    $r = ceil($r);
	} else {
	    $r = floor($r);
	}
	$total_leaves += $r;
	print "    leaves: $r\n";
	push(@leaves, $r);
	my $b = $r;
	if ($generations > 2){
	    $b = floor(log($r) / log($generations-1));
	}
	print "    branch: $b\n";
	push(@branching, $b);
    }
    print "\n";
    print "total leaves: $total_leaves\n";
    print "\n";





     # for each subsequent generation, EXCEPT the last generation
    for (my $g = 1; $g < $generations-1; $g++)
    {
	my %type;
	if ($typed_parents > 0)
	{
	    # count the number of individuals in this generation
	    my $total_people = 0;
	    if ($g == $chain_generation-1){
		for (my $l = 0; $l < $lineages;  $l++)
		{
		    my $parent_pairs = ($#{$parents->[$l]}+1)/2;
		    for (my $p = 0;  $p < $parent_pairs;  $p++)
		    {
			my $branches = $branching[$l];
			if ($p == $parent_pairs-1){ 
			    if ($l == 0){
				$branches = $branching[$l]-1;
			    } elsif ($l == $lineages-1){
				$branches = $branching[$l]-1;
			    } else {
				$branches = $branching[$l]-2;
			    }
			}
			$total_people += 2*$branches;
			if ($p == $parent_pairs-1 && $l != $lineages-1){
			    $total_people += 2;
			}
		    }
		}
	    }
	    # select the individuals to type
	    my $num_type = floor($total_people * $typed_parents);
	    my $count_selected = 0;
	    while ($count_selected < $num_type)
	    {
		my $notNew = 1;
		while ($notNew)
		{
		    my $indiv = int(rand($total_people));
		    $indiv += $individual_id;
		    if (defined($type{$indiv})){
			$notNew = 1;
		    } else {
			$notNew = 0;
			$type{$indiv} = 2;
		    }
		}
		$count_selected++;
	    }
	}


	for (my $l = 0; $l < $lineages;  $l++)
	{
	    my $new_parents = [];
	    my $parent_pairs = ($#{$parents->[$l]}+1)/2;
	    for (my $p = 0;  $p < $parent_pairs;  $p++)
	    {
		print "gen $g, lin $l, pair $p\n" if ($DEBUG);
		my $father = $parents->[$l][$p*2];
		my $mother = $parents->[$l][$p*2+1];

		my $branches = $branching[$l];
		if ($g == $chain_generation-1 && $p == $parent_pairs-1)
		{ 
		    if ($l == 0){
			$branches = $branching[$l]-1;
		    } elsif ($l == $lineages-1){
			$branches = $branching[$l]-1;
		    } else {
			$branches = $branching[$l]-2;
		    }
		}

		# branches not involved in chaining
		for (my $child=0;  $child < $branches; $child++)
		{
		    my $isTyped = 1;
		    $isTyped = 2 if (defined($type{$individual_id}));
		    print OUT "$family $individual_id $father $mother 1 $isTyped\n";
		    print "$family $individual_id $father $mother 1 $isTyped\n" if ($DEBUG);
		    $individual_id++;
		    $isTyped = 1;
		    $isTyped = 2 if (defined($type{$individual_id}));
		    print OUT "$family $individual_id 0 0 2 $isTyped\n";
		    print "$family $individual_id 0 0 2 $isTyped\n" if ($DEBUG);
		    $individual_id++;
		    push(@{$new_parents}, $individual_id-2, $individual_id-1);
		}

		# all but the last lineage chains on the right-most branch
		if ($g == $chain_generation-1 && $p == $parent_pairs-1 && $l != $lineages-1){
		    my $f2 = $parents->[$l+1][0];
		    my $m2 = $parents->[$l+1][1];

		    my $isTyped = 1;
		    $isTyped = 2 if (defined($type{$individual_id}));
		    print OUT "$family $individual_id $father $mother 1 $isTyped\n";
		    print "$family $individual_id $father $mother 1 $isTyped\n" if ($DEBUG);
		    $individual_id++;

		    my $isTyped = 1;
		    $isTyped = 2 if (defined($type{$individual_id}));
		    $f2 = 0 if (!defined($f2));
		    $m2 = 0 if (!defined($m2));
		    print OUT "$family $individual_id $f2 $m2 2 $isTyped\n";
		    print "$family $individual_id $f2 $m2 2 $isTyped\n" if ($DEBUG);
		    $individual_id++;
		    push(@{$new_parents}, $individual_id-2, $individual_id-1);
		}
	    }
	    $parents->[$l] = $new_parents;
	}
    }




    # LAST generation
    my $g = $generations-1;
    my $total_sum_of_leaves = 0;
    for (my $l = 0; $l < $lineages;  $l++)
    {
	my $parent_pairs = ($#{$parents->[$l]}+1)/2;
	if ($parent_pairs > 0)
	{

	    # calculate the number of leaves for each pair of parents
	    my $leaves_overlapped = 0;
	    if ($l != $lineages-1)
	    {
		$leaves_overlapped = $overlap*$leaves[$l];
		if ($leaves_overlapped < 1){
		    $leaves_overlapped = 1;
		} else {
		    $leaves_overlapped = floor($leaves_overlapped);
		}
	    }

	    my $lineage_leaves = $leaves[$l];
	    my $pair_leaves = ($lineage_leaves - $leaves_overlapped)/$parent_pairs;
	    if ($pair_leaves < 1){
		$pair_leaves = 1;
	    } else {
		$pair_leaves = floor($pair_leaves);
	    }

	    my $equiv_pairs = $parent_pairs-2;
	    $equiv_pairs = 0 if ($equiv_pairs < 0);
	    my $test_sum = ($equiv_pairs)*$pair_leaves + $lineage_leaves - $leaves_overlapped - $equiv_pairs*$pair_leaves + $leaves_overlapped;
	    if ($test_sum != $lineage_leaves)
	    {
		print "ERROR: leaves did not sum properly for lineage $l.\n";
		print "       Summed to $test_sum, but should have been $lineage_leaves.\n";
		print "       ($equiv_pairs)*$pair_leaves + $lineage_leaves - $leaves_overlapped + $leaves_overlapped;\n";
		exit(-345);
	    }

	    # keep track of the number of leaves appearing so far in the lineage
	    my $sum_of_leaves = 0;
	    for (my $p = 0;  $p < $parent_pairs;  $p++)
	    {
		print "gen $g, lin $l, pair $p, s_of_l $sum_of_leaves, pl $pair_leaves\n" if ($DEBUG);
		my $father = $parents->[$l][$p*2];
		my $mother = $parents->[$l][$p*2+1];

		if ($chain_generation == $generations-1)
		{
		    # last pair of parents is the chained pair
		    # the second to the last pair needs to have the 'extra' leaves
		    # (i.e. $lineage_leaves - $leaves_overlapped - $sum_of_leaves)
		    # all the other pairs have $pair_leaves number of leaves.

		    my $num_children = 0;
		    if ($p == $parent_pairs-1)
		    {
			if ($l != $lineages-1){
			    $num_children = $leaves_overlapped;
			} else {
			    $num_children = $lineage_leaves - $leaves_overlapped - $sum_of_leaves;
			}
		    }
		    elsif ($p == $parent_pairs-2)
		    {
			if ($l != $lineages-1)
			{
			    $num_children = $lineage_leaves - $leaves_overlapped - $sum_of_leaves;
			}
			else
			{
			    $num_children = $pair_leaves;
			}
		    }
		    elsif ($p < $parent_pairs-2)
		    {
			$num_children = $pair_leaves;
		    }


		    for (my $child = 0;  $child < $num_children; $child++)
		    {
			print OUT "$family $individual_id $father $mother 2 2\n";
			print "$family $individual_id $father $mother 2 2\n" if ($DEBUG);
			$individual_id++;
		    }
		    $sum_of_leaves += $num_children;
		}
		
	    }

	    if ($sum_of_leaves != $lineage_leaves)
	    {
		print "WARNING: Branching factor was insufficient to guarantee the \n";
		print "         correct number of leaves for lineage $l\n";
		print "         Try increasing the total number of leaves.\n";
	    }
	    $total_sum_of_leaves += $sum_of_leaves;
	}
    }

    print "\nTotal Number of Actual Leaves: $total_sum_of_leaves\n\n";

}
