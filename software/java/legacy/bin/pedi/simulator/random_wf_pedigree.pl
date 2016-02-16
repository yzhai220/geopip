#!/usr/bin/perl

use strict;
use POSIX;


my $DEBUG = 0;


if ($#ARGV < 5)
{
    print "$0 [n pop size] [g generations] [t fraction typed parents] [s random seed] [m is monogamous] [output file]\n";
    print "\n";
    print "Create a random WF family with the structure:\n";
    print "   n half the population size\n";
    print "   g generations\n";
    print "   t fraction of leaf parents that are typed\n";
    print "   s random seed\n";
    print "   m (bool) is monogamous\n";
    print "where the leafs are all typed";
    print "\n";
    exit(-1);
}

my $pop_size = shift(@ARGV);
my $generations = shift(@ARGV);
my $typed_parents = shift(@ARGV);
my $random_seed = shift(@ARGV);
my $is_monogamous = shift(@ARGV);
my $output_ped = shift(@ARGV);

srand($random_seed);

open(OUT, ">$output_ped") || die("Could not write to $output_ped");



# separate genders
# an array of array of array of array
# first index is gender
# second index is generation
# third index is individual
# fourth index is the parent, male or female
my $pedigree = [];
$pedigree->[0] = [];  # males
$pedigree->[1] = [];  # females

# generation 0 is leafs
for (my $g = 0;  $g < $generations;  $g++)
{
    # create array for generation
    for (my $s = 0;  $s < 2; $s++)
    {
	$pedigree->[$s][$g] = [];
    }

    # create marriages if $is_monogamous
    my $marriages = [];
    if ($is_monogamous)
    {
	my @used; # boolean to track which second parents are in marriages
	for (my $i = 0; $i < $pop_size;  $i++)
	{
	    $used[$i] = 0;
	}


	for (my $i = 0; $i < $pop_size;  $i++)
	{
	    $marriages->[$i] = [];

	    # first parent
	    $marriages->[$i][0] = 0;
	    if ($g < $generations -1)
	    {
		$marriages->[$i][0] = $i;
	    }
	    
	    # second parent
	    $marriages->[$i][1] = 0;
	    if ($g < $generations -1)
	    {
		my $parent = int(rand($pop_size));
		while ($used[$parent])
		{
		    $parent = int(rand($pop_size));
		}
		$marriages->[$i][1] = $parent;
		$used[$parent] = 1;
	    }

	}
    } # end if $is_monogamous

    for (my $s = 0;  $s < 2; $s++)
    {
	for (my $i = 0; $i < $pop_size;  $i++)
	{
	    $pedigree->[$s][$g][$i] = [];

	    if (!$is_monogamous)
	    {
		for (my $p = 0;  $p < 2;  $p++)
		{
		    if ($g < $generations -1)
		    {
			my $parent = int(rand($pop_size));
			$pedigree->[$s][$g][$i][$p] = $parent;
		    }
		    else
		    {
			$pedigree->[$s][$g][$i][$p] = 0;
		    }
		}
	    }
	    else
	    {
		my $m = int(rand($pop_size));
		for (my $p = 0;  $p < 2;  $p++)
		{
		    $pedigree->[$s][$g][$i][$p] = $marriages->[$m][$p];
		}
	    }
	}
    }
}


# Do a dfs traversal from the leaves up

# hash of names
# maps gender:generation:individual to name
my %names;

# entries are gender:generation:individual
my @queue; 

# put leaves in queue
my $indiv_counter = 1;
for (my $s = 0;  $s < 2; $s++)
{
    for (my $i = 0; $i < $pop_size;  $i++)
    {
	my $name = $indiv_counter;
	$indiv_counter++;
	my $sname = "$s:0:$i";
	$names{$sname} = $name;

	my $f;
	my $sf = "0:1:$pedigree->[$s][1][$i][0]";
	if (defined($names{$sf})){
	    $f = $names{$sf};
	} else {
	    $f = $indiv_counter;
	    $names{$sf} = $f;
	    $indiv_counter++;

	    # push the parent
	    push(@queue, "0:1:".($pedigree->[$s][1][$i][0]).":$f");
	}

	my $m;
	my $sm = "1:1:$pedigree->[$s][1][$i][1]";
	if (defined($names{$sm})){
	    $m = $names{$sm};
	} else {
	    $m = $indiv_counter;
	    $names{$sm} = $m;
	    $indiv_counter++;

	    # push the parent
	    push(@queue, "1:1:".($pedigree->[$s][1][$i][1]).":$m");
	}

	print OUT "1 $name $f $m ".($s+1)." 2\n";
    }
}

# dfs
while ($#queue != -1)
{
    my $string = shift(@queue);
    my ($s,$g,$i,$name) = split(/:/, $string);

    if ($pedigree->[$s][$g][$i][0] != 0)
    {

	my $f;
	my $sf = "0:".($g+1).":$pedigree->[$s][$g][$i][0]";
	if (defined($names{$sf})){
	    $f = $names{$sf};
	} else {
	    $f = $indiv_counter;
	    $names{$sf} = $f;
	    $indiv_counter++;

	    # push the parent
	    push(@queue, "0:".($g+1).":".($pedigree->[$s][$g][$i][0]).":$f");
	}

	my $m;
	my $sm = "1:".($g+1).":$pedigree->[$s][$g][$i][1]";
	if (defined($names{$sm})){
	    $m = $names{$sm};
	} else {
	    $m = $indiv_counter;
	    $names{$sm} = $m;
	    $indiv_counter++;

	    # push the parent
	    push(@queue, "1:".($g+1).":".($pedigree->[$s][$g][$i][1]).":$m");
	}



	my $u  = rand(1);
	
	if ($u < $typed_parents) {
	    print OUT "1 $name $f $m ".($s+1)." 2\n";
	} else {
	    print OUT "1 $name $f $m ".($s+1)." 1\n";
	}
	


    }
    else
    {
	my $u  = rand(1);
	
	if ($u < $typed_parents) {
	    print OUT "1 $name 0 0 ".($s+1)." 2\n";
	} else {
	    print OUT "1 $name 0 0 ".($s+1)." 1\n";
	}
    }
}


close(OUT);
