#! /usr/bin/perl -w

use strict;
use warnings;

# Define essential oils
my @essential_oils = (1000001, 1000002, 1000003, 1000004, 1000005, 
                      1000006, 1000007, 1000008);

# Open output file
open(my $fh, '>', "pair1") or die "Cannot open output file: $!";

# Generate all unique pairs
for my $i (0 .. $#essential_oils) {
    for my $j ($i+1 .. $#essential_oils) {
        print $fh "$essential_oils[$i]\t$essential_oils[$j]\t0\n";
    }
}

# Close file
close($fh);

print "Pair file 'pair1' created successfully!\n";
