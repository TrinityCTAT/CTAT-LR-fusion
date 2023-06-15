#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 split_align_stats_file MAX_FUZZY_OVERLAP\n\n";

my $split_align_stats_file = $ARGV[0] or die $usage;
my $MAX_FUZZY_OVERLAP = $ARGV[1] or die $usage;


main: {

    open (my $fh, $split_align_stats_file) or die "Error, cannot open file $split_align_stats_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $max_delta = pop @x;
        
        if ($max_delta <= $MAX_FUZZY_OVERLAP) {
            print join("\t", @x[0..3]) . "\n";
        }
    }
    close $fh;

    exit(0);
}


        
