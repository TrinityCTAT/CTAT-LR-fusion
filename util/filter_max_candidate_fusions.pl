#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;


#########################################################################################################
#
# --input_file <str>          input fusions file
#
# --output_file <str>         output fusions file
#
# --max_candidates <int>      maximum number of candidates to explore for fusion contig modeling.
#
###########################################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $input_file;
my $output_file;
my $max_candidates = -1;

my $ALT_MAX_EXON_DELTA = 1000;

&GetOptions ( 'help|h' => \$help_flag,
              'input_file=s' => \$input_file,
              'output_file=s' => \$output_file,
              "max_candidates=i" => \$max_candidates,
    );

if ($help_flag) {
    die $usage;
}


unless (defined($input_file) && defined($output_file) && $max_candidates > 0) {
    die $usage;
}


main: {
    
    open(my $fh, $input_file) or die "Error, cannot open $input_file";
    open(my $ofh, ">$output_file") or die "Error, cannot open $output_file";

    my $counter = -1;
    while(my $line = <$fh>) {
        $counter++; 
        if ($counter <= $max_candidates) {
            # 0 header
            # 1 first fusion
            # ...

            print $ofh $line;
        }
        else {
            last;
        }
    }

    $counter--;
    
    print STDERR "-wrote $counter fusion candidates to $output_file\n";
    close $fh;
    close $ofh;


    exit(0);
}
