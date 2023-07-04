#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;

my $usage = <<__EOUSAGE__;

#########################################
#
#  --fusions <string>    :   LR  fusions_file
#
#  --num_LR_total <int>  :   total number of long reads in the sample
#
# --output_file <str>    output filename
#
#########################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions_filename;
my $num_LR_total;
my $output_filename;

&GetOptions ( 'help|h' => \$help_flag,
	      'fusions=s' => \$fusions_filename,
	      'num_LR_total=i' => \$num_LR_total,
	      'output_file=s' => \$output_filename,
    );


if ($help_flag) {
    die $usage;
}

unless ($fusions_filename && $num_LR_total && $output_filename) {
    die $usage;
}


open(my $fh, $fusions_filename) or die "Error, cannot open file: $fusions_filename";
my $delim_reader = new DelimParser::Reader($fh, "\t");
my @column_headers = $delim_reader->get_column_headers();

open(my $ofh, ">$output_filename") or die "Error, cannot write to $output_filename";
my $delim_writer = new DelimParser::Writer($ofh, "\t", [@column_headers, "LR_FFPM"]);

while (my $row = $delim_reader->get_row()) {
    my $LR_FFPM = sprintf("%.3f", $row->{num_LR} / $num_LR_total * 1e6);
    $row->{LR_FFPM} = $LR_FFPM;
    $delim_writer->write_row($row);
}

exit(0);


    
