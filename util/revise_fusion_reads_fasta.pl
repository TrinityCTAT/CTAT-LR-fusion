#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Fasta_reader;
use Data::Dumper;

my $usage = "\n\tusage: $0 revised_FI_listing orig_FI_listing_with_reads orig_chim_reads_fasta > revised_chim_reads.fasta\n\n";

my $revised_FI_listing_file = $ARGV[0] or die $usage;
my $orig_FI_listing_with_reads_file = $ARGV[1] or die $usage;
my $orig_chim_reads_fasta_file = $ARGV[2] or die $usage;


my %fusions_want;
{
    open(my $fh, $revised_FI_listing_file) or die "Error, cannot open file: $revised_FI_listing_file";
    my $delim_reader = new DelimParser::Reader($fh, "\t");
    while(my $row = $delim_reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'} or die "Error, no #FusionName entry for row: " . Dumper($row);
        $fusions_want{$fusion_name} = 1;
    }
}

my %reads_want;
{
    open(my $fh, $orig_FI_listing_with_reads_file) or die "Error, cannot open file: $orig_FI_listing_with_reads_file";
    my $delim_reader = new DelimParser::Reader($fh, "\t");
    while(my $row = $delim_reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'} or die "Error, no #FusionName entry for row: " . Dumper($row);
        if ($fusions_want{$fusion_name}) {
            my $reads_str = $row->{'reads'} or die "Error, no reads extracted for row: " . Dumper($row);
            foreach my $readname (split(',', $reads_str)) {
                $reads_want{$readname} = 1;
            }
        }
    }
}

# write new reads fasta
my $fasta_reader = new Fasta_reader($orig_chim_reads_fasta_file);
while(my $fasta_entry = $fasta_reader->next()) {
    my $acc = $fasta_entry->get_accession();
    if ($reads_want{$acc}) {
        print $fasta_entry->get_FASTA_format();
        delete $reads_want{$acc};
    }
}
    
if (%reads_want) {
    die "Error, did not extract reads for : " . Dumper(\%reads_want);
}

exit(0);
