#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use DelimParser;


my $usage = <<__EOUSAGE__;

##################################################################
#
#  --fusions <string>   :    fusions report
#
#  --reads_fasta <string> :  reads fasta file
#
#  --reads_output <string> : fusion evidence reads output filename
#
##################################################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions;
my $reads_fasta;
my $reads_output;


&GetOptions ( 'help|h' => \$help_flag,
	      'fusions=s' => \$fusions,
	      'reads_fasta=s' => \$reads_fasta,
	      'reads_output=s' => \$reads_output,
	      
    );

if ($help_flag) {
    die $usage;
}

unless ($fusions && $reads_fasta && $reads_output) {
    die $usage;
}

main: {

    my %fusions;
    open(my $fh, $fusions) or die "Error, cannot open file $fusions";
    my $delim_reader = new DelimParser::Reader($fh, "\t");

    my %read_to_fusion;
    while (my $row = $delim_reader->get_row()) {
	my $fusion_name = $row->{"#FusionName"};
	my @LR_accs = split(/,/, $row->{'LR_accessions'});
	foreach my $LR_acc (@LR_accs) {
	    $read_to_fusion{$LR_acc} = $fusion_name;
	}
    }
    close $fh;


    my %accs_to_capture = %read_to_fusion;
    
    open(my $ofh, ">$reads_output") or die "Error, cannot write to $reads_output";
    
    my $fasta_reader = new Fasta_reader($reads_fasta);
    while (my $fasta_entry = $fasta_reader->next()) {
	my $acc = $fasta_entry->get_accession();
	if (my $fusion_name = $read_to_fusion{$acc}) {
	    my $sequence = $fasta_entry->get_sequence();
	    print  $ofh ">$fusion_name|$acc\n$sequence\n";
	    delete $accs_to_capture{$acc};
	}
    }

    if (%accs_to_capture) {
	die "Error, missing capture for the following fusion reads: " . Dumper(\%accs_to_capture);
    }   
	    

    close $ofh;

    print STDERR "-extracted fusion reads into: $reads_output\n";
    
    exit(0);

    
}
