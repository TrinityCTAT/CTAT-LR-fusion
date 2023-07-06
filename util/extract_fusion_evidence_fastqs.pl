#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fastq_reader;
use Process_cmd;
use DelimParser;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);  


my $usage = <<__EOUSAGE__;

######################################################################################
#
#  --fusions <string>           fusion predictions 'final' output file (not abridged)
#
#  Reads:
#
#  --fastq <string>             reads in fastq file #
#
#
######################################################################################


__EOUSAGE__

    ;


my $help_flag;

my $fusion_results_file;
my $fastq;


&GetOptions( 'help|h' => \$help_flag,
             
             'fusions=s' => \$fusion_results_file,
             
             'fastq=s' => \$fastq,

             
    );

if ($help_flag) {
    die $usage;
}

unless ($fusion_results_file && $fastq) {
    die $usage;
}


main: {

    ## get the core fragment names:
    my %LR_read_name_to_fusion_name;
    
    open (my $fh, $fusion_results_file) or die "Error, cannot open file $fusion_results_file";
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'} or die "Error, cannot get fusion name from " . Dumper($row);
        
        my $fusion_name_simple = $fusion_name;

	my $left_brkpt = $tab_reader->get_row_val($row, "LeftBreakpoint");
	my $right_brkpt = $tab_reader->get_row_val($row, "RightBreakpoint");
	my $brkpt_info = "${left_brkpt}^${right_brkpt}";
	$fusion_name .= "^$brkpt_info";
            
        my $junction_reads_list_txt = $tab_reader->get_row_val($row, "LR_accessions");
        
        foreach my $junction_read (split(/,/, $junction_reads_list_txt)) {
            
            $LR_read_name_to_fusion_name{$junction_read} = $fusion_name;
        }
    }
        
    &write_fastq_files($fastq, \%LR_read_name_to_fusion_name);
    
    print STDERR "\nDone.\n\n";
    
    exit(0);
    
}


####
sub write_fastq_files {
    my ($fastq_file, $LR_read_name_to_fusion_name_href) = @_;

    
    my %reads_to_capture = %{$LR_read_name_to_fusion_name_href};
    
    print STDERR "-searching fq file: $fastq_file\n";

    my $fastq_reader = new Fastq_reader($fastq_file);
        
    while (my $fq_record = $fastq_reader->next()) {
	
        my $read_name = $fq_record->get_core_read_name();
	
        if (my $fusion_name = $LR_read_name_to_fusion_name_href->{$read_name}) { 
	
	    my $record_text = $fq_record->get_fastq_record();
	    chomp $record_text;
	                
            my @lines = split(/\n/, $record_text);
            
	    my ($_1, $_2, $_3, $_4) = @lines;
                
	    $_3 = "+$fusion_name"; # encode the fusion name in the 3rd line, which is otherwise useless anyway
	    
	    print join("\n", ($_1, $_2, $_3, $_4)) . "\n";

	    delete $reads_to_capture{$read_name};
        }
    }
    
        
    if (%reads_to_capture) {
        confess "Error, failed to capture fusion evidence reads: " . Dumper(\%reads_to_capture);
    }
    
    return;
}


