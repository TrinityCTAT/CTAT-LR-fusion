#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use List::Util qw(shuffle);
use Data::Dumper;

my $usage = <<__EOUSAGE__;

############################################################################################
#
#  --FI_LR_sam <string>             : sam alignment file for the long reads
#
#  --LR_fusion_report <string>      :  fusion report indicating long reads as evidence.
#
#  --max_alignments_per_fusion <int>  : limit to max of <int> LR read alignments per fusion.
#
#############################################################################################


__EOUSAGE__

    ;



my $FI_LR_sam;
my $LR_fusion_report;
my $max_alignments_per_fusion = 0;

my $help_flag;
&GetOptions ( 'help|h' => \$help_flag,
              'FI_LR_sam=s' => \$FI_LR_sam,
              'LR_fusion_report=s' => \$LR_fusion_report,
              'max_alignments_per_fusion=i' => \$max_alignments_per_fusion,
    );


if ($help_flag) {
    die $usage;
}

unless ($LR_fusion_report && $FI_LR_sam) {
    die $usage;
}


main: {
    
    my %fusion_to_LR_accs;
    open(my $fh, $LR_fusion_report) or die "Error, cannot open file $LR_fusion_report";
    my $reader = new DelimParser::Reader($fh, "\t");
    while (my $row = $reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'} or die "Error, cannot get fusion name from row: " . Dumper($row);
        my $LR_accs = $row->{'LR_accessions'};
        my @LR_fusion_reads = split(/,/, $LR_accs);
        if ($max_alignments_per_fusion && scalar(@LR_fusion_reads) > $max_alignments_per_fusion) {
            @LR_fusion_reads = shuffle(@LR_fusion_reads);
            @LR_fusion_reads = @LR_fusion_reads[0..($max_alignments_per_fusion-1)];
        }
        foreach my $LR_fusion_read (@LR_fusion_reads) {
            $fusion_to_LR_accs{$fusion_name}->{$LR_fusion_read} = 1; 
        }
    }
    close $fh;
    if ($FI_LR_sam =~ /\.bam$/) {
        open($fh, "samtools view -h $FI_LR_sam | ") or die "Error, cannot read $FI_LR_sam via samtools";
    }
    else {
        open($fh, $FI_LR_sam) or die "Error, cannot open file: $FI_LR_sam";
    }
    while(<$fh>) {
        my $line = $_;
        if ($line =~ /^\@/) {
            print $line;
            next;
        }
        my @x = split("\t", $line);
        my $read_acc = $x[0];
        my $fusion_name = $x[2];
        if (exists $fusion_to_LR_accs{$fusion_name}->{$read_acc}) {
            print $line;
        }
    }
        
    exit(0);
    
}

