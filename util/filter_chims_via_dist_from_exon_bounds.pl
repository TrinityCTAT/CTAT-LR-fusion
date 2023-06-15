#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;

my $usage = "usage: $0 gmap.gff3.summary minJ minJS min_novel_J\n\n";

my $gmap_summary_file = $ARGV[0] or die $usage;
my $minJ = $ARGV[1];
my $minJS = $ARGV[2];
my $min_novel_J = $ARGV[3];


my $MIN_ANCHOR_ENTROPY = 1.0;

unless (defined ($minJ) && defined($minJS) && defined($min_novel_J)) {
    die $usage;
}

my $filt_file = "$gmap_summary_file.filt";
open (my $ofh, ">$filt_file") or die "Error, cannot write to $filt_file";

open (my $fh, $gmap_summary_file) or die $!;
my $delim_parser = new DelimParser::Reader($fh, "\t");

my $delim_writer = new DelimParser::Writer(*STDOUT, "\t", [$delim_parser->get_column_headers()]);
my $delim_writer_filt = new DelimParser::Writer($ofh, "\t", [$delim_parser->get_column_headers(), "filtered"]);

while (my $row = $delim_parser->get_row()) {
    
    my $J = $row->{JunctionReadCount};
    my $S = $row->{SpanningFragCount};
    
    my $left_brkpt_entropy = $row->{left_brkpt_entropy};
    my $right_brkpt_entropy = $row->{right_brkpt_entropy};
    
    my $geneA = $row->{LeftGene};
    my $geneB = $row->{RightGene};
    
    my $splice_type = $row->{SpliceType};

    if ($geneA eq $geneB) {
        $row->{filtered} = "NO SELFIES";
        $delim_writer_filt->write_row($row);
        next; 
    }
    
    my $record_pass = 0;
    
    if (defined($J) && defined($S) && $J =~ /^\d+$/ && $S =~ /^\d+$/) {
        
        ## Examining Illumina PE read support
        
        my $sumJS = $J + $S;
        if ($splice_type eq "ONLY_REF_SPLICE"  && $J >= $minJ && $sumJS >= $minJS) {
            # reference junction
            $record_pass = 1;
        }
        elsif ( $splice_type ne "ONLY_REF_SPLICE" && $J >= $min_novel_J  && $sumJS >= $minJS) {
            $record_pass = 1;
        }
        
        ## but...  if only junction reads, the anchors must meet the min entropy requirement.
        if ($J > 0 && $S == 0 && ($left_brkpt_entropy < $MIN_ANCHOR_ENTROPY || $right_brkpt_entropy < $MIN_ANCHOR_ENTROPY) ) {
            $row->{filtered} = "Fails to meet min entropy requirement at junction anchor region";
            $delim_writer_filt->write_row($row);
            next;
        }
            
    }
    
        
    if ($record_pass) {
        $delim_writer->write_row($row);
    }
    else {
        $row->{filtered} = "Fails to meet minJ($J) >= $minJ or minJS_sum($J + $S) >= $minJS requirement or not ref splice and J:$J < min_novel_J: $min_novel_J";
        $delim_writer_filt->write_row($row);
    }
    
}


exit(0);
