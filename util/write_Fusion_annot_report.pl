#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "usage: $0 gmap.map.gff3.chims_described.w_read_support [LONG_READS_ONLY_FLAG]\n\n";

my $chim_file = $ARGV[0] or die $usage;
my $long_reads_only_flag = $ARGV[1] || 0;


=inputfmt

0       Locus_38_Transcript_1/1_Confidence_1.000_Length_212
1       2
2       [chr17:(95-212)39451038-39451155 (+) 100.00%];[chr17:(1-94)55141310-55159832 (-) 100.00%]
3       MED1;0;95;chr17:39451038;STXBP4;0;94;chr17:55141310;MED1--STXBP4
4       1
5       0
6       CATAAGGTAAAC
7       1.79
8       TCGGTTTCCCCC
9       1.46
10      .
11      .

=cut

main: {

    
    my @header_cols = ("#FusionName", "JunctionReadCount", "SpanningFragCount", "trans_acc", "trans_brkpt", 
                       "LeftGene", "LeftBreakpoint", 
                       "RightGene", "RightBreakpoint", 
                       "SpliceType", 
                       "left_brkpt_anchorseq", "left_brkpt_entropy", "right_brkpt_anchorseq", "right_brkpt_entropy");
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", [@header_cols]);
        
    
    my %fusion_token_to_rows;

    open (my $fh, $chim_file) or die "Error, cannot open file $chim_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $trinity_acc = $x[0];
        my $fusion_info = $x[3];
        my $J = $x[4];
        my $S = $x[5];
    
        my $left_brkpt_anchorseq = $x[6];
        my $left_brkpt_entropy = $x[7];
        my $right_brkpt_anchorseq = $x[8];
        my $right_brkpt_entropy = $x[9];
        
        if (! defined($J)) { $J = "."; }
        if (! defined($S)) { $S = "."; }

        my ($geneA, $deltaA, $trans_brkptA, $breakpointA, 
            $geneB, $deltaB, $trans_brkptB, $breakpointB, $fusion_name) = split(/;/, $fusion_info);
        
        
        my ($chrA, $coordA) = split(/:/, $breakpointA);
        my ($chrB, $coordB) = split(/:/, $breakpointB);
        
        my $trans_brkpt = join("-", $trans_brkptA, $trans_brkptB);
        
        my ($junction_type) = ($deltaA != 0 || $deltaB != 0) ? "INCL_NON_REF_SPLICE" : "ONLY_REF_SPLICE";
    

        my $row = { "#FusionName" => $fusion_name,
                    'JunctionReadCount' => $J,
                    'SpanningFragCount' => $S,
                    'trans_acc' => $trinity_acc,
                    'trans_brkpt' => $trans_brkpt,
                    
                    'LeftGene' => $geneA,
                    'LeftBreakpoint' => "$chrA:$coordA",
                    
                    'RightGene' => $geneB,
                    'RightBreakpoint' => "$chrB:$coordB",

                    'SpliceType' => $junction_type,

                    'left_brkpt_anchorseq' => $left_brkpt_anchorseq,
                    'left_brkpt_entropy' => $left_brkpt_entropy,
                    
                    'right_brkpt_anchorseq' => $right_brkpt_anchorseq,
                    'right_brkpt_entropy' => $right_brkpt_entropy,
                    
        };
        
        if ($long_reads_only_flag) {
            my $fusion_token = join("^", $fusion_name, $chrA, $coordA, $chrB, $coordB);
            push (@{$fusion_token_to_rows{$fusion_token}}, $row);
        }
        else {
            $tab_writer->write_row($row);
        }
    }

    if ($long_reads_only_flag) {
        # aggregate long reads info.
        
        my @agg_rows;
        foreach my $fusion_token (keys %fusion_token_to_rows) {
            my @rows = @{$fusion_token_to_rows{$fusion_token}};
            my $repr_row = shift @rows;
            for my $row (@rows) {
                $repr_row->{JunctionReadCount} += $row->{JunctionReadCount};
                $repr_row->{trans_acc} .= ";" . $row->{trans_acc};
            }
            push (@agg_rows, $repr_row);
        }
            
        @agg_rows = reverse sort {$a->{JunctionReadCount}<=>$b->{JunctionReadCount}} @agg_rows;
        
        foreach my $row (@agg_rows) {
            $tab_writer->write_row($row);
        }
    }
    
    close $fh;

    exit(0);
    
}
