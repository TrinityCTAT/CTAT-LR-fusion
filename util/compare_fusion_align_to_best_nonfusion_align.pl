#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GMAP_gff3_parser;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 gmap.map.gff3.chims_described gmap.map.gff3.chims_described.fasta.gmap.gff3 EXTEND_LEN\n\n\n";

my $chims_described_file = $ARGV[0] or die $usage;
my $chims_gmap_gff3_file = $ARGV[1] or die $usage;
my $EXTEND = $ARGV[2] or die $usage;


main: {

    my @spans = &GMAP_gff3_parser::parse_GMAP_gff3_alignments($chims_gmap_gff3_file);

    my %trans_to_align_info = &convert_spans_to_align_summary(@spans);
    
    open (my $fh, $chims_described_file) or die "Error, cannot open file $chims_described_file";

    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        my $target_trans = $x[0];
        my $chim_info = $x[3];
        my @c = split(/;/, $chim_info);
        
        my $brkpt_left = $c[2];
        my $brkpt_right = $c[6];

        my $brkpt = int( ($brkpt_left + $brkpt_right)/2 + 0.5);
        
        my $left_align_href = $trans_to_align_info{$target_trans . "____left"};
        my $delta_left = "NA";
        if ($left_align_href) {
            my @delta_left_brkpts = ( abs($left_align_href->{range_lend} - $brkpt_left), abs($left_align_href->{range_rend} - $brkpt_left) );
            @delta_left_brkpts = sort {$a<=>$b} @delta_left_brkpts;
            $delta_left = shift @delta_left_brkpts;
        }
                
        my $right_align_href = $trans_to_align_info{$target_trans . "____right"};
        my $delta_right = "NA";
        
        if ($right_align_href) {
            
            my @delta_right_brkpts = ( abs( $right_align_href->{range_lend} - $EXTEND), abs($right_align_href->{range_rend} - $EXTEND) );
            @delta_right_brkpts = sort {$a<=>$b} @delta_right_brkpts;
            $delta_right = shift @delta_right_brkpts;
        }

        my @deltas;
        if ($delta_left =~ /\d/) {
            push (@deltas, $delta_left);
        }
        if ($delta_right =~ /\d/) {
            push (@deltas, $delta_right);
        }
        my $max_delta = "NA";
        if (@deltas) {
            @deltas = sort {$a<=>$b} @deltas;
            $max_delta = pop @deltas;
        }
                    
        
        print join("\t", @x,
                   #$align_href->{chr},
                   #$align_href->{range_lend},
                   #$align_href->{range_rend},
                   #$align_href->{orient},
                   #$align_href->{per_id},
                   $delta_left,
                   $delta_right,
                   $max_delta,
            ) . "\n";
    }
    close $fh;

}



####
sub convert_spans_to_align_summary {
    my (@spans) = @_;

    my %trans_id_to_summary_info;

    foreach my $span (@spans) {
        
        my $target = $span->{target};

        if (exists $trans_id_to_summary_info{$target}) {
            die "Error, $target already processed ... should only have one alignment per target";
        }
        $trans_id_to_summary_info{$target} = $span;
    }

    return(%trans_id_to_summary_info);
}
