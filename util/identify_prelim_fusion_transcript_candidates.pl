#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Fastq_reader;
use File::Basename;
use Process_cmd;
use Pipeliner;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use List::Util qw(min max);

my $usage = <<__EOUSAGE__;


#########################################################################################################
#
# --chims_described <string>   chims.described file.
#
# --max_exon_delta <int>       maximum dist from ref exon boundary
#
# --output_prefix <string>    prefix name for output files (prefix).transcripts.fa and (prefix).FI_listing
#
# --min_num_LR <int>          minimum number of long reads required as evidence.
# 
# --min_FFPM <float>          min fusion expression for candidates to pursue
#
# --num_total_reads <int>     number of total reads (used for FFPM calculation)
#
# --max_foldback_frac <float> maximum fraction of reads that can be fold-backs (default: 0.5, set to 1.0 to disable filtering)
#
###########################################################################################################


__EOUSAGE__

    ;


my $chims_described;
my $MAX_EXON_DELTA;
my $help_flag;
my $output_prefix;
my $min_FFPM;
my $num_total_reads;
my $min_num_LR = 0;
my $max_foldback_frac = 0.5;

my $ALT_MAX_EXON_DELTA = 1000;

&GetOptions ( 'help|h' => \$help_flag,
              'chims_described=s' => \$chims_described,
              'max_exon_delta=i' => \$MAX_EXON_DELTA,
              'output_prefix=s' => \$output_prefix,
              'min_FFPM=f' => \$min_FFPM,
              'num_total_reads=i' => \$num_total_reads,
              'min_num_LR=i' => \$min_num_LR,
              'max_foldback_frac=f' => \$max_foldback_frac,
    );

if ($help_flag) {
    die $usage;
}

unless ($chims_described) {
    print STDERR "\n\nERROR - must specify --chims_desacribed <string>\n";
    die $usage;
}

unless($output_prefix) {
    print STDERR "\n\nERROR - must specify --output_prefix <string>\n";
    die $usage;
}


unless (defined($MAX_EXON_DELTA)) {
    print STDERR " - must specify --max_exon_delta <int>\n";
    die $usage;
}

unless (defined($num_total_reads) && $num_total_reads > 0) {
    print STDERR " - must specify --num_total_reads <int>\n";
    die $usage;
}

unless (defined ($min_FFPM) ) {
    print STDERR "\n\nError, must set --min_FFPM ";
    die $usage;
}



$chims_described = &ensure_full_path($chims_described);


unless (-s $chims_described) {
    confess "Error, cannot locate file $chims_described";
}


main: {
    
    my @fusion_candidates = &parse_chims($chims_described);
    
    @fusion_candidates = reverse sort {$a->{num_reads} <=> $b->{num_reads}} @fusion_candidates;

    my $num_fusion_candidates = scalar(@fusion_candidates);
    print STDERR "Pre-FFPM-filtering of prelim phase-1 candidates: $num_fusion_candidates fusion pairs\n";
    
    &write_candidates_summary("$output_prefix.preliminary_candidates_info_from_chims_described", \@fusion_candidates);
        
    @fusion_candidates = &filter_chims(\@fusion_candidates, $num_total_reads, $min_num_LR, $min_FFPM, $MAX_EXON_DELTA, $max_foldback_frac);

    $num_fusion_candidates = scalar(@fusion_candidates);
    
    print STDERR "Post-FFPM-filtering of prelim phase-1 candidates: $num_fusion_candidates fusion pairs\n";

    &write_candidates_summary("$output_prefix.preliminary_candidates_info_from_chims_described.read_support_filtered", \@fusion_candidates);

    exit(0);
        
}
        



####
sub parse_chims {
    my ($chims_described_file) = @_;
    
    my %fusion_pairs;

    open (my $fh, $chims_described_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; } # header or comment
        chomp;
        my $line = $_;

        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        my $fusion_info = $x[3];
        
        ## No mitochondrial targets
        if ($fusion_info =~ /chrM:/) {
            next;
        }
        
        my ($geneA, $deltaA, $trans_brkptA, 
            $chrA_n_coordA,
            $geneB, $deltaB, $trans_brkptB, 
            $chrB_n_coordB,
            $fusion_name, $foldback_flag) = split(/;/, $fusion_info);


        # foldback_flag may not be defined in older chims_described files
        $foldback_flag = "" unless defined($foldback_flag);
        
        # Track fold-back reads (overlapping alignments on opposite strands)
        # These are potential artifacts but we count them rather than filter them out
        my $is_foldback = ($foldback_flag eq "FOLDBACK") ? 1 : 0;
        
        my $fusion_info_struct = $fusion_pairs{$fusion_name};
        
        if (! defined $fusion_info_struct) {
            
            # init struct
            $fusion_info_struct = $fusion_pairs{$fusion_name} = { fusion_name => $fusion_name,
                                                                  deltaA => [],
                                                                  deltaB => [],
                                                                  read_names => [],
                                                                  num_reads => 0,
                                                                  num_foldback_reads => 0,
                                                                  trans_brkpt_delta => [],
            };
        }
        
        
        push (@{$fusion_info_struct->{deltaA}}, $deltaA);
        push (@{$fusion_info_struct->{deltaB}}, $deltaB);
        push (@{$fusion_info_struct->{read_names}}, $trans_acc);
        # Track fold-back reads separately to help identify potential artifacts
        $fusion_info_struct->{num_reads}++;
        
        if ($is_foldback) {
            $fusion_info_struct->{num_foldback_reads}++;
        }

        my $trans_brkpt_delta = abs($trans_brkptB - $trans_brkptA);
        push(@{$fusion_info_struct->{trans_brkpt_delta}}, $trans_brkpt_delta);
        
        
    }
    close $fh;
    

    my @fusion_candidates = values %fusion_pairs;


    ## compute_mean_deltas:
    foreach my $fusion_info_struct (@fusion_candidates) {

        $fusion_info_struct->{median_deltaA} = &compute_median_val(@{$fusion_info_struct->{deltaA}});
        $fusion_info_struct->{median_deltaB} = &compute_median_val(@{$fusion_info_struct->{deltaB}});

        $fusion_info_struct->{min_deltaA} = min(@{$fusion_info_struct->{deltaA}});
        $fusion_info_struct->{min_deltaB} = min(@{$fusion_info_struct->{deltaB}});

        $fusion_info_struct->{median_trans_brkpt_delta} = &compute_median_val(@{$fusion_info_struct->{trans_brkpt_delta}});
        $fusion_info_struct->{min_trans_brkpt_delta} = min(@{$fusion_info_struct->{trans_brkpt_delta}});
        
    }
    
    
    return(@fusion_candidates);

}

####
sub compute_mean_val {
    my @vals = @_;
    my $num_vals = scalar @vals;
    my $sum = 0;
    foreach my $val (@vals) {
        $sum += $val;
    }

    my $mean_val = $sum / $num_vals;
    
    return($mean_val);
}

####
sub compute_median_val {
    my @vals = @_;
    my $num_vals = scalar @vals;

    if ($num_vals == 1) {
        return($vals[0]);
    }

    my $midpt = int($num_vals/2);
    if ($num_vals % 2 == 0) {

        return(compute_mean_val($vals[$midpt-1], $vals[$midpt]));
    }
    else {
        return($vals[$midpt]);
    }
    
}


####
sub filter_chims {
    my ($fusion_candidates_aref, $num_total_reads, $min_num_LR, $min_FFPM, $MAX_EXON_DELTA, $max_foldback_frac) = @_;


    my @fusion_candidates;
    
    foreach my $fusion_candidate (@$fusion_candidates_aref) {
        
        my $num_reads = $fusion_candidate->{num_reads};
        my $ffpm = $num_reads / $num_total_reads * 1e6;

        my $min_alt_exon_delta = min($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        my $max_alt_exon_delta = max($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        
        # Calculate fold-back fraction
        my $foldback_frac = ($num_reads > 0) ? $fusion_candidate->{num_foldback_reads} / $num_reads : 0;
        
        if ($ffpm >= $min_FFPM
            &&
            $num_reads >= $min_num_LR
            &&
            # Filter out fusions with too many fold-back reads
            $foldback_frac <= $max_foldback_frac
            &&
            (
             # deltas within range on both sides
             ($fusion_candidate->{min_deltaA} <= $MAX_EXON_DELTA && $fusion_candidate->{min_deltaB} <= $MAX_EXON_DELTA)
                ||

             # deltas within range on either side but need more than one read as evidence.
                ( ($min_alt_exon_delta <= $MAX_EXON_DELTA && $max_alt_exon_delta <= $ALT_MAX_EXON_DELTA) &&  $num_reads > 1 )
              )
            ) {
            push(@fusion_candidates, $fusion_candidate);
        }
    }
    
    return(@fusion_candidates);

}


####
sub write_candidates_summary {
    my ($prelim_candidates_summary_outfile, $fusion_candidates_aref) = @_;

    
    
    open(my $ofh, ">$prelim_candidates_summary_outfile") or die $!;
        
    my $header = join("\t", "#FusionName",
                      "median_deltaA", "median_deltaB",
                      "min_deltaA", "min_deltaB",
                      "median_trans_brkpt_delta",
                      "min_trans_brkpt_delta",
                      "num_reads",
                      "num_foldback_reads");
    
    print $ofh "$header\n";
    
        

    foreach my $fusion_info_struct (@$fusion_candidates_aref) {


            my $outline = join("\t", $fusion_info_struct->{fusion_name}, 

                        int($fusion_info_struct->{median_deltaA} + 0.5),
                        int($fusion_info_struct->{median_deltaB} + 0.5),

                        $fusion_info_struct->{min_deltaA},
                        $fusion_info_struct->{min_deltaB},
                        
                        int($fusion_info_struct->{median_trans_brkpt_delta} + 0.5),
                        $fusion_info_struct->{min_trans_brkpt_delta},

                        $fusion_info_struct->{num_reads},
                        $fusion_info_struct->{num_foldback_reads});
            
            print $ofh "$outline\n";

    
    }

    close $ofh;

    return;
}

