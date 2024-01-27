#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use File::Basename;
use Process_cmd;
use Pipeliner;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use List::Util qw(min max);

my $usage = <<__EOUSAGE__;


#########################################################################################################
#
# --trans_fasta <string>     transcripts.fasta for long reads
#
# --chims_described <string>   chims.descripbed file.
#
# --max_exon_delta <int>       maximum dist from ref exon boundary
#
# --output_prefix <string>    prefix name for output files (prefix).transcripts.fa and (prefix).FI_listing
#
# --min_FFPM <float>          min fusion expression for candidates to pursue
#
# --skip_read_extraction      dont extract the fusion reads, just generate the prelim report.
#
###########################################################################################################


__EOUSAGE__

    ;


my $trans_fasta;
my $chims_described;
my $MAX_EXON_DELTA;
my $help_flag;
my $output_prefix;
my $min_FFPM;
my $SKIP_READ_EXTRACTION = 0;


my $ALT_MAX_EXON_DELTA = 1000;

&GetOptions ( 'help|h' => \$help_flag,
              'trans_fasta=s' => \$trans_fasta,
              'chims_described=s' => \$chims_described,
              'max_exon_delta=i' => \$MAX_EXON_DELTA,
              'output_prefix=s' => \$output_prefix,
              'min_FFPM=f' => \$min_FFPM,
              'skip_read_extraction' => \$SKIP_READ_EXTRACTION,
    );

if ($help_flag) {
    die $usage;
}

unless ($trans_fasta) {
    print STDERR "\n\nERROR - must specify --trans_fasta <string> \n";
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

unless (defined ($min_FFPM) ) {
    print STDERR "\n\nError, must set --min_FFPM ";
    die $usage;
}


$trans_fasta = &ensure_full_path($trans_fasta);
$chims_described = &ensure_full_path($chims_described);


foreach my $file ($trans_fasta, $chims_described) {
    unless (-s $file) {
        confess "Error, cannot locate file $file";
    }
}


main: {


    # get total read count.
    my $TOTAL_READS = `grep '>' $trans_fasta | wc -l `;
    if ($?) { 
        die "Error, cmd: \"grep '>' $trans_fasta | wc -l  \" died with ret $?";
    }
    chomp $TOTAL_READS;
    $TOTAL_READS = int($TOTAL_READS);
    unless ($TOTAL_READS > 0) {
        die "Error, could not count number of reads from file: $trans_fasta (shouldnt happen....)";
    }
    
    
    my @fusion_candidates = &parse_chims($chims_described);
    
    @fusion_candidates = reverse sort {$a->{num_reads} <=> $b->{num_reads}} @fusion_candidates;
    
    my $prelim_candidates_summary_outfile = "$output_prefix.preliminary_candidates_info";
    &write_candidates_summary($prelim_candidates_summary_outfile, \@fusion_candidates);
        
    @fusion_candidates = &filter_chims(\@fusion_candidates, $TOTAL_READS, $min_FFPM, $MAX_EXON_DELTA);
    
    

    
    
    # store the total number of reads.
    open (my $ofh, ">$trans_fasta.LR_read_count.txt") or die $!;
    print $ofh "$TOTAL_READS\n";
    close $ofh;
    
    
    my $num_fusion_candidates = scalar(@fusion_candidates);
    my $num_fusion_candidate_reads = 0;
    my %reads_want;
    foreach my $fusion_candidate (@fusion_candidates) {
        $num_fusion_candidate_reads += $fusion_candidate->{num_reads};
        unless($SKIP_READ_EXTRACTION) {
            foreach my $read (@{$fusion_candidate->{read_names}}) {
                $reads_want{$read} = 1;
            }
        }
    }
    
    print STDERR "Post-FFPM-filtering of prelim phase-1 candidates: $num_fusion_candidates fusion pairs involving $num_fusion_candidate_reads reads.\n";

    &write_candidates_summary("$output_prefix.FI_listing", \@fusion_candidates);


    if ($SKIP_READ_EXTRACTION) {
        print STDERR "-skipping read extraction and stopping here. See prelim candidates: $output_prefix.FI_listing\n\n";
        exit(0);
    }

    open(my $ofh_fasta, ">$output_prefix.transcripts.fa") or die $!;
    
    my $fasta_reader = new Fasta_reader($trans_fasta);
    while (my $seq_obj = $fasta_reader->next()) {
        my $accession = $seq_obj->get_accession();

        if (exists $reads_want{$accession}) {
            my $sequence = $seq_obj->get_sequence();
            print $ofh_fasta ">$accession\n$sequence\n";
            delete $reads_want{$accession};
        }
    }
    close $ofh_fasta;
    
    if (%reads_want) {
        confess "Error, missing some reads during extraction: " . Dumper(\%reads_want);
    }
    
    print STDERR "-done. See files: $output_prefix.transcripts.fa and $output_prefix.FI_listing\n";
        
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
            $fusion_name) = split(/;/, $fusion_info);


        
        my $fusion_info_struct = $fusion_pairs{$fusion_name};
        
        if (! defined $fusion_info_struct) {
            
            # init struct
            $fusion_info_struct = $fusion_pairs{$fusion_name} = { fusion_name => $fusion_name,
                                                                  deltaA => [],
                                                                  deltaB => [],
                                                                  read_names => [],
                                                                  num_reads => 0 };
        }
        
        
        push (@{$fusion_info_struct->{deltaA}}, $deltaA);
        push (@{$fusion_info_struct->{deltaB}}, $deltaB);
        push (@{$fusion_info_struct->{read_names}}, $trans_acc);
        $fusion_info_struct->{num_reads}++;
        
    }
    close $fh;
    

    my @fusion_candidates = values %fusion_pairs;


    ## compute_mean_deltas:
    foreach my $fusion_info_struct (@fusion_candidates) {
        $fusion_info_struct->{mean_deltaA} = &compute_mean_val(@{$fusion_info_struct->{deltaA}});
        $fusion_info_struct->{mean_deltaB} = &compute_mean_val(@{$fusion_info_struct->{deltaB}});

        $fusion_info_struct->{median_deltaA} = &compute_median_val(@{$fusion_info_struct->{deltaA}});
        $fusion_info_struct->{median_deltaB} = &compute_median_val(@{$fusion_info_struct->{deltaB}});

        $fusion_info_struct->{min_deltaA} = min(@{$fusion_info_struct->{deltaA}});
        $fusion_info_struct->{min_deltaB} = min(@{$fusion_info_struct->{deltaB}});
        
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
    my ($fusion_candidates_aref, $TOTAL_READS, $min_FFPM, $MAX_EXON_DELTA) = @_;


    my @fusion_candidates;
    
    foreach my $fusion_candidate (@$fusion_candidates_aref) {
        
        my $num_reads = $fusion_candidate->{num_reads};
        my $ffpm = $num_reads / $TOTAL_READS * 1e6;

        my $min_alt_exon_delta = min($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        my $max_alt_exon_delta = max($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        
        if ($ffpm >= $min_FFPM
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
    print $ofh join("\t", "#FusionName",
                    "median_deltaA", "median_deltaB",
                    "min_deltaA", "min_deltaB",
                    "mean_deltaA", "mean_deltaB",
                    "num_reads") . "\n";

    foreach my $fusion_info_struct (@$fusion_candidates_aref) {

        print $ofh join("\t", $fusion_info_struct->{fusion_name}, 

                        int($fusion_info_struct->{median_deltaA} + 0.5),
                        int($fusion_info_struct->{median_deltaB} + 0.5),

                        $fusion_info_struct->{min_deltaA},
                        $fusion_info_struct->{min_deltaB},
                        
                        int($fusion_info_struct->{mean_deltaA} + 0.5),
                        int($fusion_info_struct->{mean_deltaB} + 0.5),

                        $fusion_info_struct->{num_reads}) . "\n";
    }

    close $ofh;

    return;
}
