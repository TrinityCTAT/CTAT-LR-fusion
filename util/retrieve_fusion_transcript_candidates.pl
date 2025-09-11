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
# --reads <string>             fastA or fastQ file containing reads
#
# --chims_described <string>   chims.descripbed file.
#
# --max_exon_delta <int>       maximum dist from ref exon boundary
#
# --output_prefix <string>    prefix name for output files (prefix).transcripts.fa and (prefix).FI_listing
#
# --min_num_LR <int>          minimum number of long reads required as evidence.
# 
# --min_FFPM <float>          min fusion expression for candidates to pursue
#
# --skip_read_extraction      dont extract the fusion reads, just generate the prelim report.
#
# --num_total_reads <int>     number of total reads (used for FFPM calculation)
#
###########################################################################################################


__EOUSAGE__

    ;


my $reads_file;
my $chims_described;
my $MAX_EXON_DELTA;
my $help_flag;
my $output_prefix;
my $min_FFPM;
my $SKIP_READ_EXTRACTION = 0;
my $num_total_reads;
my $min_num_LR = 0;

my $ALT_MAX_EXON_DELTA = 1000;

&GetOptions ( 'help|h' => \$help_flag,
              'reads=s' => \$reads_file,
              'chims_described=s' => \$chims_described,
              'max_exon_delta=i' => \$MAX_EXON_DELTA,
              'output_prefix=s' => \$output_prefix,
              'min_FFPM=f' => \$min_FFPM,
              'skip_read_extraction' => \$SKIP_READ_EXTRACTION,
              'num_total_reads=i' => \$num_total_reads,
              'min_num_LR=i' => \$min_num_LR,
    );

if ($help_flag) {
    die $usage;
}

unless ($reads_file) {
    print STDERR "\n\nERROR - must specify --reads <string> \n";
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


$reads_file = &ensure_full_path($reads_file);
$chims_described = &ensure_full_path($chims_described);


foreach my $file ($reads_file, $chims_described) {
    unless (-s $file) {
        confess "Error, cannot locate file $file";
    }
}


main: {
    
    my @fusion_candidates = &parse_chims($chims_described);
    
    @fusion_candidates = reverse sort {$a->{num_reads} <=> $b->{num_reads}} @fusion_candidates;
    
    my $prelim_candidates_summary_outfile = "$output_prefix.preliminary_candidates_info";
    &write_candidates_summary($prelim_candidates_summary_outfile, \@fusion_candidates);
        
    @fusion_candidates = &filter_chims(\@fusion_candidates, $num_total_reads, $min_num_LR, $min_FFPM, $MAX_EXON_DELTA);
    
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

    my $reads_file_type = &get_reads_file_type($reads_file);

    my $reader = ($reads_file_type eq "FASTA") ? new Fasta_reader($reads_file) : new Fastq_reader($reads_file);
    
    while (my $seq_obj = $reader->next()) {
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
                                                                  num_reads => 0,
                                                                  trans_brkpt_delta => [],
            };
        }
        
        
        push (@{$fusion_info_struct->{deltaA}}, $deltaA);
        push (@{$fusion_info_struct->{deltaB}}, $deltaB);
        push (@{$fusion_info_struct->{read_names}}, $trans_acc);
        $fusion_info_struct->{num_reads}++;

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
    my ($fusion_candidates_aref, $num_total_reads, $min_num_LR, $min_FFPM, $MAX_EXON_DELTA) = @_;


    my @fusion_candidates;
    
    foreach my $fusion_candidate (@$fusion_candidates_aref) {
        
        my $num_reads = $fusion_candidate->{num_reads};
        my $ffpm = $num_reads / $num_total_reads * 1e6;

        my $min_alt_exon_delta = min($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        my $max_alt_exon_delta = max($fusion_candidate->{min_deltaA}, $fusion_candidate->{min_deltaB});
        
        if ($ffpm >= $min_FFPM
            &&
            $num_reads >= $min_num_LR
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
    open(my $ofh_wreads, ">$prelim_candidates_summary_outfile.with_reads") or die $!;
    
    my $header = join("\t", "#FusionName",
                      "median_deltaA", "median_deltaB",
                      "min_deltaA", "min_deltaB",
                      "median_trans_brkpt_delta",
                      "min_trans_brkpt_delta",
                      "num_reads");
    
    print $ofh "$header\n";
    print $ofh_wreads "$header\treads\n";

        

    foreach my $fusion_info_struct (@$fusion_candidates_aref) {


            my $outline = join("\t", $fusion_info_struct->{fusion_name}, 

                        int($fusion_info_struct->{median_deltaA} + 0.5),
                        int($fusion_info_struct->{median_deltaB} + 0.5),

                        $fusion_info_struct->{min_deltaA},
                        $fusion_info_struct->{min_deltaB},
                        
                        int($fusion_info_struct->{median_trans_brkpt_delta} + 0.5),
                        $fusion_info_struct->{min_trans_brkpt_delta},

                        $fusion_info_struct->{num_reads});
            
            print $ofh "$outline\n";

            my $reads = join(",", @{$fusion_info_struct->{read_names}});

            print $ofh_wreads "$outline\t$reads\n";
            

    }

    close $ofh;

    return;
}


####
sub get_reads_file_type {
    my ($reads_file) = @_;

    if ($reads_file =~ /\.(fasta|fa)(\.gz)?$/i) {
        return("FASTA");
    }
    elsif ($reads_file =~ /\.(fastq|fq)(\.gz)?/i) {
        return("FASTQ");
    }
    else {
        confess "Error, cannot determine if $reads_file is fastA or fastQ format based on the filename";
    }
}



