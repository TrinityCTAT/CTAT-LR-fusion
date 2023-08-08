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
###########################################################################################################


__EOUSAGE__

    ;





my $trans_fasta;
my $chims_described;
my $MAX_EXON_DELTA;
my $help_flag;
my $output_prefix;
my $min_FFPM;

&GetOptions ( 'help|h' => \$help_flag,
              'trans_fasta=s' => \$trans_fasta,
              'chims_described=s' => \$chims_described,
              'max_exon_delta=i' => \$MAX_EXON_DELTA,
              'output_prefix=s' => \$output_prefix,
              'min_FFPM=f' => \$min_FFPM,
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
    
    my %fusion_pair_read_counts;
    my %fusion_pair_to_read_names;
    
    my %chims = &parse_chims($chims_described, $MAX_EXON_DELTA, \%fusion_pair_read_counts, \%fusion_pair_to_read_names);

    print STDERR "Pre-FFPM-filtering of prelim phase-1 candidates: " . scalar(keys %fusion_pair_read_counts) . " fusion pairs involving " . scalar(keys %chims) . " reads.\n";
    
    &filter_chims_via_FFPM(\%chims, \%fusion_pair_read_counts, \%fusion_pair_to_read_names, $TOTAL_READS, $min_FFPM);
    
    print STDERR "Post-FFPM-filtering of prelim phase-1 candidates: " . scalar(keys %fusion_pair_read_counts) . " fusion pairs involving " . scalar(keys %chims) . " reads.\n";
    
    open(my $ofh_fasta, ">$output_prefix.transcripts.fa") or die $!;
        
    my $fasta_reader = new Fasta_reader($trans_fasta);

    my $read_counter = 0;
    while (my $seq_obj = $fasta_reader->next()) {
        
        $read_counter += 1;

        my $accession = $seq_obj->get_accession();
        
        if (exists $chims{$accession}) {
            
            my $sequence = $seq_obj->get_sequence();
            
            print $ofh_fasta ">$accession\n$sequence\n";
            
            delete $chims{$accession};

        }
    }
    close $ofh_fasta;
    
    
    unless ($read_counter == $TOTAL_READS) {
        die "Error, $read_counter != $TOTAL_READS, so mismatch between read counting results. (bug) ";
    }
    

    if (%chims) {
        die "ERROR, missing sequences for reads: " . Dumper(\%chims);
    }

    # write the fusion listing
    open(my $ofh_FI_list, ">$output_prefix.FI_listing") or die $!;
    foreach my $fusion_pair (sort {$fusion_pair_read_counts{$b} <=> $fusion_pair_read_counts{$a}} keys %fusion_pair_read_counts) {
        print $ofh_FI_list "$fusion_pair\t$fusion_pair_read_counts{$fusion_pair}\n";
    }
    close $ofh_FI_list;

    # store the total number of reads.
    open (my $ofh, ">$trans_fasta.LR_read_count.txt") or die $!;
    print $ofh "$read_counter\n";
    close $ofh;
    

    print STDERR "-done. See files: $output_prefix.transcripts.fa and $output_prefix.FI_listing\n";

    exit(0);
        
}
        



####
sub parse_chims {
    my ($chims_described_file, $MAX_EXON_DELTA, $fusion_pair_read_counts_href, $fusion_pair_to_read_names_href) = @_;

    my %chims;
    
    open (my $fh, $chims_described_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; } # header or comment
        chomp;
        my $line = $_;

        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        my $fusion_info = $x[3];

        my ($geneA, $deltaA, $trans_brkptA, 
            $chrA_n_coordA,
            $geneB, $deltaB, $trans_brkptB, 
            $chrB_n_coordB,
            $fusion_name) = split(/;/, $fusion_info);


	unless ($deltaA <= $MAX_EXON_DELTA && $deltaB <= $MAX_EXON_DELTA) { next; }
	
        my $brkpt_range = join("-", sort ($trans_brkptA, $trans_brkptB));
        
        push (@{$chims{$trans_acc}}, { line => $line,
                                       brkpt_range => $brkpt_range,
              }
            );

        $fusion_pair_read_counts_href->{$fusion_name}++;
        push (@{$fusion_pair_to_read_names_href->{$fusion_name}}, $trans_acc);
    }
    close $fh;
    
    return(%chims);
}


####
sub filter_chims_via_FFPM {
    my ($chims_href, $fusion_pair_read_counts_href, $fusion_pair_to_read_names_href, $TOTAL_READS, $min_FFPM) = @_;

    foreach my $fusion_pair (keys %$fusion_pair_to_read_names_href) {
        my $read_support = $fusion_pair_read_counts_href->{$fusion_pair};
        my $ffpm = $read_support / $TOTAL_READS * 1e6;
        if ($ffpm < $min_FFPM) {
            ## remove candidates.
            my @transcripts = @{$fusion_pair_to_read_names_href->{$fusion_pair}};
            delete $fusion_pair_read_counts_href->{$fusion_pair};
            foreach my $transcript (@transcripts) {
                delete $chims_href->{$transcript};
            }
        }
    }
    
    return;
}
