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
# --chims_described <string>   chims described file

# --reads <string>             fastA or fastQ file containing reads
#
# --fusions <string>           fusion candidates input file
#
# --output_prefix <string>    prefix name for output files (prefix).transcripts.fa and (prefix).FI_listing
#
# --skip_read_extraction      dont extract the fusion reads, just generate the prelim report.
#
###########################################################################################################


__EOUSAGE__

    ;


my $chims_described_file;
my $reads_file;
my $fusions_input_file;
my $help_flag;
my $output_prefix;
my $SKIP_READ_EXTRACTION = 0;

&GetOptions ( 'help|h' => \$help_flag,
              'chims_described=s' => \$chims_described_file,
              'reads=s' => \$reads_file,
              'fusions=s' => \$fusions_input_file,
              'output_prefix=s' => \$output_prefix,
              'skip_read_extraction' => \$SKIP_READ_EXTRACTION,
    );

if ($help_flag) {
    die $usage;
}

unless($chims_described_file) {
    print STDERR "\n\nERROR - must specify --chims_described <string> \n";
}

unless ($reads_file) {
    print STDERR "\n\nERROR - must specify --reads <string> \n";
    die $usage;
}

unless ($fusions_input_file) {
    print STDERR "\n\nERROR - must specify --fusions <string>\n";
    die $usage;
}

unless($output_prefix) {
    print STDERR "\n\nERROR - must specify --output_prefix <string>\n";
    die $usage;
}

$reads_file = &ensure_full_path($reads_file);
$fusions_input_file = &ensure_full_path($fusions_input_file);

foreach my $file ($reads_file, $fusions_input_file) {
    unless (-s $file) {
        confess "Error, cannot locate file $file";
    }
}


main: {

    my %fusion_targets = &parse_fusion_targets($fusions_input_file);
    
    my @fusion_candidates = &parse_chims($chims_described_file, \%fusion_targets);
    
    @fusion_candidates = reverse sort {$a->{num_reads} <=> $b->{num_reads}} @fusion_candidates;

    # prep reads info
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

    print STDERR "Prelim phase-1 candidates: $num_fusion_candidates fusion pairs involving $num_fusion_candidate_reads reads.\n";

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
    my ($chims_described_file, $fusion_targets_href) = @_;
    
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

        if (! exists $fusion_targets_href->{$fusion_name}) {
            next;
        }
        
        my $fusion_info_struct = $fusion_pairs{$fusion_name};
        
        if (! defined $fusion_info_struct) {
            
            # init struct
            $fusion_info_struct = $fusion_pairs{$fusion_name} = { fusion_name => $fusion_name,
                                                                  read_names => [],
                                                                  num_reads => 0,
                                                                  
            };
        }
        
        
        push (@{$fusion_info_struct->{read_names}}, $trans_acc);
        $fusion_info_struct->{num_reads}++;

        
    }
    close $fh;
    

    my @fusion_candidates = values %fusion_pairs;

    return(@fusion_candidates);

}



####
sub write_candidates_summary {
    my ($prelim_candidates_summary_outfile, $fusion_candidates_aref) = @_;

    
    
    open(my $ofh, ">$prelim_candidates_summary_outfile") or die $!;
    open(my $ofh_wreads, ">$prelim_candidates_summary_outfile.with_reads") or die $!;
    
    my $header = join("\t", "#FusionName",
                      "num_reads");
    
    print $ofh "$header\n";
    print $ofh_wreads "#FusionName\tread\n";

    foreach my $fusion_info_struct (@$fusion_candidates_aref) {

        my $fusion_name = $fusion_info_struct->{fusion_name};
        my $num_reads = $fusion_info_struct->{num_reads};
        
        print $ofh "$fusion_name\t$num_reads\n";
        my @reads  = @{$fusion_info_struct->{read_names}};
        foreach my $read (@reads) {
            print $ofh_wreads "$fusion_name\t$read\n";
        }
    }
    
    close $ofh;
    close $ofh_wreads;
    
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



####
sub parse_fusion_targets {
    my ($fusions_file) = @_;

    my %fusion_targets;
    
    open(my $fh, $fusions_file) or die "Error, cannot open file $fusions_file";
    my $header = <$fh>;
    
    while(<$fh>) {
        chomp;
        my @vals = split(/\t/);
        my $fusion_name = $vals[0];
        $fusion_targets{$fusion_name} = 1;
    }
    close $fh;

    return(%fusion_targets);
}
