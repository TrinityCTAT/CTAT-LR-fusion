#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Pipeliner;

my $usage = "\n\n\tusage: $0 mm2.map.gff3.chims_described mm2.map.gff3.chims_described.fasta EXTEND_LENGTH genome_lib_dir min_per_id\n\n";

my $chims_described_file = $ARGV[0] or die $usage;
my $chims_fasta_file = $ARGV[1] or die $usage;
my $EXTEND = $ARGV[2] or die $usage;
my $genome_lib_dir = $ARGV[3] or die $usage;
my $min_per_id = $ARGV[4] or die $usage;

## configuration:
my $GENOME = "$genome_lib_dir/ref_genome.fa";
my $MM2_DB_DIR = "$genome_lib_dir";
my $MM2_DB_NAME = "ref_genome.fa.mm2";
my $MM2_idx = "$MM2_DB_DIR/$MM2_DB_NAME/ref_genome.fa.mmi";
my $REF_GTF = "$genome_lib_dir/ref_annot.gtf";
my $MM2_splice_file = "$REF_GTF.mm2.splice.bed";

## make cmd line args
my $max_intron_length = 100000;
my $CPU = 4;


my $UTILDIR = $FindBin::Bin;

main: {

    $min_per_id = $min_per_id/100;
    
    my %transcript_to_breakpoint = &parse_chimera_preds($chims_described_file);

    my $fasta_reader = new Fasta_reader($chims_fasta_file);
    my %trans_seqs = $fasta_reader->retrieve_all_seqs_hash(%transcript_to_breakpoint);

    my $chim_frag_file = "$chims_fasta_file.split.fa";
    open (my $ofh, ">$chim_frag_file");
    
    foreach my $trans (keys %transcript_to_breakpoint) {
        
        my $sequence = $trans_seqs{$trans} or die "Error, no sequence for trans: $trans";
        
        my $breakpoint_range = $transcript_to_breakpoint{$trans} or die "Error, no breakpoint range for $trans";
        
        my ($brk_left, $brk_right) = sort {$a<=>$b} split(/-/, $breakpoint_range);

        my $seq_range_left = substr($sequence, 0, $brk_left + $EXTEND) or die "Error, no substr for range left";
        my $seq_range_right = substr($sequence, $brk_right - $EXTEND) or die "Error, no substr for range right";

        print $ofh ">$trans" . "____left\n"
            . "$seq_range_left\n"
            . ">$trans" . "____right\n"
            . "$seq_range_right\n";
        
    }
    close $ofh;


    my $pipeliner = new Pipeliner(-verbose => 1);
    ## run MM2, capture all top hits within reason.
    my $mm2_output_prefix = "$chim_frag_file.mm2";
    my $cmd = "minimap2 -ax splice --junc-bed $MM2_splice_file -O6,24 -B4 -L -t $CPU -cs -ub -G $max_intron_length $MM2_idx $chim_frag_file > $mm2_output_prefix.sam";
    $pipeliner->add_commands(new Command($cmd, "mm2_chim_frags.ok"));
    
    $cmd = "$UTILDIR/SAM_to_gff3.minimap2.pl  $mm2_output_prefix.sam >  $mm2_output_prefix.gff3";
    $pipeliner->add_commands(new Command($cmd, "mm2_chim_frags_gff3.ok"));
    
    $pipeliner->run();
    
    
    exit(0);
}


####
sub parse_chimera_preds {
    my ($chims_described_file) = @_;

    my %trans_to_brk;

    open (my $fh, $chims_described_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        my $trans_acc = $x[0];
        my $info = $x[3];
        my @pts = split(/;/, $info);
        my $brk_left = $pts[2];
        my $brk_right = $pts[6];

        $trans_to_brk{$trans_acc} = "$brk_left-$brk_right";

    }

    close $fh;

    return(%trans_to_brk);
}
        
