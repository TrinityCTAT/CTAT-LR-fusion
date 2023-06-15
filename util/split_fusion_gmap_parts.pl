#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Pipeliner;

my $usage = "\n\n\tusage: $0 gmap.map.gff3.chims_described gmap.map.gff3.chims_described.fasta EXTEND_LENGTH genome_lib_dir min_per_id\n\n";

my $chims_described_file = $ARGV[0] or die $usage;
my $chims_fasta_file = $ARGV[1] or die $usage;
my $EXTEND = $ARGV[2] or die $usage;
my $genome_lib_dir = $ARGV[3] or die $usage;
my $min_per_id = $ARGV[4] or die $usage;

## configuration:
my $GENOME = "$genome_lib_dir/ref_genome.fa";
my $GMAP_DB_DIR = "$genome_lib_dir";
my $GMAP_DB_NAME = "ref_genome.fa.gmap";


main: {

    $min_per_id = $min_per_id/100;
    
    my %transcript_to_breakpoint = &parse_chimera_preds($chims_described_file);

    my $fasta_reader = new Fasta_reader($chims_fasta_file);
    my %trans_seqs = $fasta_reader->retrieve_all_seqs_hash();

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

    ## run GMAP, capture all top hits within reason.
    my $gmap_output_file = "$chim_frag_file.gmap.gff3";
    my $cmd = "gmap -D $GMAP_DB_DIR -d $GMAP_DB_NAME $chim_frag_file -f 3 -n 1 -t 4 --min-identity $min_per_id > $gmap_output_file";
    
    my $pipeliner = new Pipeliner(-verbose => 1);
    $pipeliner->add_commands(new Command($cmd, "$gmap_output_file.ok"));

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
        
