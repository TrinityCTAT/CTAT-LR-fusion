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

my $usage = "\n\n\tusage: $0 trans.fasta gmap.map.gff3.chims_described\n\n";

my $trans_fasta = $ARGV[0] or die $usage;
my $chims_described = $ARGV[1] or die $usage;

$trans_fasta = &ensure_full_path($trans_fasta);
$chims_described = &ensure_full_path($chims_described);


foreach my $file ($trans_fasta, $chims_described) {
    unless (-s $file) {
        confess "Error, cannot locate file $file";
    }
}


main: {
    
    my %chims = &parse_chims($chims_described);
    
    
    my $fasta_reader = new Fasta_reader($trans_fasta);

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $accession = $seq_obj->get_accession();
        
        if (exists $chims{$accession}) {
            
            my $sequence = $seq_obj->get_sequence();
            
            print ">$accession\n$sequence\n";
            
            delete $chims{$accession};

        }
    }
    

    if (%chims) {
        die "ERROR, missing sequences for accessions: " . keys %chims;
    }

    exit(0);
        

    
}
        



####
sub parse_chims {
    my ($chims_described_file) = @_;

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

        my $brkpt_range = join("-", sort ($trans_brkptA, $trans_brkptB));
        
        push (@{$chims{$trans_acc}}, { line => $line,
                                       brkpt_range => $brkpt_range,
              }
            );
    }
    close $fh;

    return(%chims);
}

