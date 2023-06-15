#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Process_cmd;

my $usage = "\n\n\tusage: $0 gmap.map.gff3.chims_described.fasta gmap.map.gff3.chims_described.w_read_support.J1.JS2.NJ3 gmap.map.gff3.chims_described.fasta.bowtie2.bam\n\n";


my $chims_fasta = $ARGV[0] or die $usage;
my $fusion_preds = $ARGV[1] or die $usage;
my $bam_file = $ARGV[2] or die $usage;

main: {

    unless (-e "$chims_fasta.fai") {
        &process_cmd("samtools faidx $chims_fasta");
    }

    my %junction_reads;
    my %spanning_frags;
    my %trans_breakpoint;
    &parse_fusion_evidence($fusion_preds, \%junction_reads, \%spanning_frags, \%trans_breakpoint);

    &write_fusion_evidence_bams($chims_fasta, $bam_file, \%junction_reads, \%spanning_frags);
    
    &write_breakpoint_gff(\%trans_breakpoint);


    exit(0);

}

####
sub write_breakpoint_gff {
    my ($trans_breakpoint_href) = @_;
    
    open (my $ofh, ">breakpoints.gff") or die $!;

    foreach my $trans_acc (keys %$trans_breakpoint_href) {
        
        my ($breakpoint, $fusion_name) = split(/\t/, $trans_breakpoint_href->{$trans_acc});
        
        my ($brk_lend, $brk_rend) = sort {$a<=>$b} split(/-/, $breakpoint);

        print $ofh join("\t", $trans_acc, "breakpoint", "match", $brk_lend, $brk_rend, ".", "+", ".", "ID=\"$fusion_name|$trans_acc.brkpt:$brk_lend-$brk_rend\"") . "\n";
    }

    
    close $ofh;

    return;
}



####
sub write_fusion_evidence_bams {
    my ($chims_fasta, $bam_file, $junction_reads_href, $spanning_frags_href) = @_;

    my $junction_frags_sam_file = "junction_reads.sam";
    my $spanning_frags_sam_file = "spanning_frags.sam";

    open (my $ofh_junc, ">$junction_frags_sam_file") or die $!;
    open (my $ofh_span, ">$spanning_frags_sam_file") or die $!;

    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {
        my $trans_acc = $sam_entry->get_scaffold_name();
        my $frag_name = $sam_entry->get_read_name();
        
        my $token = join("$;", $trans_acc, $frag_name);

        if (exists $junction_reads_href->{$token}) {
            print $ofh_junc $sam_entry->get_original_line() . "\n";
        }
        elsif (exists $spanning_frags_href->{$token}) {
            print $ofh_span $sam_entry->get_original_line() . "\n";
        }
    }

    close $ofh_junc;
    close $ofh_span;

    
    ## convert to coord-sorted bam files.
    
    my $cmd = "samtools view -T $chims_fasta -Sb $junction_frags_sam_file | samtools sort -o $junction_frags_sam_file.csorted.bam";
    &process_cmd($cmd);
    
    &process_cmd("samtools index $junction_frags_sam_file.csorted.bam");


    $cmd =  "samtools view -T $chims_fasta -Sb $spanning_frags_sam_file | samtools sort -o $spanning_frags_sam_file.csorted.bam"; 
    &process_cmd($cmd);
    
    &process_cmd("samtools index $spanning_frags_sam_file.csorted.bam");
    
    unlink($junction_frags_sam_file, $spanning_frags_sam_file);
    
    return;

}




####
sub parse_fusion_evidence {
    my ($fusion_preds, $junction_reads_href, $spanning_frags_href, $trans_breakpoint_href) = @_;

    open (my $fh, $fusion_preds) or die "Error, cannot open file $fusion_preds";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        
        my $fusion_info = $x[3];
        my @pts = split(/;/, $fusion_info);
        my $brkpt = join("-", $pts[2], $pts[6]);
        my $fusion_name = pop @pts;
        $trans_breakpoint_href->{$trans_acc} = "$brkpt\t$fusion_name";

        my $junction_read_list = $x[10];
        my $spanning_frag_list = $x[11];
        
        

        foreach my $junction_read (split(/,/, $junction_read_list)) {
            my $token = join("$;", $trans_acc, $junction_read);

            $junction_reads_href->{$token} = 1;
        }

        foreach my $spanning_frag (split(/,/, $spanning_frag_list) ) {
            my $token = join("$;", $trans_acc, $spanning_frag);
            
            $spanning_frags_href->{$token} = 1;
        }

    }
    close $fh;

    
    return;
}
