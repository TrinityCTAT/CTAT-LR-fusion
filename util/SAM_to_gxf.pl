#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Carp;
use List::Util qw(min max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);




my $usage = <<__EOUSAGE__;

####################################################################
#
#  --sam <string>      sam or bam file with transcript alignments 
# 
#  --format <string>   gff3 | gtf 
#
#  --allow_non_primary    report non-primary alignments (by default, only reports primary alignments)
#
#  --debug|d
#
####################################################################


__EOUSAGE__

    ;


my $sam_file;
my $format;
my $help_flag;
my $DEBUG = 0;
my $allow_non_primary = 0;

&GetOptions ( 'help|h' => \$help_flag,
              'format=s' => \$format,
              'sam=s' => \$sam_file,
              'allow_non_primary' => \$allow_non_primary,
              'debug|d' => \$DEBUG);


if ($help_flag) {
    die $usage;
}

unless ($sam_file && $format) {
    
    die $usage;

}

unless ($format =~ /^(gff3|gtf)$/) {
    die "Sorry, only gff3 or gtf is allowed for --format";
}


main: {

    my %PATH_COUNTER;
    
	my $sam_reader = new SAM_reader($sam_file);

	while ($sam_reader->has_next()) {

		my $sam_entry = $sam_reader->get_next();

        if ( (! $sam_entry->is_primary() ) && ( ! $allow_non_primary) )  {
            next;
        }
        
		if ($sam_entry->is_query_unmapped()) {
			next;
		}

		my $read_name = $sam_entry->get_read_name();
        if ($read_name =~ /\.p\d$/) {
            # not the first path reported.
            next;
        }
                
        my $sequence = $sam_entry->get_sequence();
        if ($sequence eq "*") {
            next;
        }
        
        my $sam_line = $sam_entry->get_original_line();

        if ($DEBUG) {
            print "$sam_line\n";
        }


        my $NM = 0; # full edit distance (includes mismatches and indels)
                
        if ($sam_line =~ /NM:i:(\d+)/) {
            $NM = $1;
        }
        else {
            die "Error, couldn't extract num mismatches from sam line: $sam_line";
        }
        my $cigar_align = $sam_entry->get_cigar_alignment();
        my $num_indel_nts = 0;
        while($cigar_align =~ /(\d+)[DI]/g) {
            $num_indel_nts += $1;
        }

        my $num_mismatches = $NM - $num_indel_nts;

        if ($num_mismatches < 0) {
            confess "Error, calculated negative mismatch count from: NM:$NM, indel:$num_indel_nts, cigar: $cigar_align";
        }
        

		my $scaff_name = $sam_entry->get_scaffold_name();
		
		my $strand = $sam_entry->get_query_strand();

        
		my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();

        my $min_coord = 0 + 'inf';
        my $max_coord = 0 - 'inf';

        my $align_len = 0;
        {
            foreach my $coordset (@$genome_coords_aref) {
                my $seglen = abs($coordset->[1] - $coordset->[0]) + 1;
                $align_len += $seglen;

                $min_coord = min($min_coord, $coordset->[0], $coordset->[1]);
                $max_coord = max($max_coord, $coordset->[0], $coordset->[1]);

                if ($DEBUG) {
                    print STDERR join("\t", $coordset->[0], $coordset->[1], "seglen: $seglen") . "\n";
                }
            }
        }


        
        if ($DEBUG) {
            print STDERR "num_mismatches: $num_mismatches, align_length: $align_len\n";
        }
        
        my $per_id = sprintf("%.1f", 100 - $num_mismatches/$align_len * 100); 

        my $align_counter = "$read_name.p" . ++$PATH_COUNTER{$read_name};


        my @genome_n_trans_coords;

                
        while (@$genome_coords_aref) {
            my $genome_coordset_aref = shift @$genome_coords_aref;
            my $trans_coordset_aref = shift @$query_coords_aref;

            my ($genome_lend, $genome_rend) = @$genome_coordset_aref;
            my ($trans_lend, $trans_rend) = sort {$a<=>$b} @$trans_coordset_aref;

            push (@genome_n_trans_coords, [ $genome_lend, $genome_rend, $trans_lend, $trans_rend ] );

        }

        ## merge neighboring features if within a short distance unlikely to represent an intron.
        my @merged_coords;
        push (@merged_coords, shift @genome_n_trans_coords);

        my $MERGE_DIST = 10;
        while (@genome_n_trans_coords) {
            my $coordset_ref = shift @genome_n_trans_coords;
            my $last_coordset_ref = $merged_coords[$#merged_coords];
            
            if ($coordset_ref->[0] - $last_coordset_ref->[1] <= $MERGE_DIST) {
                # merge it.
                $last_coordset_ref->[1] = $coordset_ref->[1];

                if ($strand eq "+") {
                    $last_coordset_ref->[3] = $coordset_ref->[3];
                } else {
                    $last_coordset_ref->[2] = $coordset_ref->[2];
                }
            }
            else {
                # not merging.
                push (@merged_coords, $coordset_ref);
            }
        }

        #my $trans_align_len = 0;
        #oreach my $coordset_ref (@merged_coords) {
        #    my ($genome_lend, $genome_rend, $trans_lend, $trans_rend) = @$coordset_ref;
        #    $trans_align_len += $trans_rend - $trans_lend + 1;
        #}
        #
        #if ($DEBUG) {
        #    print "interval-based alignment length: $trans_align_len\n";
        #
        
        
        if ($format eq "gtf") {
            # need a transcript record.
            @merged_coords = sort {$a->[0]<=>$b->[0]} @merged_coords;
            my $genome_lend = $merged_coords[0]->[0];
            my $genome_rend = $merged_coords[-1]->[1];
            my $trans_lend = $merged_coords[0]->[2];
            my $trans_rend = $merged_coords[-1]->[3];
            
            print join("\t",
                       $scaff_name,
                       "minimap2",
                       "transcript",
                       $genome_lend, $genome_rend,
                       $per_id,
                       $strand,
                       ".",
                       "gene_id \"$align_counter.g\"; transcript_id \"$align_counter.mrna\"; Target \"$read_name $trans_lend $trans_rend\";"
                ) . "\n"; 
        }
        
        foreach my $coordset_ref (@merged_coords) {
            my ($genome_lend, $genome_rend, $trans_lend, $trans_rend) = @$coordset_ref;

            
            my $feat_type = ($format eq "gff3") ? "cDNA_match" : "exon";
            
            my $info_field = ($format eq "gff3") ?
                "ID=$align_counter;Parent=$align_counter.mrna;Target=$read_name $trans_lend $trans_rend" :
                "gene_id \"$align_counter.g\"; transcript_id \"$align_counter.mrna\"; Target \"$read_name $trans_lend $trans_rend\";";
            

            print join("\t",
                       $scaff_name,
                       "minimap2",
                       $feat_type,
                       $genome_lend, $genome_rend,
                       $per_id,
                       $strand,
                       ".",
                       $info_field
                ) . "\n";
        }
        print "\n";
        
        
        
	}


	exit(0);
}
