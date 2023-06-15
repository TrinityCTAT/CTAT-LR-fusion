#!/usr/bin/env perl

package GMAP_gff3_parser;

use strict;
use warnings;


=returns_list_of_align_structs:

 {
            'chr' => 'chr7',
            'per_id' => '100.00',
            'range_lend' => '1',
            'align_id' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464.path1',
            'range_rend' => '2746',
            'rend' => '781105',
            'target' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464',
            'lend' => '588842',
            'orient' => '+',
            'exons' => [
                         {
                           'chr' => 'chr7',
                           'per_id' => 100,
                           'range_lend' => '1',
                           'align_id' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464.path1',
                           'range_rend' => '1398',
                           'target' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464',
                           'rend' => 590239,
                           'lend' => 588842,
                           'orient' => '+'
                         },
                         {
                           'chr' => 'chr7',
                           'per_id' => 100,
                           'range_lend' => '1399',
                           'align_id' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464.path1',
                           'range_rend' => '1480',
                           'target' => 'Locus_1667_Transcript_16/18_Confidence_0.500_Length_3464',
                           'rend' => 591107,
                           'lend' => 591026,
                           'orient' => '+'
                         },
                    ...


=cut

   


####
sub parse_GMAP_gff3_alignments {
    my ($gmap_gff3_file) = @_;
        
    my %target_to_aligns;
    print STDERR "-loading alignment data\n";
    open (my $fh, $gmap_gff3_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }

        chomp;
        my ($chr, $filename, $type, $lend, $rend, $per_id, $orient, $dot, $info) = split(/\t/);
        
        my %info_hash;
        foreach my $keyval (split(/;/, $info)) {
            my ($key, $val) = split(/=/, $keyval);
            $info_hash{$key} = $val;
        }

        my $alignment_ID = $info_hash{ID};
        
        my ($target, $range_lend, $range_rend) = split(/\s+/, $info_hash{Target});

        push (@{$target_to_aligns{$target}->{$alignment_ID}}, { chr => $chr,
                                                                target => $target,
                                                                align_id => $alignment_ID,
                                                                lend => $lend,
                                                                rend => $rend,
                                                                orient => $orient,
                                                                per_id => $per_id,
                                                                
                                                                range_lend => $range_lend,
                                                                range_rend => $range_rend,
                                                            });
        
    }
    close $fh;


    my @span_structs;

    foreach my $target (keys %target_to_aligns) {
        
        ## create spans
        my @span_ids = keys %{$target_to_aligns{$target}};
                
        my @spans;
        foreach my $span_id (@span_ids) {
            my $exon_hits_aref = $target_to_aligns{$target}->{$span_id};
            my $span_struct = &convert_to_span($exon_hits_aref);
            push (@span_structs, $span_struct);
        }
    }
    
    return(@span_structs);
}


####
sub convert_to_span {
    my ($exons_aref) = @_;

    my @chr_coords;
    my @span_coords;
    
    my $orient;
    my $chr;
    my $target;
    my $align_id;
    
    my $sum_per_id = 0;
    my $sum_len = 0;


    foreach my $exon (@$exons_aref) {
        
        # chr_coords
        push (@chr_coords, $exon->{lend});
        push (@chr_coords, $exon->{rend});
       
        # transcript_coords
        push (@span_coords, $exon->{range_lend});
        push (@span_coords, $exon->{range_rend});

        my $len = abs($exon->{rend} - $exon->{lend} + 1);
        $sum_len += $len;
        $sum_per_id += $len * $exon->{per_id};
        

        my $exon_chr = $exon->{chr};
        if ($chr && $exon_chr ne $chr) {
            die "inconsistent chr assignments";
        }
        else {
            $chr = $exon_chr;
        }
        
        my $exon_orient = $exon->{orient};
        if ($orient && $exon_orient ne $orient) {
            die "inconsistent exon orient";
        }
        else {
            $orient = $exon_orient;
        }

        $target = $exon->{target};
        $align_id = $exon->{align_id};
        
    }
    
    @chr_coords = sort {$a<=>$b} @chr_coords;
    @span_coords = sort {$a<=>$b} @span_coords;
    
    my $chr_coords_lend = $chr_coords[0];
    my $chr_coords_rend = $chr_coords[$#chr_coords];
    
    my $span_coords_lend = $span_coords[0];
    my $span_coords_rend = $span_coords[$#span_coords];
    

    my $per_id = $sum_per_id / $sum_len;

    my $span_struct = { lend => $chr_coords_lend,
                        rend => $chr_coords_rend,
                        orient => $orient,
                        chr => $chr,

                        target => $target,
                        align_id => $align_id,

                        range_lend => $span_coords_lend,
                        range_rend => $span_coords_rend,

                        exons => $exons_aref,
                    
                        per_id => sprintf("%.2f", $per_id),
                        
                    };

    return($span_struct);
}

1; #EOM
