#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;
use List::Util qw(min max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);




my $usage = <<__EOUSAGE__;

###########################################################
#
#  --FI_gtf <string>        :  FI contigs gtf filename
#
#  --LR_gff3 <string>       :  LR alignments in gff3 format.
#
#  --output_prefix <string> :  prefix for output files
#
#  --snap_dist <int>        :  if breakpoint is at most this distance from a reference exon boundary, position gets snapped to the splice site.
#
###########################################################


__EOUSAGE__

    ;


my $help_flag;
my $FI_gtf_filename;
my $LR_gff3_filename;
my $output_prefix;
my $SNAP_dist;



&GetOptions ( 'help|h' => \$help_flag,
              'FI_gtf=s' => \$FI_gtf_filename,
              'LR_gff3=s' => \$LR_gff3_filename,
              'output_prefix=s' => \$output_prefix,
              'snap_dist=i' => \$SNAP_dist,
    );

if ($help_flag) {
    die $usage;
}

unless ($FI_gtf_filename && $LR_gff3_filename && $output_prefix && defined($SNAP_dist) ) {
    die $usage;
}



my $DONOR_TYPE = "DONOR";
my $ACCEPTOR_TYPE = "ACCEPTOR";
my $NA_TYPE = "NA";

main: {
    
    my %orig_coord_info;
    my %scaffold_to_gene_coordsets = &parse_FI_gtf_filename($FI_gtf_filename, \%orig_coord_info);
    
    # organize original coordinate info
    my %scaffold_to_orig_coords = &organize_original_coordinate_info(\%orig_coord_info);
    
    #print STDERR Dumper(\%scaffold_to_gene_coordsets);
    
    my %scaffold_to_LR_coords = &parse_LR_alignment_gff3_file($LR_gff3_filename);
    
    my %LR_fusion_trans_ids;

    foreach my $scaffold (keys %scaffold_to_gene_coordsets) {
        
        my @genes = keys %{$scaffold_to_gene_coordsets{$scaffold}};

        if (scalar @genes != 2) {
            die "Error, dont have only two genes for scaffold: $scaffold: " . Dumper(\@genes);
        }

        my ($geneA_coords_href, $geneB_coords_href) = &get_gene_coords($scaffold, $scaffold_to_gene_coordsets{$scaffold});
 
        my $geneA_max = max(keys %$geneA_coords_href);
        my $geneB_min = min(keys %$geneB_coords_href);
        
        unless (exists $scaffold_to_LR_coords{$scaffold}) { next; }

        my @LR_accs = keys %{$scaffold_to_LR_coords{$scaffold}};
        foreach my $LR_acc (@LR_accs) {
            my @LR_coordsets = sort {$a->[0]<=>$b->[0]} @{$scaffold_to_LR_coords{$scaffold}->{$LR_acc}};
            
            # ignore singletons
            if (scalar(@LR_coordsets) < 2) { next; } # at least 2 sets of coordinates, indicating an intron
            
            my $min_LR_coord = $LR_coordsets[0]->[0];
            my $max_LR_coord = $LR_coordsets[$#LR_coordsets]->[1];
            

            if ($min_LR_coord < $geneA_max && $max_LR_coord > $geneB_min) { # spans both genes 
                
                my ($break_left, $break_right) = &get_breakpoint_coords(\@LR_coordsets, $geneA_max, $geneB_min);
                
                $LR_fusion_trans_ids{$LR_acc} = "$scaffold:$break_left-$break_right";
            }

        }
    }

    &report_LR_fusions($LR_gff3_filename, \%LR_fusion_trans_ids, \%orig_coord_info, \%scaffold_to_orig_coords, $output_prefix);
    
    exit(0);
}


####
sub shared_coordinate {
    my ($coordsA_href, $coordsB_href) = @_;

    foreach my $coord (keys %$coordsA_href) {
        
        if ($coordsB_href->{$coord}) {
            return(1);
        }
    }

    return(0);
}



####
sub get_gene_coords {
    my ($scaffold, $genes_to_coords_href) = @_;

    my @gene_coords_hrefs;
    foreach my $gene (split(/--/, $scaffold)) {
        
        my @coords = @{$genes_to_coords_href->{$gene}};

        my %gene_coords;
        foreach my $coordpair (@coords) {
            my ($lend, $rend) = @$coordpair;
            $gene_coords{$lend} = 1;
            $gene_coords{$rend} = 1;
        }
        push (@gene_coords_hrefs, \%gene_coords);
    }

    return(@gene_coords_hrefs);

}


####
sub parse_FI_gtf_filename {
    my ($FI_gtf_filename, $orig_coord_info_href) = @_;

    my %scaff_to_gene_to_coords;

    my %scaff_to_trans_to_coords;

    open (my $fh, $FI_gtf_filename) or die "Error, cannot open file $FI_gtf_filename";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $scaffold_id = $x[0];
        my $type = $x[2];
        
        unless ($type eq 'exon') { next; }
        
        my $info = $x[8];
        my $gene_id = "";
        
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_id = $1;
        }
        elsif ($info =~ /FI_gene_label \"([^\"]+)\"/) {
            $gene_id = $1;
            my @x = split(/\^/, $gene_id);
            $gene_id = $x[0];
        }
        else {
            die "Error, not able to extract gene_name or FI_gene_label value from $info";
        }
        

        my $transcript_id;
        if ($info =~ /transcript_id \"([^\"]+)\"/) {
            $transcript_id = $1;
        }
        else {
            die "Error, cannot extract transcript id from $info";
        }
        
        

        my ($lend, $rend) = ($x[3], $x[4]);
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, [$lend, $rend]);
        
        push (@{$scaff_to_trans_to_coords{$scaffold_id}->{$transcript_id}}, [$lend, $rend]);
        
        
        my $strand = $x[6];
        unless ($strand eq '+') {
            confess "Error, FI contigs should always have annotations on the + strand only";
        }
        

        # get original coordinate mapping info
        $info =~ /orig_coord_info \"([^,]+),(\d+),(\d+),([+-])\"/ or die "Error, cannot parse original coordinate info from $info";
        my $orig_chr = $1;
        my $orig_lend = $2;
        my $orig_rend = $3;
        my $orig_orient = $4;
        
        my ($orig_end5, $orig_end3) = ($orig_orient eq '+') ? ($orig_lend, $orig_rend) : ($orig_rend, $orig_lend);
        
        $orig_coord_info_href->{$scaffold_id}->{$lend} = { chrom => $orig_chr,
                                                           coord => $orig_end5,
                                                           orient => $orig_orient,
                                                           contig_coord => $lend,
                                                           splice_junc => $NA_TYPE, # update later
        };

        $orig_coord_info_href->{$scaffold_id}->{$rend} = { chrom => $orig_chr,
                                                           coord => $orig_end3,
                                                           orient => $orig_orient,
                                                           contig_coord => $rend,
                                                           splice_junc => $NA_TYPE, # update later
        };
        
        

    }
    close $fh;

    
    ###############################
    ## Update splice junction info:
    
    foreach my $scaffold_id (keys %scaff_to_trans_to_coords) {
        
        foreach my $transcript_id (keys %{$scaff_to_trans_to_coords{$scaffold_id}}) {
            
            my @coordsets = sort {$a->[0]<=>$b->[0]} @{$scaff_to_trans_to_coords{$scaffold_id}->{$transcript_id}};
            
            # assign splice juncs:
            for (my $i = 0; $i < $#coordsets; $i++) {
                my $left_junc = $coordsets[$i]->[1];
                my $right_junc = $coordsets[$i+1]->[0];
                
                $orig_coord_info_href->{$scaffold_id}->{$left_junc}->{splice_junc} = $DONOR_TYPE;
                $orig_coord_info_href->{$scaffold_id}->{$right_junc}->{splice_junc} = $ACCEPTOR_TYPE;
            }
        }
    }
    
    return(%scaff_to_gene_to_coords);
}




####
sub parse_LR_alignment_gff3_file {
    my ($LR_gff3_filename) = @_;

    
    my %scaffold_to_trans_coords;
    
    open (my $fh, $LR_gff3_filename) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; } # comment line
        unless (/\w/) { next; }
        
        chomp;
        my @x = split(/\t/);
        my $scaff = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $info = $x[8];

        my $LR_id;
        if ($info =~ /ID=([^;]+)\.p\d+;/) {

            $LR_id = $1;
        }
        else {
            die "Error, cannot find LR ID from $info";
        }

        push (@{$scaffold_to_trans_coords{$scaff}->{$LR_id}}, [$lend, $rend]);
        

    }
    close $fh;
    
    return(%scaffold_to_trans_coords);

}

####
sub report_LR_fusions {
    my ($LR_gff3_filename, $LR_ids_href, $orig_coord_info_href, $scaffold_to_orig_coords_href, $output_prefix) = @_;
    
    my $chimeric_trans_gff3_filename = "$output_prefix.gff3";
    open(my $chimeric_trans_gff3_ofh, ">$chimeric_trans_gff3_filename") or die "Error, cannot write to file: $chimeric_trans_gff3_filename";
    
    my $LR_fusion_breakpoint_summary_filename = "$output_prefix.breakpoint_info.tsv";
    open(my $LR_breakpoint_summary_ofh, ">$LR_fusion_breakpoint_summary_filename") or die "Error, cannot write to $LR_fusion_breakpoint_summary_filename";
    
    my %scaff_breakpoint_to_read_support;
    foreach my $LR_id (keys %$LR_ids_href) {
        my $scaff_breakpoint = $LR_ids_href->{$LR_id};
        print "#LRFusionTranscript:\t$LR_id\t$scaff_breakpoint\n";
        push (@{$scaff_breakpoint_to_read_support{$scaff_breakpoint}}, $LR_id);
    }
        
    ## generate fusion breakpoint summary report
    my @fusion_structs;
    foreach my $breakpoint (keys %scaff_breakpoint_to_read_support) {
        my @LR_reads = @{$scaff_breakpoint_to_read_support{$breakpoint}};
        my $num_reads = scalar(@LR_reads);
        my ($scaffold, $breakpoint_coords) = split(/:/, $breakpoint);
        my ($break_lend, $break_rend) = split(/-/, $breakpoint_coords);
        
        my $scaffold_orig_coord_info_href = $orig_coord_info_href->{$scaffold};
        my $scaffold_coordinate_mappings_aref = $scaffold_to_orig_coords_href->{$scaffold};

        my ($adj_break_lend, $left_genome_breakpoint, $left_ref_splice_mapping) = &infer_genome_breakpoint_from_local_coord($break_lend, 
                                                                                                                            $scaffold_orig_coord_info_href, 
                                                                                                                            $scaffold_coordinate_mappings_aref,
                                                                                                                            $DONOR_TYPE,
                                                                                                                            $SNAP_dist);
        
        # in case it snapped:
        $break_lend = $adj_break_lend;

        my ($adj_break_rend, $right_genome_breakpoint, $right_ref_splice_mapping)  = &infer_genome_breakpoint_from_local_coord($break_rend, 
                                                                                                                               $scaffold_orig_coord_info_href, 
                                                                                                                               $scaffold_coordinate_mappings_aref,
                                                                                                                               $ACCEPTOR_TYPE,
                                                                                                                               $SNAP_dist);
        
        # in case it snapped:
        $break_rend = $adj_break_rend;
        
        my $splice_type = ($left_ref_splice_mapping == 1 && $right_ref_splice_mapping == 1) ? "ONLY_REF_SPLICE" : "INCL_NON_REF_SPLICE";
        
        my ($left_gene, $right_gene) = split(/--/, $scaffold);

        push (@fusion_structs,
              { fusion_name => $scaffold,
                LeftGene => $left_gene,
                RightGene => $right_gene,
                LeftLocalBreakpoint => $break_lend,
                RightLocalBreakpoint => $break_rend,
                LeftBreakpoint => $left_genome_breakpoint,
                RightBreakpoint => $right_genome_breakpoint,
                num_LR => $num_reads,
                LR_accessions => \@LR_reads,
                SpliceType => $splice_type,
              } );
    }


    @fusion_structs = &merge_identical_breakpoints(@fusion_structs);
    

    @fusion_structs = reverse sort {$a->{num_LR} <=> $b->{num_LR}} @fusion_structs;

    
    

    
    print $LR_breakpoint_summary_ofh join("\t", "#FusionName", "num_LR", 
                                          "LeftGene",
                                          "LeftLocalBreakpoint", 
                                          "LeftBreakpoint", 
                                          "RightGene",
                                          "RightLocalBreakpoint",
                                          "RightBreakpoint", 
                                          "SpliceType",
                                          "LR_accessions") . "\n";

    foreach my $fusion (@fusion_structs) {

        print $LR_breakpoint_summary_ofh join("\t", 
                                              $fusion->{fusion_name}, 
                                              $fusion->{num_LR},
                                              $fusion->{LeftGene},
                                              $fusion->{LeftLocalBreakpoint},
                                              $fusion->{LeftBreakpoint},
                                              $fusion->{RightGene},
                                              $fusion->{RightLocalBreakpoint},
                                              $fusion->{RightBreakpoint},
                                              $fusion->{SpliceType},
                                              join(",", @{$fusion->{LR_accessions}})) . "\n";
                
    }
    

    ## extract the chimeric alignments from the gff3 file.
    open (my $fh, $LR_gff3_filename) or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my $line = $_;
        my @x = split(/\t/);
        my $info = $x[8];
        
        if ($info =~ /ID=([^;]+)\.p\d+;/) {

            my $LR_id = $1;
            
            if ($LR_ids_href->{$LR_id}) {
                print $chimeric_trans_gff3_ofh $line;
            }
        }
        else {
            die "Error, cannot find LR ID from $info";
        }
        
    }
    close $fh;


    close($LR_breakpoint_summary_ofh);
    close($chimeric_trans_gff3_ofh);
    

    return;
}

####
sub get_breakpoint_coords {
    my ($LR_coordsets_aref, $geneA_max, $geneB_min) = @_;
    
    for (my $i = 0; $i < $#$LR_coordsets_aref; $i++) {
        
        my $segment_left_aref = $LR_coordsets_aref->[$i];
        my $segment_right_aref = $LR_coordsets_aref->[$i+1];

        # check adjacent exon coordinate boundaries to see if they are closest to the different gene boundaries
        # as expected for a fusion breakpoint.

        my $left_end = $segment_left_aref->[1];
        my $right_end = $segment_right_aref->[0];
        
        if (&is_closer($left_end, $geneA_max, $geneB_min) && &is_closer($right_end, $geneB_min, $geneA_max)) {
            return($left_end, $right_end);
        }
    }

    confess ("Error, not finding a proper fusion breakpoint for : " . Dumper($LR_coordsets_aref) . " with gene bounds {$geneA_max, $geneB_min} ");
    
}

####
sub is_closer {
    # determine if coord_A is closer to coord_B than to coord_C
    my ($coord_A, $coord_B, $coord_C) = @_;

    if (abs($coord_A - $coord_B) < abs($coord_A - $coord_C)) {
        return(1);
    }
    else {
        return(0);
    }
}


####
sub organize_original_coordinate_info { 
    my ($orig_coord_info_href) = @_;
    
    my %scaffold_to_orig_coords;

    foreach my $scaffold (keys %$orig_coord_info_href) {
        my @coordinates = keys %{$orig_coord_info_href->{$scaffold}};
        
        my @coord_structs;
        foreach my $coordinate (@coordinates) {
            my $struct = $orig_coord_info_href->{$scaffold}->{$coordinate};
            
            push (@coord_structs, $struct);

        }
        @coord_structs = sort {$a->{contig_coord}<=>$b->{contig_coord}} @coord_structs;

        $scaffold_to_orig_coords{$scaffold} = \@coord_structs;

    }

    return(%scaffold_to_orig_coords);
}



sub infer_genome_breakpoint_from_local_coord {
    my ($break_coord, $scaffold_orig_coord_info_href, $scaffold_coordinate_mappings_aref, $SPLICE_SITE_TYPE, $SNAP_dist) = @_;

    my $struct = $scaffold_orig_coord_info_href->{$break_coord};
    if (defined($struct) && $struct->{splice_junc} eq $SPLICE_SITE_TYPE) {
        my ($chrom, $coord, $orient) = ($struct->{chrom}, $struct->{coord}, $struct->{orient});
        
        #print STDERR "-no snap, found ref splice match for $break_coord\n";
        return($break_coord, "$chrom:$coord:$orient", 1);
    }
    else {
        
        # map to closest coordinate.
        my $smallest_delta = undef;
        my $closest_struct = undef;
        foreach my $struct (@$scaffold_coordinate_mappings_aref) {
            my $contig_coord = $struct->{contig_coord};
            my $delta = $break_coord - $contig_coord;
            
            if ( (! defined($smallest_delta)) 
                 || 
                 ($closest_struct->{splice_junc} ne $SPLICE_SITE_TYPE && $struct->{splice_junc} eq $SPLICE_SITE_TYPE) # always prefer matching splice junc type
                 ||
                 ($closest_struct->{splice_junc} eq $struct->{splice_junc} && abs($delta) < abs($smallest_delta) ) # same splice type but closer.
                ) {
                
                $smallest_delta = $delta;
                $closest_struct = $struct;
            }
        }

        my $chrom = $closest_struct->{chrom};
        my $orient = $closest_struct->{orient};
        my $coord;
        if (abs($smallest_delta) <= $SNAP_dist && $closest_struct->{splice_junc} eq $SPLICE_SITE_TYPE) {
            $coord = $closest_struct->{coord};
            my $adj_break = $closest_struct->{contig_coord};
            
            #print STDERR "Snapping $break_coord -> $adj_break\n";
            
            return($adj_break, "$chrom:$coord:$orient", 1); # now at splice site
        }
        else {
            if ($orient eq '+') {
                $coord = $closest_struct->{coord} + $smallest_delta;
            }
            else {
                $coord = $closest_struct->{coord} - $smallest_delta;
            }
            #print STDERR "No snap, here\'s delta info: $chrom:$coord:$orient\tbreak_coord: $break_coord\tsmallest_delta: $smallest_delta\tcontig_info: " . Dumper($closest_struct);

            return($break_coord, "$chrom:$coord:$orient", 0);
        }
        
        
        
    }
    
}


####
sub merge_identical_breakpoints {
    my (@fusion_structs) = @_;

    my %fusion_token_to_consolidated_fusions;

    foreach my $fusion (@fusion_structs) {
        my $token = join("$;", 
                         $fusion->{fusion_name},
                         $fusion->{LeftLocalBreakpoint},
                         $fusion->{RightLocalBreakpoint},
                         $fusion->{LeftBreakpoint},
                         $fusion->{RightBreakpoint},
                         $fusion->{SpliceType});

        if (my $existing_fusion_struct = $fusion_token_to_consolidated_fusions{$token}) {

            #print STDERR "-merging fusion structs for " . Dumper($fusion) . " with " . Dumper($existing_fusion_struct);
            
            # add info to existing fusion:
            $existing_fusion_struct->{num_LR} += $fusion->{num_LR};
            push (@{$existing_fusion_struct->{LR_accessions}}, @{$fusion->{LR_accessions}});
        }
        else {
            $fusion_token_to_consolidated_fusions{$token} = $fusion;
        }
    }
    
    return (values %fusion_token_to_consolidated_fusions);
}

