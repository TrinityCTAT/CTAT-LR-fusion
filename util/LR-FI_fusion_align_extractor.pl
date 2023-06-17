#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 genePairContig.gtf LR_gmap.gff3 output_prefix\n\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $gff3_align_file = $ARGV[1] or die $usage;
my $output_prefix = $ARGV[2] or die $usage;

main: {
    
    my %orig_coord_info;
    my %scaffold_to_gene_coordsets = &parse_gtf_file($gtf_file, \%orig_coord_info);
    
    # organize original coordinate info
    my %scaffold_to_orig_coords = &organize_original_coordinate_info(\%orig_coord_info);
    
    #print STDERR Dumper(\%scaffold_to_gene_coordsets);
    
    my %scaffold_to_LR_coords = &parse_gff3_file($gff3_align_file);
    
    my %LR_fusion_trans_ids;

    foreach my $scaffold (keys %scaffold_to_gene_coordsets) {
        
        my @genes = keys %{$scaffold_to_gene_coordsets{$scaffold}};

        if (scalar @genes != 2) {
            die "Error, dont have only two genes for scaffold: $scaffold: " . Dumper(\@genes);
        }

        my ($geneA_coords_href, $geneB_coords_href) = &get_gene_coords($scaffold, $scaffold_to_gene_coordsets{$scaffold});
 
        #print "GeneA: " . Dumper($geneA_coords_href) 
        #    . "GeneB: " . Dumper($geneB_coords_href);
        
       
        my @trin_accs = keys (%{$scaffold_to_LR_coords{$scaffold}});
        foreach my $trin_acc (@trin_accs) {
            my $trin_coords_href = $scaffold_to_LR_coords{$scaffold}->{$trin_acc};
            
            # ignore singletons
            if (scalar (keys %$trin_coords_href) < 4) { next; } # at least 2 sets of coordinates, indicating an intron
            
            #print "LR: $trin_acc " . Dumper($trin_coords_href);

            if (&shared_coordinate($geneA_coords_href, $trin_coords_href)
                &&
                &shared_coordinate($geneB_coords_href, $trin_coords_href) ) {

                my ($break_left, $break_right) = &get_breakpoint_coords($geneA_coords_href, $geneB_coords_href, $trin_coords_href);

                $LR_fusion_trans_ids{$trin_acc} = "$scaffold:$break_left-$break_right";
            }

        }
    }

    &report_LR_fusions($gff3_align_file, \%LR_fusion_trans_ids, \%orig_coord_info, \%scaffold_to_orig_coords, $output_prefix);
    
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
sub parse_gtf_file {
    my ($gtf_file, $orig_coord_info_href) = @_;

    my %scaff_to_gene_to_coords;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
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
        
        my ($lend, $rend) = ($x[3], $x[4]);
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, [$lend, $rend]);
        


        # get original coordinate mapping info
        $info =~ /orig_coord_info \"([^,]+),(\d+),(\d+),([+-])\"/ or die "Error, cannot parse original coordinate info from $info";
        my $orig_chr = $1;
        my $orig_lend = $2;
        my $orig_rend = $3;
        my $orig_orient = $4;

        $orig_coord_info_href->{$scaffold_id}->{$lend} = { chrom => $orig_chr,
                                                           coord => $orig_lend,
                                                           orient => $orig_orient,
                                                           contig_coord => $lend,
        };

        $orig_coord_info_href->{$scaffold_id}->{$rend} = { chrom => $orig_chr,
                                                           coord => $orig_rend,
                                                           orient => $orig_orient,
                                                           contig_coord => $rend,
        };
        
        

    }
    close $fh;

    
    return(%scaff_to_gene_to_coords);
}




####
sub parse_gff3_file {
    my ($gff3_align_file) = @_;

    
    my %scaffold_to_trans_coords;
    
    open (my $fh, $gff3_align_file) or die $!;
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

        $scaffold_to_trans_coords{$scaff}->{$LR_id}->{$lend} = 1;
        $scaffold_to_trans_coords{$scaff}->{$LR_id}->{$rend} = 1;

    }
    close $fh;
    
    return(%scaffold_to_trans_coords);

}

####
sub report_LR_fusions {
    my ($gff3_align_file, $LR_ids_href, $orig_coord_info_href, $scaffold_to_orig_coords_href, $output_prefix) = @_;
    
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

        my ($left_genome_breakpoint, $left_ref_splice_mapping) = &infer_genome_breakpoint_from_local_coord($break_lend, $scaffold_orig_coord_info_href, $scaffold_coordinate_mappings_aref);
        my ($right_genome_breakpoint, $right_ref_splice_mapping)  = &infer_genome_breakpoint_from_local_coord($break_rend, $scaffold_orig_coord_info_href, $scaffold_coordinate_mappings_aref);

        my $splice_type = ($left_ref_splice_mapping == 1 && $right_ref_splice_mapping == 1) ? "ONLY_REF_SPLICE" : "INCL_NON_REF_SPLICE";
        
        push (@fusion_structs,
              { fusion_name => $scaffold,
                LeftLocalBreakpoint => $break_lend,
                RightLocalBreakpoint => $break_rend,
                LeftBreakpoint => $left_genome_breakpoint,
                RightBreakpoint => $right_genome_breakpoint,
                num_LR => $num_reads,
                LR_accessions => \@LR_reads,
                SpliceType => $splice_type,
              } );
    }

    @fusion_structs = reverse sort {$a->{num_LR} <=> $b->{num_LR}} @fusion_structs;

    
    print $LR_breakpoint_summary_ofh join("\t", "fusion_name", "num_LR", 
               "LeftLocalBreakpoint", "RightLocalBreakpoint",
               "LeftBreakpoint", "RightBreakpoint", "SpliceType",
               "LR_accessions") . "\n";

    foreach my $fusion (@fusion_structs) {

        print $LR_breakpoint_summary_ofh join("\t", 
                                              $fusion->{fusion_name}, 
                                              $fusion->{num_LR},
                                              $fusion->{LeftLocalBreakpoint},
                                              $fusion->{RightLocalBreakpoint},
                                              $fusion->{LeftBreakpoint},
                                              $fusion->{RightBreakpoint},
                                              $fusion->{SpliceType},
                                              join(",", @{$fusion->{LR_accessions}})) . "\n";
                
    }
    

    ## extract the chimeric alignments from the gff3 file.
    open (my $fh, $gff3_align_file) or die $!;
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
    my ($geneA_coords_href, $geneB_coords_href, $trin_coords_href) = @_;

    ## get left breakpoint
    my @left_shared_coords;

    foreach my $coord (keys %$geneA_coords_href) {
        if ($trin_coords_href->{$coord}) {
            push (@left_shared_coords, $coord);
        }
    }
    
    @left_shared_coords = sort {$a<=>$b} @left_shared_coords;
    my $left_breakpoint = pop @left_shared_coords;
    unless ($left_breakpoint) {
        confess "Error, no left breakpoint";
    }
    
    ## get right breakpoint
    my @right_shared_coords;

    foreach my $coord (keys %$geneB_coords_href) {
        if ($trin_coords_href->{$coord}) {
            push (@right_shared_coords, $coord);
        }
    }
    @right_shared_coords = sort {$a<=>$b} @right_shared_coords;

    my $right_breakpoint = shift @right_shared_coords;

    unless ($right_breakpoint) {
        confess "Error, no right breakpoint";
    }
       
    return($left_breakpoint, $right_breakpoint);
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
    my ($break_lend, $scaffold_orig_coord_info_href, $scaffold_coordinate_mappings_aref) = @_;

    
    if (my $struct = $scaffold_orig_coord_info_href->{$break_lend}) {
        my ($chrom, $coord, $orient) = ($struct->{chrom}, $struct->{coord}, $struct->{orient});
        
        return("$chrom:$coord:$orient", 1);
    }
    else {
        
        # map to closest coordinate.
        my $smallest_delta = undef;
        my $closest_struct = undef;
        foreach my $struct (@$scaffold_coordinate_mappings_aref) {
            my $contig_coord = $struct->{contig_coord};
            my $delta = $contig_coord - $break_lend;
            if ( (! defined($smallest_delta)) || abs($delta) < abs($smallest_delta)) {
                $smallest_delta = $delta;
                $closest_struct = $struct;
            }
        }

        my $chrom = $closest_struct->{chrom};
        my $orient = $closest_struct->{orient};
        my $coord;
        if ($orient eq '+') {
            $coord = $closest_struct->{coord} + $smallest_delta;
        }
        else {
            $coord = $closest_struct->{coord} - $smallest_delta;
        }
        return("$chrom:$coord:$orient", 0);
        
    }
    
}
