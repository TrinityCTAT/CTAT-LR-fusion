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
use SAM_reader;

my $usage = "\n\n\tusage: $0 trans.fasta gmap.map.gff3.chims_described left.fq right.fq\n\n";

my $trans_fasta = $ARGV[0] or die $usage;
my $chims_described = $ARGV[1] or die $usage;
my $left_fq_file = $ARGV[2] or die $usage;
my $right_fq_file = $ARGV[3] or die $usage;

my $ANCHOR = 12;
my $MIN_PERCENT_IDENTITY = 98;
my $MIN_ENTROPY = 1.5;

$trans_fasta = &ensure_full_path($trans_fasta);
$chims_described = &ensure_full_path($chims_described);
if ($left_fq_file ne "NA") {
    $left_fq_file = &ensure_full_path($left_fq_file);
}
if ($right_fq_file ne "NA") {
    $right_fq_file = &ensure_full_path($right_fq_file);
}

foreach my $file ($trans_fasta, $chims_described, $left_fq_file, $right_fq_file) {
    if ($file ne "NA") {
        unless (-s $file) {
            confess "Error, cannot locate file $file";
        }
    }
}


main: {
    
    my %chims = &parse_chims($chims_described);
    
    my $chim_candidates_fasta_filename = $chims_described . ".fasta";
    
    print STDERR "-extracting chim candidate seqs\n";
    my %chim_seqs = &extract_chim_candidate_seqs($trans_fasta, $chim_candidates_fasta_filename, \%chims);
    
    print STDERR "-getting breakpoint region entropy\n";
    my %chim_brkpt_entropy = &get_breakpoint_region_entropy(\%chims, \%chim_seqs);
    
    print STDERR "-computing trans_seq_entropy\n";
    my %trans_seq_entropy = &compute_trans_seq_entropy(\%chim_seqs);

    
    my %fusion_support;

    if ($left_fq_file ne "NA") {
        print STDERR "-aligning reads using bowtie2\n";
        my $bowtie2_bam = &align_reads_using_bowtie2($chim_candidates_fasta_filename, $left_fq_file, $right_fq_file);
    
        print STDERR "-capturing fusion support.\n";
        %fusion_support = &capture_fusion_support($bowtie2_bam, \%chims, \%trans_seq_entropy);

    }
    else {
        print STDERR "-no short reads provided. Skipping short read eval\n";
    }
    
        
    ## generate output, include junction and spanning frag support info:
    
    foreach my $target_trans_id (keys %chims) {
       
        my @junction_reads;
        my @spanning_frags;
        
        # init with just single long read support.
        my $J = 1;
        my $S = 0;
        
        if ($left_fq_file ne "NA") {  

            if (exists $fusion_support{$target_trans_id}) {
                my $info_href = $fusion_support{$target_trans_id};
                if (exists $info_href->{spanning}) {
                    @spanning_frags = keys %{$info_href->{spanning}};
                }
                if (exists $info_href->{junction}) {
                    @junction_reads = keys %{$info_href->{junction}};
                }
                
            }
            
            $J = scalar(@junction_reads);
            $S = scalar(@spanning_frags);
        }
        

        my $spanning_frag_list = join(",", @spanning_frags) || ".";
        my $junction_frag_list = join(",", @junction_reads) || ".";


        
        my $chim_info_aref = $chims{$target_trans_id};

        my $chim_brkpt_href = $chim_brkpt_entropy{$target_trans_id};
        
        
        foreach my $chim_info (@$chim_info_aref) {

            my $line = $chim_info->{line};


                    
            print join("\t", $line, $J, $S, 
                       
                       $chim_brkpt_href->{left_brkpt_anchor},
                       $chim_brkpt_href->{left_brkpt_entropy},
 
                       $chim_brkpt_href->{right_brkpt_anchor},
                       $chim_brkpt_href->{right_brkpt_entropy},
                       
                       $junction_frag_list, 
                       $spanning_frag_list) . "\n";
        }
    }
    
    exit(0);

}


####
sub capture_fusion_support {
    my ($bam_file, $chims_href, $seq_entropy_href) = @_;


    my %fusion_support;

    # sort according to read name and contig
    # these should all be perfect pairs

    my $pipeliner = new Pipeliner(-verbose => 2);

    my $sorted_alignments_file = "$bam_file.sorted_by_read_n_contig.sam.gz";
    my $cmd = "samtools view $bam_file | sort -T . -k1,1 -k3,3 | gzip -c - > $sorted_alignments_file";

    $pipeliner->add_commands(new Command($cmd, "$sorted_alignments_file.ok"));
    
    $pipeliner->run();

    
    my $sam_reader = new SAM_reader($sorted_alignments_file);

    while ($sam_reader->has_next()) {
        my $sam_entryA = $sam_reader->get_next();
        my $sam_entryB = $sam_reader->get_next();
        
        my $target_trans_id = $sam_entryA->get_scaffold_name();
        unless ($sam_entryB->get_scaffold_name() eq $target_trans_id) {
            confess "Error, trans IDs dont match for supposed pairs.";
        }

        my $frag_name = $sam_entryA->get_core_read_name();
        if ($frag_name ne $sam_entryB->get_core_read_name()) {
            confess "Error, core frag names dont match up: $frag_name vs. " . $sam_entryB->get_core_read_name();
        }
        

        unless ($sam_entryA->is_first_in_pair() xor $sam_entryB->is_first_in_pair()) {
            #print STDERR "error, have unpaired pair... skipping.\n";
            #print STDERR $sam_entryA->get_original_line() . "\n";
            #print STDERR $sam_entryB->get_original_line() . "\n";
            #print STDERR "\n";
           
            # note, there are some cases of multimapping pairs at the same gene (repetitive regions), ignore these too...

            next;
        }

        unless (&minimal_percent_identity($sam_entryA) && &minimal_percent_identity($sam_entryB)) {
            next; 
        }
        
        
        
        my $brkpt_range = $chims_href->{$target_trans_id}->[0]->{brkpt_range} or die "Error, no breakpoint range for transcript: $target_trans_id"; # brkpt is constant for all annotated entries of this transcript
        my ($break_left, $break_right) = split(/-/, $brkpt_range);


        my ($trans_coords_A_aref, @trash1) = $sam_entryA->get_alignment_coords();
        my ($trans_coords_B_aref, @trash2) = $sam_entryB->get_alignment_coords();

        # sort them according by coordinate
        ($trans_coords_A_aref, $trans_coords_B_aref) = sort {$a->[0]->[0] <=> $b->[0]->[0]} ($trans_coords_A_aref, $trans_coords_B_aref);

        my ($A_lend, $A_rend) = ($trans_coords_A_aref->[0]->[0], $trans_coords_A_aref->[$#$trans_coords_A_aref]->[1]);
        my ($B_lend, $B_rend) = ($trans_coords_B_aref->[0]->[0], $trans_coords_B_aref->[$#$trans_coords_B_aref]->[1]);
        
        if ($A_lend < $break_left && $B_rend > $break_right) {

            #################################
            ## fragment overlaps breakpoint.
            ##################################
            
            ## determine if it's a spanning pair or a fusion junction.
            if ($A_rend < $break_left && $B_lend > $break_right) {

                ####################
                ## a spanning pair:
                ####################
                
                my $seq_entropy_aref = $seq_entropy_href->{$target_trans_id} or die "Error, no seq entropy stored for transcript: $target_trans_id";

                #print STDERR "Found spanning fragment for $target_trans_id\n";
                if (&average_align_entropy($trans_coords_A_aref, $seq_entropy_aref) >= $MIN_ENTROPY
                    &&
                    &average_align_entropy($trans_coords_B_aref, $seq_entropy_aref) >= $MIN_ENTROPY) {
                    
                    $fusion_support{$target_trans_id}->{spanning}->{$frag_name}++;
                }
            }
            else {

                ##########################################################
                ## see if any alignment overlaps the point of the junction
                ##########################################################

                foreach my $align_seg (@$trans_coords_A_aref, @$trans_coords_B_aref) {
                    my ($lend, $rend) = @$align_seg;
                    
                    ## ensure the alignment meets the anchor requirement.

                    #                     brktp
                    #    ------------------|------------------------ 
                    #           <-- anchor on each side --->
                    
                    
                    if ($lend <= ($break_left - $ANCHOR) && $rend >= ($break_right + $ANCHOR)) {
                        # overlaps junction breakpoint
                        #print STDERR "Found a JUNCTION read for $target_trans_id\n";
                        $fusion_support{$target_trans_id}->{junction}->{$frag_name}++;
                        last;
                    }
                }
            }
        }
    }
    
    return(%fusion_support);

}


####
sub align_reads_using_bowtie2 {
    my ($trans_fasta, $left_fq_file, $right_fq_file) = @_;

    my $trans_dirname = dirname($trans_fasta);
    my $bowtie2_target = "$trans_dirname/bowtie2_target";
    
    my $cmd = "ln -sf $trans_fasta $bowtie2_target.fa";
    &process_cmd($cmd);

    my $pipeliner = new Pipeliner(-verbose=>2);
    $cmd = "bowtie2-build $bowtie2_target.fa $bowtie2_target > /dev/null";
    $pipeliner->add_commands(new Command($cmd, "$bowtie2_target.build.ok"));

    $cmd = "bash -c \"set pipefail -o && bowtie2 -k10 -p 4 --no-mixed --no-discordant --very-fast --end-to-end -x $bowtie2_target -1 $left_fq_file -2 $right_fq_file "
        . " | samtools view -F 4 -Sb - | samtools sort -@ 4 -o $trans_fasta.bowtie2.bam\"";
    $pipeliner->add_commands(new Command($cmd, "$trans_dirname/bowtie2_align.ok"));
    
    $pipeliner->run();
    
    return("$trans_fasta.bowtie2.bam");
    
}



####
sub extract_chim_candidate_seqs {
    my ($trans_fasta_filename, $output_filename, $chims_href) = @_;

    my %chim_seqs;

    my $fasta_reader = new Fasta_reader($trans_fasta_filename);
    open (my $ofh, ">$output_filename") or die "Error, cannot write to file: $trans_fasta_filename";
    
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $accession = $seq_obj->get_accession();
        
        if (exists $chims_href->{$accession}) {
            
            my $sequence = $seq_obj->get_sequence();
            
            print $ofh ">$accession\n$sequence\n";
            
            $chim_seqs{$accession} = $sequence;
            
        }
    }
    
    close $ofh;
        
            
    return (%chim_seqs);
    
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


####
sub get_breakpoint_region_entropy {
    my ($chims_href, $seqs_href) = @_;

    my %chim_brkpt_entropy;

    foreach my $acc (keys %$chims_href) {
        my $brkpt_range = $chims_href->{$acc}->[0]->{brkpt_range};
        my $sequence = $seqs_href->{$acc};

        my ($left_pt, $right_pt) = sort {$a<=>$b} split(/-/, $brkpt_range);
        
        my $left_seq_range = substr($sequence, $left_pt - $ANCHOR, $ANCHOR);
        my $right_seq_range = substr($sequence, $right_pt+1, $ANCHOR);

        my $left_brkpt_entropy = sprintf("%.2f", &compute_entropy($left_seq_range));
        my $right_brkpt_entropy = sprintf("%.2f", &compute_entropy($right_seq_range));

        $chim_brkpt_entropy{$acc}->{left_brkpt_entropy} = $left_brkpt_entropy;
        $chim_brkpt_entropy{$acc}->{right_brkpt_entropy} = $right_brkpt_entropy;
        
        $chim_brkpt_entropy{$acc}->{left_brkpt_anchor} = $left_seq_range;
        $chim_brkpt_entropy{$acc}->{right_brkpt_anchor} = $right_seq_range;
        
        
        
    }

    return (%chim_brkpt_entropy);
}

####
sub compute_entropy {
    my ($sequence) = @_;

    my @chars = split(//, $sequence);
    my %char_counter;

    foreach my $char (@chars) {
        $char_counter{$char}++;
    }
    
    my $num_chars = scalar(@chars);

    my $entropy = 0;
    foreach my $char (keys %char_counter) {
        my $count = $char_counter{$char};
        my $p = $count / $num_chars;

        $entropy += $p * (  log(1/$p) / log(2) );
    }

    return($entropy);
}

####
sub compute_trans_seq_entropy {
    my ($chim_seqs_href) = @_;

    my %trans_entropy;
    
    foreach my $trans_acc (keys %$chim_seqs_href) {
        my $sequence = $chim_seqs_href->{$trans_acc};
        
        for (my $i = 0; $i <= length($sequence) - $ANCHOR; $i++) {
            my $subseq = substr($sequence, $i, $ANCHOR);
            my $entropy = &compute_entropy($subseq);
            $trans_entropy{$trans_acc}->[$i] = $entropy;
        }
    }

    return(%trans_entropy);
}


####
sub minimal_percent_identity {
    my ($sam_entry) = @_;

    my $line = $sam_entry->get_original_line();

    my $seq_len = length($sam_entry->get_sequence());
    
    $line =~ /\s+NM:i:(\d+)/ or die "Error, cannot extract num mismatches from $line";
    my $num_mismatches = $1;

    my $percent_identity = 100 - ($num_mismatches/$seq_len * 100);

    if ($percent_identity > $MIN_PERCENT_IDENTITY) {
        return(1);
    }
    else {
        return(0);
    }

}
                               
            
####
sub average_align_entropy {
    my ($align_coords_aref, $seq_entropy_aref) = @_;

    my @entropies;
    foreach my $align_seg (@$align_coords_aref) {
        my ($lend, $rend) = @$align_seg;
        for (my $i = $lend; $i <= $rend - $ANCHOR; $i++) {
            push (@entropies, $seq_entropy_aref->[$i]);
        }
    }

    my $sum = 0;
    foreach my $entropy (@entropies) {
        $sum += $entropy;
    }

    my $avg_entropy = $sum/scalar(@entropies);

    return($avg_entropy);
}


