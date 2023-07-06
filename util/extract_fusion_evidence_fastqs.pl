#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fastq_reader;
use Process_cmd;
use DelimParser;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);  


my $usage = <<__EOUSAGE__;

######################################################################################
#
#  --fusions <string>           fusion predictions 'final' output file (not abridged)
#
#  Reads:
#
#  --fastq <string>             reads in fastq file #
#
#  --output_prefix <string>     output prefix
#
# * optional:
#
#  --by_fusion <string>         dirname; write per fusion ev files.
#
#  --by_isoform <string>        dirname;  define fusions according to breakpoint instead of more generically by gene pair
#                               Provide the name of the directory to store the fusion files.  The fq read files will be
#                               named according to the fusion^isoform_brkpt name.
#
######################################################################################


__EOUSAGE__

    ;


my $help_flag;

my $fusion_results_file;
my $fastq;
my $output_prefix;
my $BY_ISOFORM;
my $BY_FUSION;

&GetOptions( 'help|h' => \$help_flag,
             
             'fusions=s' => \$fusion_results_file,
             
             'fastq=s' => \$fastq,

             'output_prefix=s' => \$output_prefix,
             
             'by_fusion=s' => \$BY_FUSION,

             'by_isoform=s' => \$BY_ISOFORM,
             
             
    );

if ($help_flag) {
    die $usage;
}

unless ($fusion_results_file && $fastq && $output_prefix) {
    die $usage;
}


if ($BY_ISOFORM) {
    
    if (-d $BY_ISOFORM) {
        die "Error, directory $BY_ISOFORM already exists.  Choose a new name or remove the current directory before rerunning.";
    }
    my $ret = system("mkdir -p $BY_ISOFORM");
    if ($ret) {
        die "Error, couldn't create directory $BY_ISOFORM ";
    }
}

if ($BY_FUSION) {
    if (-d $BY_FUSION) {
        die "Error, directory $BY_FUSION already exists.  Choose a new name or remove the current directory before rerunning.";
    }
    my $ret = system("mkdir -p $BY_FUSION");
    if ($ret) {
        die "Error, couldn't create directory $BY_FUSION ";
    }
}    


main: {

    ## get the core fragment names:
    my %LR_read_name_to_fusion_name;
    
    open (my $fh, $fusion_results_file) or die "Error, cannot open file $fusion_results_file";
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'} or die "Error, cannot get fusion name from " . Dumper($row);
        
        my $fusion_name_simple = $fusion_name;

        if ($BY_ISOFORM) {
            my $left_brkpt = $tab_reader->get_row_val($row, "LeftBreakpoint");
            my $right_brkpt = $tab_reader->get_row_val($row, "RightBreakpoint");
            my $brkpt_info = "${left_brkpt}-${right_brkpt}";
            $fusion_name .= "^$brkpt_info";
        }
        
        my $junction_reads_list_txt = $tab_reader->get_row_val($row, "LR_accessions");
        
        foreach my $junction_read (split(/,/, $junction_reads_list_txt)) {
            
            $LR_read_name_to_fusion_name{$junction_read} = $fusion_name;
        }
    }
        
    &write_fastq_files($fastq, $output_prefix, \%LR_read_name_to_fusion_name);
    
    print STDERR "\nDone.\n\n";
    
    exit(0);
    
}


####
sub write_fastq_files {
    my ($fastq_file, $output_prefix, $LR_read_name_to_fusion_name_href) = @_;

    my $output_fastq_file = "$output_prefix.fusion_evidence_reads.fq";

    my @samples = split(/,/, $sample_names);
    
    open (my $ofh, ">$output_fastq_file") or die "Error, cannot write to $output_fastq_file";

    my %reads_to_capture = keys %{$LR_read_name_to_fusion_name_href};
    
    print STDERR "-searching fq file: $fastq_file\n";

    my $fastq_reader = new Fastq_reader($fastq_file);
        
    while (my $fq_record = $fastq_reader->next()) {
            
        my $core_read_name = $fq_record->get_read_name();
        #print STDERR "[$core_read_name]\n";
            
        my $record_text = $fq_record->get_fastq_record();
        chomp $record_text;
        
        if (my $simple_fusion_name = $LR_read_name_to_fusion_name_href->{re$core_frag_to_simple_fusion_href->{$core_read_name}) {
            
            

        if ($BY_FUSION) {
        
                open(my $fusion_reads_ofh, ">>$BY_FUSION/${simple_fusion_name}${output_fastq_file_suffix}.fq") or die "Error, cannot append to file $BY_FUSION/${simple_fusion_name}${output_fastq_file_suffix}.fq";
                print $fusion_reads_ofh $record_text . "\n";
                close $fusion_reads_ofh;
            }
        }
        
        my $fusion_instances_with_read_aref = $core_frag_name_to_fusion_name_href->{$core_read_name};
        
        if (ref $fusion_instances_with_read_aref) {
            
            my @fusion_instances = @$fusion_instances_with_read_aref;
            
            
            my @lines = split(/\n/, $record_text);
            
            my $reported_in_full_file = 0;
            
            foreach my $fusion_instance (@fusion_instances) {    
                
                my ($_1, $_2, $_3, $_4) = @lines;
                
                $_3 = "+$fusion_instance"; # encode the fusion name in the 3rd line, which is otherwise useless anyway
                
                unless ($reported_in_full_file) {
                    # report in the full file only once.
                    print $ofh join("\n", ($_1, $_2, $_3, $_4)) . "\n";
                    $reported_in_full_file = 1;
                }
                
                if ($BY_ISOFORM) {
                    my ($fusion_instance_filename, $rest) = split(/\|[JS][12]?\|/, $fusion_instance);
                    
                    $fusion_instance_filename =~ s/\W+/_/g; 
                    $fusion_instance_filename =~ s/_$//;
                    
                    open(my $isoform_ofh, ">>$BY_ISOFORM/${fusion_instance_filename}${output_fastq_file_suffix}.fq") or die "Error, cannot append to file $BY_ISOFORM/${fusion_instance}${output_fastq_file_suffix}.fq";
                    print $isoform_ofh join("\n", ($_1, $_2, $_3, $_4)) . "\n";
                    close $isoform_ofh;
                }
            }
            delete $reads_to_capture{$core_read_name} if exists $reads_to_capture{$core_read_name};
        }
    }
    
    

    print STDERR "\nDone writing to $output_fastq_file\n\n";
    
    if (%reads_to_capture) {
        confess "Error, failed to capture fusion evidence reads: " . Dumper(\%reads_to_capture);
    }
    
    return;
}


####
sub append_reads_to_fusion {
    my ($fusion_name, $core_frag_name_to_fusion_name_href, $reads_href) = @_;

    foreach my $frag_name (keys %$reads_href) {
        $core_frag_name_to_fusion_name_href->{$frag_name} = $fusion_name ;
    }
    
    return;
}

####
sub parse_core_frag_names {
    my ($comma_delim_read_list_txt) = @_;

    my %core_frag_names;

    my @read_names = split(/,/, $comma_delim_read_list_txt);
    foreach my $read_name (@read_names) {
        $read_name =~ s|/[12]$||;
        $core_frag_names{$read_name} = 1;
    }

    
    return(%core_frag_names);
}


####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @sample_names;
    my @left_fqs;
    my @right_fqs;
    
    open(my $fh, $samples_file) or die "Error, cannot open file: $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($sample_name, $left_fq, $right_fq) = @x;

        push (@sample_names, $sample_name);
        
        push (@left_fqs, $left_fq);
        
        if ($right_fq) {
            push (@right_fqs, $right_fq);
        }
    }

    my $sample_names_txt = join(",", @sample_names);
    my $left_fq_files = join(",", @left_fqs);
    my $right_fq_files = join(",", @right_fqs);


    return($sample_names_txt, $left_fq_files, $right_fq_files);
}
