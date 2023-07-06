#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fastq_reader;


my $usage = <<__EOUSAGE__;


Extracts fusion evidence reads based on read annotations in 3rd line of fq file

Writes each fusion to a file named according to the fusion type.


#####################################################
#
# --fastq <string>     :  fastq file
#
# --output_dir <string>    /path/to/output/directory
#
#######################################################


__EOUSAGE__

    ;


my $help_flag;
my $fastq;
my $output_dir;
   

&GetOptions ( 'h' => \$help_flag,
              'fastq=s' => \$fastq,
              'output_dir=s' => \$output_dir,
    );


unless ($fastq && $output_dir) {
    die $usage;
}

main: {
    
    unless (-d $output_dir) {
        mkdir($output_dir) or die "Error, cannot mkdir $output_dir";
    }

    my %fusion_to_reads = &capture_reads_by_fusion($fastq);
    &write_fusion_reads(\%fusion_to_reads, $output_dir);

    print STDERR "\n\nDone.\n\n";
    
    exit(0);

}

####
sub capture_reads_by_fusion {
    my ($fq_file) = @_;

    my %fusion_to_reads;
    
    my $fq_reader = new Fastq_reader($fq_file);

    while(my $fq_record = $fq_reader->next()) {

        my $fq_record_text = $fq_record->get_fastq_record();

        my @lines = split(/\n/, $fq_record_text);
        my $fusion_info = $lines[2];

        $fusion_info =~ s/^\+//;
        my @pts = split(/\^/, $fusion_info);
        my $fusion_name = $pts[0];

        push (@{$fusion_to_reads{$fusion_name}}, $fq_record_text);
    }

    return(%fusion_to_reads);
    
}


####
sub write_fusion_reads {
    my ($fusion_to_reads_href, $output_dir) = @_;
    
    my @fusions = keys %$fusion_to_reads_href;
    foreach my $fusion (@fusions) {
        my @fq_records = @{$fusion_to_reads_href->{$fusion}};

        my $fq_filename = "$output_dir/${fusion}.fastq";
        print STDERR "-writing $fq_filename\n";
        open(my $ofh, ">$fq_filename") or die "Error, cannot write to file: $fq_filename";
        print $ofh join("", @fq_records);
        close $ofh;
    }

    return;
}
