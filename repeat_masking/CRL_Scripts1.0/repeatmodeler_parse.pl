#!/usr/bin/perl -w
# This script will parse the output of RepeatModeler into the unknown elements and those with identification.
# Megan Bowman
# 09 April 2014

use strict;
use Getopt::Long();
use Bio::Seq;
use Bio::SeqIO; 

my $usage = "\nUsage: $0 --fastafile <fasta file name> --unknowns <name of unknowns output file> --identities <name of output file with identities>\n";

my ($fastafile, $in, $out, $unknowns, $identities);

Getopt::Long::GetOptions('fastafile=s' => \$fastafile,
			 'unknowns=s' => \$unknowns,
			 'identities=s' => \$identities);     

if (!-e $fastafile) {
    die "$fastafile does not exist!\n";
} 

if (-e $unknowns) {
    die "$unknowns already exists!\n";
}

if (-e $identities) {
    die "$identities already exists!\n";
}

my ($id, $desc, $seq, $name, $out1, $out2, $seqobj);

$in  = Bio::SeqIO->new(-file => "$fastafile", -format => 'Fasta');
$out1 = Bio::SeqIO->new(-file => ">$unknowns", -format => 'Fasta');
$out2 = Bio::SeqIO->new(-file => ">$identities", -format => 'Fasta');

while ($seqobj = $in ->next_seq()) {
  $id = $seqobj->display_id();
  if ($id =~ /rnd-\d+_family-\d+#.+\/.+$/) {
    $desc = $seqobj->desc();
    $seq = $seqobj->seq();
    $out2->write_seq($seqobj);
    next;
  } 
  if ($id =~ /rnd-\d+_family-\d+#Unknown$/) {
    $desc = $seqobj->desc();
    $seq = $seqobj->seq(); 
    $out1->write_seq($seqobj);
    next;
  } 
  if ($id =~ /rnd-\d+_family-\d+#.+$/) {
    $desc = $seqobj->desc();
    $seq = $seqobj->seq();
    $out2->write_seq($seqobj);
    next;
  }   
  else {
    print "$id\n";
  }
}


#rnd-1_family-392#LTR/Gypsy
#rnd-1_family-8#LTR/Gypsy
#rnd-1_family-707#DNA/CMC-EnSpm


