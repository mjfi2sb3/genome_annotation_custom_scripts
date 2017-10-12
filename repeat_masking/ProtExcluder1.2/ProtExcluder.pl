#! /usr/bin/perl -w

$usage = "Protexcluder.pl -f bpsofflankingtoremove(default=50)  blastx/blastpfile fastafile \n";

# to exclude the portion matching protein subject in a nucleotide sequence file

if (@ARGV < 2) {die "$usage";}

use Getopt::Std;

getopts("f:");

$Len  = defined $opt_f ? $opt_f : 50;

`rm -f $ARGV[1].ssi`;

`/local/data/apps/ProtExcluder1.2/matchtract.pl $ARGV[0] > $ARGV[0].mt`;

`/local/data/apps/ProtExcluder1.2/countaanu.pl $ARGV[0].mt > $ARGV[0].mtca`;

`/local/data/apps/ProtExcluder1.2/rmlowcomplexitymathc.pl  $ARGV[0].mtca 3 60 >  $ARGV[0].mtca_3_60`;

`/local/data/apps/ProtExcluder1.2/blastformatProt.pl $ARGV[0] > $ARGV[0].f`;

`/local/data/apps/ProtExcluder1.2/rmlowcomfromBF.pl $ARGV[0].mtca_3_60 $ARGV[0].f > $ARGV[0].fnolow`;

`sort -k 6,6 -k 3,3n $ARGV[0].fnolow > $ARGV[0].fnolows`;

`/local/data/apps/ProtExcluder1.2/mergequeryBF.pl $ARGV[0].fnolows $Len > $ARGV[0].fnolowm50`;

`/local/data/apps/ProtExcluder1.2/unmatchedregionBF.pl $ARGV[0].fnolowm50 $Len > $ARGV[0].fnolowm50MSP`;

`/local/data/apps/ProtExcluder1.2/mspesl-sfetch.pl $ARGV[1] $ARGV[0].fnolowm50MSP 0 $ARGV[0].fnolowm50seq`;

`/local/data/apps/ProtExcluder1.2/mergeunmatchedregion.pl $ARGV[0].fnolowm50seq > $ARGV[0].fnolowm50seqm`;

`/local/data/apps/ProtExcluder1.2/GCcontent.pl $ARGV[0].fnolowm50seqm > $ARGV[0].fnolowm50seqmGC`;

`/local/data/apps/ProtExcluder1.2/rmshortseq_noN.pl $ARGV[0].fnolowm50seqmGC $ARGV[0].fnolowm50seqm 50 > $ARGV[0].fnolowm50seqmns`;

`/local/data/apps/ProtExcluder1.2/getanycolumnuni.pl $ARGV[0].fnolow 6 > $ARGV[0].fnolowlist`;

`/local/data/apps/ProtExcluder1.2/rmlistedseq.pl $ARGV[0].fnolowlist $ARGV[1] >  $ARGV[1]nPr`;

`cat $ARGV[1]nPr $ARGV[0].fnolowm50seqmns > temp`;

`/local/data/apps/ProtExcluder1.2/fasta-reformat.pl temp 50 > $ARGV[1]noProtFinal`;
