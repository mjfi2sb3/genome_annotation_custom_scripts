#!/bin/bash -x
# this piepline is written according to the steps outlined at http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced
#

root_dir=$PWD;
genome=$root_dir"/S.Pimp.5kb.v1.0.fasta";
pref="s.pimp";

# 1. MITEs (miniature inverted repeat transposable elements) 
module load perl muscle/priv.3.8.31 mdust ncbi-blast/2.2.26 mite-hunter;
mkdir -p mite;
cd mite;
perl `which MITE_Hunter_manager.pl` -i $genome -g $pref -n 40 -S 12345678
cat *Step8_*.fa > MITE.lib
cd ..;

# 2. LTR (long terminal repeat) retrotransposons 
# 2.1. Collection of relatively recent LTR retrotranposons
mkdir -p ltr;
cd ltr;
module purge;
module load applications-extra genometools perl;
# 2.1.1. Collection of candidate elements with LTRs that are 99% or more in similarity using LTRharvest
gt suffixerator -db ${genome} -indexname seqfileindex -tis -suf -lcp -des -ssp -dna
gt ltrharvest -index seqfileindex -out seqfile.out99 -outinner seqfile.outinner99 -gff3 seqfile.gff99 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10  > seqfile.result99;

# 2.1.2. Using LTRdigest to find elements with PPT (poly purine tract) or PBS (primer binding site)
gt gff3 -sort seqfile.gff99 > seqfile.gff99.sort;
gt ltrdigest -trnas /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/gtrnadb/eukaryota/eukaryotic-tRNAs.fa seqfile.gff99.sort seqfileindex > seqfile.gff99.dgt;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step1.pl --gff seqfile.gff99.dgt;

#2.1.3. Further filtering of the candidate elements
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile seqfile.out99 --resultfile seqfile.result99 --sequencefile $genome --removed_repeats CRL_Step2_Passed_Elements.fasta

mkdir fasta_files;
mv  Repeat_*.fasta fasta_files/.;
mv  CRL_Step2_Passed_Elements.fasta  fasta_files/.;
cd  fasta_files;
module load muscle;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step3.pl --directory $PWD --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25;
mv  CRL_Step3_Passed_Elements.fasta  ../.;
cd  ..;

#2.1.4. Identify elements with nested insertions
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/ltr_library.pl --resultfile seqfile.result99  --step3 CRL_Step3_Passed_Elements.fasta --sequencefile ${genome};
cat lLTR_Only.lib ../mite/MITE.lib > repeats_to_mask_LTR99.fasta
module purge;
module load private perl/5.20.1  ncbi-rmblast gcc hmmer trf ncbi-rmblast/2.2.28 repeatmasker/priv.4.0.7;
RepeatMasker -pa 40 -lib repeats_to_mask_LTR99.fasta   -nolow  -dir .  seqfile.outinner99;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/cleanRM.pl  seqfile.outinner99.out  seqfile.outinner99.masked  >  seqfile.outinner99.unmasked
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/rmshortinner.pl  seqfile.outinner99.unmasked  50 > seqfile.outinner99.clean;
module purge;
module load ncbi-blast+;
makeblastdb -in /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017DNA  -dbtype prot;
blastx  -query  seqfile.outinner99.clean -db /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017DNA -num_threads 40 -evalue 1e-10 -num_descriptions 10 -out seqfile.outinner99.clean_blastx.out.txt
perl  /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/outinner_blastx_parse.pl --blastx  seqfile.outinner99.clean_blastx.out.txt  --outinner  seqfile.outinner99;

# 2.1.5 Building examplars
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step4.pl -step3 CRL_Step3_Passed_Elements.fasta --resultfile seqfile.result99 --innerfile passed_outinner_sequence.fasta  --sequencefile ${genome};
makeblastdb -in lLTRs_Seq_For_BLAST.fasta  -dbtype nucl;
blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads 40 -out lLTRs_Seq_For_BLAST.fasta.out;
makeblastdb -in  Inner_Seq_For_BLAST.fasta  -dbtype nucl
blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta  -evalue 1e-10 -num_descriptions 1000 -out  Inner_Seq_For_BLAST.fasta.out;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out --step3 CRL_Step3_Passed_Elements.fasta --final LTR99.lib --pcoverage 90 --pidentity 80;

mkdir 99 && mv *99* lLTR* Inner_Seq_For_BLAST.fasta* passed_outinner_sequence.fasta CRL_Step* fasta_files 99/. && ln 99/LTR99.lib;


#######################################################
# 2.2.	Collection of relatively old LTR retrotransposons
module purge;
module load applications-extra genometools perl;
gt suffixerator -db ${genome} -indexname seqfileindex -tis -suf -lcp -des -ssp -dna
gt ltrharvest -index seqfileindex -out seqfile.out85 -outinner seqfile.outinner85 -gff3 seqfile.gff85 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10  > seqfile.result85;
gt gff3 -sort seqfile.gff85 > seqfile.gff85.sort;
gt ltrdigest -trnas /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/gtrnadb/eukaryota/eukaryotic-tRNAs.fa seqfile.gff85.sort seqfileindex > seqfile.gff85.dgt;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step1.pl --gff seqfile.gff85.dgt;

perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile seqfile.out85 --resultfile seqfile.result85 --sequencefile $genome --removed_repeats CRL_Step2_Passed_Elements.fasta

mkdir fasta_files;
mv  Repeat_*.fasta fasta_files/.;
mv  CRL_Step2_Passed_Elements.fasta  fasta_files/.;
cd  fasta_files;
module load muscle;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step3.pl --directory $PWD --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25;
mv  CRL_Step3_Passed_Elements.fasta  ../.;
cd  ..;

2.1.4. Identify elements with nested insertions
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/ltr_library.pl --resultfile seqfile.result85  --step3 CRL_Step3_Passed_Elements.fasta --sequencefile ${genome};
cat lLTR_Only.lib ../mite/MITE.lib > repeats_to_mask_LTR85.fasta
module purge;
module load private perl/5.20.1  ncbi-rmblast gcc hmmer trf ncbi-rmblast/2.2.28 repeatmasker/priv.4.0.7;
RepeatMasker -pa 40 -lib repeats_to_mask_LTR85.fasta   -nolow  -dir .  seqfile.outinner85;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/cleanRM.pl  seqfile.outinner85.out  seqfile.outinner85.masked  >  seqfile.outinner85.unmasked
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/rmshortinner.pl  seqfile.outinner85.unmasked  50 > seqfile.outinner85.clean;
module purge;
module load ncbi-blast+;
makeblastdb -in /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017DNA  -dbtype prot;
blastx  -query  seqfile.outinner85.clean -db /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017DNA -num_threads 40 -evalue 1e-10 -num_descriptions 10 -out seqfile.outinner85.clean_blastx.out.txt
perl  /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/outinner_blastx_parse.pl --blastx  seqfile.outinner85.clean_blastx.out.txt  --outinner  seqfile.outinner85;

2.1.5 Building examplars
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step4.pl -step3 CRL_Step3_Passed_Elements.fasta --resultfile seqfile.result85 --innerfile passed_outinner_sequence.fasta  --sequencefile ${genome};
makeblastdb -in lLTRs_Seq_For_BLAST.fasta  -dbtype nucl;
blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads 40 -out lLTRs_Seq_For_BLAST.fasta.out;
makeblastdb -in  Inner_Seq_For_BLAST.fasta  -dbtype nucl
blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta  -evalue 1e-10 -num_descriptions 1000 -out  Inner_Seq_For_BLAST.fasta.out;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out --step3 CRL_Step3_Passed_Elements.fasta --final LTR85.lib --pcoverage 90 --pidentity 80;

mkdir 85 && mv *85* lLTR* Inner_Seq_For_BLAST.fasta* passed_outinner_sequence.fasta CRL_Step* fasta_files 85/. && ln 85/LTR85.lib;

#######################################################




module purge;
module load private perl/5.20.1  ncbi-rmblast gcc hmmer trf ncbi-rmblast/2.2.28 repeatmasker/priv.4.0.7;
awk -F "," '{print $1}' LTR85.lib > LTR85.fix.lib # shorten ids for repeatmask
RepeatMasker -pa 40 -lib LTR99.lib -dir . LTR85.lib;
mv LTR85.fix.lib.masked LTR85.lib.masked
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/remove_masked_sequence.pl  --masked_elements  LTR85.lib.masked  --outfile  FinalLTR85.lib;
cat LTR99.lib FinalLTR85.lib > allLTR.lib;

3. Collecting repetitive sequences by RepeatModeler
cat allLTR.lib ../mite/MITE.lib > allMITE_LTR.lib;
# S.Pimp.v1.0.RMviridplante.fasta was repeat masked prior using repeatmasker and Viridiplantae as --species
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/rmaskedpart.pl  /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/S.Pimp.v1.0.RMviridplante.fasta  50  >  S.Pimp.v1.0.RMviridplante.fasta.um;
RepeatMasker -pa 45 -lib allMITE_LTR.lib -dir . S.Pimp.v1.0.RMviridplante.fasta.um;
ln -s S.Pimp.v1.0.RMviridplante.fasta.um.masked seqfile.masked;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/rmaskedpart.pl  seqfile.masked  50  >  umseqfile;
module purge;
module load private perl/5.20.1  ncbi-rmblast gcc hmmer trf ncbi-rmblast/2.2.28 repeatmasker/priv.4.0.7 repeatscout recon nseg repeatmodeler/priv.1.0.9;
BuildDatabase -name umseqfiledb -engine ncbi umseqfile;
RepeatModeler -pa 40 -database umseqfiledb >& umseqfile.out;
# rm=`ls -1 -dtr RM_*|tail -1` && cd $rm;
# rm_output_dir=$PWD;
# cd ..;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/repeatmodeler_parse.pl --fastafile umseqfiledb-families.fa --unknowns repeatmodeler_unknowns.fasta --identities repeatmodeler_identities.fasta;
module purge; module load perl ncbi-blast+;
makeblastdb -in /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017  -dbtype prot;
blastx -query repeatmodeler_unknowns.fasta -db /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/protein_TE/Tpases031017 -num_threads 40 -evalue 1e-10 -num_descriptions 10 -out modelerunknown_blast_results.txt;
perl /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/CRL_Scripts1.0/transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknowns.fasta;
mv  unknown_elements.txt  ModelerUnknown.lib;
cat  identified_elements.txt  repeatmodeler_identities.fasta  > ModelerID.lib;

# 4. Exclusion of gene fragments
cat ModelerID.lib allLTR.lib ../mite/MITE.lib > KnownRepeats.lib
cat KnownRepeats.lib ModelerUnknown.lib > allRepeats.lib;

blastx -query allRepeats.lib -num_threads 40 -db /local/data/kaust_research/tomatoes/s.pimp/local.copy/rm/ltr/plant.refseq.prexp201710.faa  -evalue 1e-10 -num_descriptions 10 -out allRepeats.lib_blast_results.txt
module load gcc/4.7.4 hmmer private protexcluder;
perl /local/data/apps/ProtExcluder1.2/ProtExcluder.pl -option  allRepeats.lib_blast_results.txt  allRepeats.lib