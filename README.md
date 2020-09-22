# Custom codes for the PromethION oral phageome paper 
This page explains custom codes used in the paper "Long-read metagenomics using PromethION uncovers novel oral bacteriophages, and interaction with host bacteria " by Koji Yahara , Masato Suzuki, Aki Hirabayashi, Wataru Suda, Masahira Hattori, Yutaka Suzuki, Yusuke Okazaki 

A figure summarizing the bioinformatic workflow is available [here](http://yahara.hustle.ne.jp/linked_figures/FigS1_workflow_v1.0.pdf) and as a supplementary figure of the paper.

Most steps in the workflow use only publicly available software, whereas the following steps use the custom codes explained below
- examination of novelty of viral contigs (mostly accordingly to Paez-Espino et al (2017) Nature Protocol)
- analysis of genomic context of prophages
- analysis of remote homologs
- CRISPR spacer analysis

This page explains the custome codes according to the guildeline "For manuscripts utilizing custome algorithms or software that are central to the research and no yet described in published literature" of Nature Research.

This page is written based on "Code and Software Submission Checklist" of Nature Research.

## System requirements 
- Linux OS
- UGE/SGE for submitting parallel jobs (otherwise, please comment out all "qsub" in the code to execute each command one-by-one)
- Perl (v5.16)
- mcl (v14-137)
- Prokka (v1.13)

  


#### Any required non-standard hardware

None for the custom codes.

In the whole workflow, assembly of the large metagenome data (>30Gb per sample) requires > 1TB RAM.


## Installation guide

#### Instruction

- "git clone"  this repository and go to the directory  
- Please download [a 12GB compressed file](https://www.dropbox.com/s/hmnj61yx2mep0zb/PromethION_phageome_demo_github.tgz?dl=0) containing input (and ouptut) files for running the demo below, and uncompress it
- Please install the softwares in the "System requirements" above


#### Typical install time on a "normal" computer

Within one hour, mostly for downloading the 12GB compressed file.



## Demo

The codes are written below.   Expected output and run time are explained at each step.  The output files are in advance stored in the "output_example" folder included in the 12.5GB downloadable compressed file.

## Instruction for use

If you would like to run the custom code on your data, please modify it using your text editor and execute it so that your data file is used as input of the code.


------
##### Examination of novelty of viral contigs (genome clustering with IMG/VR 2.0 database) 

According to step 8-11 in Paez-Espino et al (2017) Nature Protocol,

```
DB_IMG=IMGVR_all_nucleotides.fna # downloaded from IMG/VR

for aa in `find . -maxdepth 2 -type d -name assembly_split_fas_cat_pilonCorrected_virsorter`
do
  #
  # 8|  Run the following Blast command:
  # 
  qsub -cwd -S /bin/bash  -l s_vmem=16G,mem_req=16G <<< "
  blastn \
    -query ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta \
    -db ${DB_IMG} \
    -outfmt '6 std qlen slen' \
    -out   ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.blout \
    -evalue 1.0e-50 \
    -perc_identity 80 \
    -num_threads 4
  "
done

# after the submitted jobs finish,

for aa in ` find . -maxdepth 2 -type d -name assembly_split_fas_cat_pilonCorrected_virsorter `
do
  #
  # 9|  Remove self-hits from the Blast output by running the following command:
  #
  cat              ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.blout |
  awk '$1 != $2' > ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.blout

  #
  # 10| Parse the output using specific cutoffs 
  #     (similarity ≥90%; covered length ≥75%; covered length requires at least one contig of >1,000-bp length) by running the following command:
  #
  cd script1/
  module  load java/1.8.0_45
  java Parse_BLAST ../${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.blout \
                 > ../${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.sim90cov75hitLeast1kb.blout
  # For simplicity, after compilation, users are advised to run the executable directly from the folder in which  
  # the .class files are located.

  cd .. 
  
  #
  # 11| Carry out single linkage clustering by running the following command:
  #
  mcl               ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.sim90cov75hitLeast1kb.blout --abc \
                 -o ${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.sim90cov75hitLeast1kb.blout.slc
done


```
(expected run time: 10 minutes)

------
##### Examination of novelty of viral contigs (together with results obtained by vConTACT v2) 

```
for aa in ` find . -maxdepth 2 -type d -name assembly_split_fas_cat_pilonCorrected_virsorter `
do
  # asssuming that output directories of vContact2 are located at ${aa}.vcontact2
  if [ ! -d "${aa}.vcontact2" ]; then
    echo "Error: ${aa}.vcontact2 doesn't exist"
    exit 1
  fi

  PREFIX="${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat"
  
  perl ./utils/check_bp_fasta_length_desc.pl -f ${PREFIX}.fasta -a > ${PREFIX}.fasta.length.txt
  
  CMD="perl ./utils/examine_VIRSorter_viral_contigs.v2.pl "
  CMD="${CMD} -p ${aa}/VIRSorter_global-phage-signal.csv"
  CMD="${CMD} -l ${PREFIX}.fasta.length.txt "
  CMD="${CMD} -c ${PREFIX}.fasta_vs_mVCs.noSelfHits.sim90cov75hitLeast1kb.blout.slc "
  
  echo ${CMD}
  eval ${CMD} 
done

# summarize for viral categories 1,2,4,5 
for aa in ` find . -maxdepth 2 -type d -name assembly_split_fas_cat_pilonCorrected_virsorter `
do    aa=${aa}/Predicted_viral_sequences/VIRSorter_cat123456.cat.fasta_vs_mVCs.noSelfHits.sim90cov75hitLeast1kb.blout.slc.info.sort
  perl ./utils/summarize_examined_VIRSorter_viral_contigs.v3.pl -f ${aa}
done > summarize_examined_VIRSorter_viral_contigs.v3.pl.out.txt

```
(expected output: genomeS2018_flye___\*/\*_virsorter/Predicted_viral_sequences/*.blout.slc.info.sort for each sample, used for Table S4.  "cluster_vs_mVCs_1column" and "refFamily" correspond to "clustered with IMG/VR v2.0" and "assinged family", respectively.)
(expected output: summarize_examined_VIRSorter_viral_contigs.v3.pl.out.txt, which corresponds to Table 2 and Figure S4)
(expected run time: within a few minutes)


------
##### Analysis of genomic context of prophages (for annotation of surrounding sequences)

```
# --------------------------------------------
# For surrounding sequences of each 
# “most confidence” and “likely” prophage
# --------------------------------------------
for aa in `ls genomeS2018_flye_*/assembly_split_fas_cat_pilonCorrected_virsorter/Predicted_viral_sequences/VIRSorter_prophages_cat-[45].gb`
do
  bb=`echo ${aa} | perl -pe 's/_virsorter\/Predicted_viral_sequences\/VIRSorter_prophages_cat-[45].gb/.fas.length/g'`
  
  CMD="perl ./utils/extract_VIRSorter_prophage_context.pl -f ${aa} -l ${bb}"
  echo ${CMD}
  eval ${CMD}
done

# for the specific contig of sample 3 (Figure 2)
VIR_ID=VIRSorter_contig_181_pilon_gene_21_gene_120-18237-95638-cat_4
CATEGORY=`echo ${VIR_ID} | perl -pe 's/^.*cat_//g'`
CONTIG_ID=`echo ${VIR_ID} | perl -pe 's/^.*(contig_[0-9]+).*$/$1/g'`

for aa in `ls ./genomeS2018_flye_MS1809_5/assembly_split_fas_cat_pilonCorrected_virsorter/Predicted_viral_sequences/VIRSorter_prophages_cat-${CATEGORY}.gb.context.host*.fas | grep -v __contig_` # left or right
do
  CMD="perl ./utils/extract_fas_by_key.pl -f ${aa} -p _${CONTIG_ID}_ "
  echo ${CMD}
  HOST_FAS=`${CMD}`
  
  CMD="perl ./utils/execute_prokka.pl -f ${HOST_FAS} -a '--locustag $ID --force'"
  echo ${CMD}
  eval ${CMD}
  
  OUTDIR=`echo ${HOST_FAS} | perl -pe 's/\.fas/_prokka/g'`
  CONTIG_ID=`echo ${HOST_FAS} | perl -pe 's/^.*(contig_[0-9]+)_.*$/$1/g'`
  
  GFF=`ls ${OUTDIR}/*.gff`
  perl -i -pe "s/contig1/${CONTIG_ID}/g" ${GFF}
  echo ${GFF}
done
```
(expected output: VIRSorter_prophages_cat-4.gb.context.hostRight__contig_181_.gff, VIRSorter_prophages_cat-4.gb.context.hostLeft_contig_181_.gff, used for Figure 2A)
(expected run time: within a few minutes)

------
##### Analysis of genomic context of prophages (distribution of difference in length between a prophage assembled from the long-reads and a corresponding viral sequence in IMG/VR v2.0 database)

```
VIRAL_LISTF=cat45_clustered_with_IMGVR.tsv 
# category 4 & 5 clusterd with a viral seqeunce in IMG/VR 
# extracted from Table S4 (Viral sequences identified by VirSorter)

perl check_cat45_clustered_with_IMGVR.pl \
  -f ${VIRAL_LISTF} \
  -m DB/IMGVR_all_Sequence_information.vOUT_N4_hostDefined23906.txt \
  -l DB/IMGVR_all_nucleotides.fna.length > cat45_clustered_with_IMGVR.Host.longestLenDiff.tsv
```
(expected output: cat45_clustered_with_IMGVR.Host.longestLenDiff.tsv, used for Figure 2B)
(expected run time: within a  minutes)

------
##### Analysis of remote homologs

```
# assuming 
#   phate_phanotate.gff produced by PHANOTATE 
#   *.hhblitsUniClust30.gz files produced by hhblits (gzipped) 
#   are in PipelineOutput and PipelineOutput_jumbo directories

for DIR in `find PipelineOutput*/ -type d -name \*out`
do
  echo ${DIR}
  qsub -cwd -S /bin/bash -l s_vmem=16G,mem_req=16G <<< "
    perl ./utils/hhblits_summary_after_phanotate.pl -d ${DIR} > ${DIR}.hhblits.gff
  "
done

grep "[[:cntrl:]]resistance[[:cntrl:]]" PipelineOutput/*.hhblits.gff  > remote_homolog_AMRgenes_percentIden.highlyAbundant.tsv

grep "[[:cntrl:]]lysogeny[[:cntrl:]]"   PipelineOutput/*.hhblits.gff  > remote_homolog_integrase_percentIden.highlyAbundant.tsv

for aa in `cat highlyAbundant28Siphoviridae.list`
do
  grep ${aa} remote_homolog_AMRgenes_percentIden.highlyAbundant.tsv
done | perl -pe 's/gff:VIRSorter/gff\tVIRSorter/g' > remote_homolog_AMRgenes_percentIden.highlyAbundant28Siphoviridae.tsv 

grep "[[:cntrl:]]resistance[[:cntrl:]]" PipelineOutput_jumbo/*.hhblits.gff | 
grep -v VIRSorter_contig_811_pilon-cat_2 | 
grep -v VIRSorter_contig_659_pilon-cat_2 > remote_homolog_AMRgenes_percentIden.jumbo.tsv
```
(expected output: remote_homolog_*.tsv in which 11th and 12rd columns indicate %identity and %aligned length explained in the main text
(expected run time: several minutes)

------
##### CRISPR spacer analysis

```
#
# assuming CAT and CRISPRdetected were conducted in advance
#

#
# combine results of CAT and VIRSorter
#
for FAS_FOR_BT in `ls genomeS2018_flye_*/assembly_split_fas_cat_pilonCorrected.fas`
do 
  BT_PREFIX=`echo ${FAS_FOR_BT} | perl -pe 's/\.f[a-z]+$//g'`
  echo ${BT_PREFIX}
  
  CONTIG_LENGTH=${FAS_FOR_BT}.length
  if [ ! -s "${CONTIG_LENGTH}" ]; then
    echo "Error ${CONTIG_LENGTH}"
    break
  fi
  
  VIRSORTER_DIR=${BT_PREFIX}_virsorter
  if [ ! -d "${VIRSORTER_DIR}" ]; then
    echo "Error ${VIRSORTER_DIR}"
    break
  fi

  CAT_NAMES=${BT_PREFIX}_CAT.contig2classification_names.txt
  if [ ! -s "${CAT_NAMES}" ]; then
    echo "Error ${CAT_NAMES}"
    break
  fi

  CMD="perl ./utils/format_CAT_with_VIRSorter.pl -s ${VIRSORTER_DIR} -c ${CAT_NAMES}" 
  echo ${CMD}
  eval ${CMD}
done
# output: *_names.classified.txt



#
# BLAST each CRISPR spacer against the viral sequences detected by VIRSorter
# (20 minutes)
#
for aa in `ls genomeS2018_flye_*/assembly_split_fas_cat_pilonCorrected_CAT.contig2classification_names.classified.txt`
do
  EACH_SAMPLE_FLYE_DIR=`echo ${aa} | perl -pe 's/\/assembly_split_fas_cat_pilonCorrected_CAT\.contig2classification_names\.classified\.txt//g'`
  
  GFF_CRISPR=${EACH_SAMPLE_FLYE_DIR}_CRISPRdetect.all.cat.gff
  VIRSORTER_GBK_DIR=${EACH_SAMPLE_FLYE_DIR}/assembly_split_fas_cat_pilonCorrected_virsorter/Predicted_viral_sequences/
  VIRSORTER_CAT_FASTA=${VIRSORTER_GBK_DIR}/VIRSorter_cat123456.cat.fasta
  
  OUTFILE=`echo ${GFF_CRISPR}.CAT.classified.txt | perl -pe 's/^.*\///g'`
  
  CMD="perl ./utils/compare_CRISPRdetect_CAT_classified.pl -c ${aa} -p ${GFF_CRISPR} -g ${VIRSORTER_GBK_DIR} -f ${VIRSORTER_CAT_FASTA}  > ${OUTFILE}"
  echo ${CMD}
  qsub -cwd -S /bin/bash -l medium -l s_vmem=32G,mem_req=32G <<< "${CMD}"
done
# output: *_CRISPRdetect.all.cat.gff.CAT.classified.txt

#
# After the submitted jobs finish, summarize
#
perl ./utils/summarize_CRISPRdetect.all.cat.gff.CAT.classified.txt.pl -p genomeS2018_flye_ -t 1245    > summarize_CRISPRdetect.all.cat.gff.CAT.classified.txt.category1245.tsv



#
# BLAST each CRISPR spacer against IMG/VR v2.0 database
# (6 hours)
#
DB_IMG=IMGVR_all_nucleotides.fna # downloaded from IMG/VR

for aa in `ls genomeS2018_flye_*_CRISPRdetect.all.cat.gff.CAT.classified.txt`
do
  OUTFILE=`echo ${aa} | perl -pe 's/^.*\///g' | perl -pe 's/\.classified\.txt/.hit.txt/g' | perl -pe 's/\.CAT\.hit\.txt/.IMGVR.hit.txt/g'` 
  
  CMD="perl ./utils/compare_CRISPRdetect_IMGVR_v2.pl -g ${aa}  -f ${DB_IMG} > ${OUTFILE}"
  echo ${CMD}
  qsub -cwd -S /bin/bash -l medium -l s_vmem=32G,mem_req=32G <<< "${CMD}"
done
# output: *_CRISPRdetect.all.cat.gff.IMGVR.hit.txt

#
# After the submitted jobs finish, summarize
#
arr=(
  bacterial
  viral_cat12
  prophage_cat45
)

echo "sample spacer_bacterial spacer_viral_cat12 spacer_prophage_cat45 protospacer_bacterial protospacer_viral_cat12 protospacer_prophage_cat45" > summarize_CRISPRdetect.all.cat.gff.IMGVR.hit.txt.category1245.tsv

for aa in `ls genomeS2018_flye_*_CRISPRdetect.all.cat.gff.IMGVR.hit.txt`
do
  OUT_LINE=${aa}

  # num. spacers
  for TYPE2 in ${arr[@]}
  do
    WC=`grep ${TYPE2} ${aa} | wc -l`
    OUT_LINE="${OUT_LINE} ${WC}"
  done

  # num. hit protospacers
  for TYPE2 in ${arr[@]}
  do
    WC=`grep ${TYPE2} ${aa} | grep Query | wc -l`
    OUT_LINE="${OUT_LINE} ${WC}"
  done

  echo ${OUT_LINE}
done >> summarize_CRISPRdetect.all.cat.gff.IMGVR.hit.txt.category1245.tsv
```

(expected output: summarize_CRISPRdetect.*.txt.category1245.tsv used for Figure 5)
(expected run time: within a hour)
