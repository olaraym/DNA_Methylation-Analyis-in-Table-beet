
# Methylation analysis of the Pacbio_seq data Using the CCS_meth pipeline.
# The CCS_meth pipeline used in this study is based on a deep learning approach and was designed to detect methylation 
# in the 5mCpG context from circular consensus sequences (Hifi reads), it uses an attention-based bidirectional Gated Recurrent Unit (GRU) 
# networks to predict DNA methylation states of cytosine in Hifi reads (Peng N. et al., 2022). The approach also consisted of a series of 
# steps which is detailed in Figure 7B. and in the following lines of codes below;

# Modification probability calling from the aligned PacBio reads using the code below;

ccsmeth call_mods --input hifi.pbmm2.bam --ref ../TB-assembly.posscaf.fa \
--model_file model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v1.ckpt \
--output ./new/cmod --threads 50 --threads_call 10 --model_type attbigru2s \
--rm_per_readsite --mode align --modbam modbam --motif CG

# Extraction of modification probabilities for CG methylation context using the line of code below;

ccsmeth call_freqt --input_path /path/to/cmod.per_redsite.tsv \ 
	--result_file /path/to/new_CG_ccsmeth --bed –sort\ 
	--prob_cf 0.0 --refsites_only \
	--motifs CG --mod_loc 0 --ref /path/to/TB-assembly.posscaf.fa \
	--contigs --threads 30

# Next, the command line below using awk was used for extracting and sorting the necessary columns 
# and information, see the line of code below;
awk '{OFS="\t"} {print $1"\t"$2"\t"$3"\t"$6"\t"$11"\t"$10}' new_CG_ccsmeth | awk '$4=="+"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}$4==\
"-"{print $1"\t"$2-1"\t"$3-1"\t"$4"\t"$5"\t"$6}' | sort -u -k1,1 -k2,2n | awk '$6>=4' > mCG_ccs.bed


