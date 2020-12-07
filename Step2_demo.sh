## prepare the following files for Step2 test
# demo/SBDS_ref/combine_mod.fa
# demo/SBDS_ref/combine_mod.info
# demo/SBDS_ref/combine_mod_region.bed
# demo/SBDS_ref/SBDS.bed
# demo/SBDS_ref/SBDSP1.bed
# demo/demo.bam
# demo/demo.bam.bai

########################
true_gene=SBDS
pseudo_gene=SBDSP1
total_gene_num=2
result='demo/' ## main directory to save results
combine_ref_dir='demo/SBDS_ref/' ## main directory for combined reference

# Step2: generate reads-region mapping content matrix for a sample
bam_file='demo/demo.bam'
sample_name='demo'

region_true=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$true_gene.bed | tail -n 1`
region_pseudo=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$pseudo_gene.bed | tail -n 1`
# for each sample
samtools view -h $bam_file $region_true $region_pseudo >$result/$sample_name.sam
samtools fasta $result/$sample_name.sam >$result/$sample_name.fa
blat $combine_ref_dir/combine_mod.fa $result/$sample_name.fa $result/$sample_name.output.psl
awk -F " " {'print "chr\t"$16"\t"$17"\t"$10"\t.\t"$9"\t"$12"\t"$13'} $result/$sample_name.output.psl | sed '1,5d' >$result/$sample_name.output.bed
bedtools intersect -a $result/$sample_name.output.bed -b $combine_ref_dir/combine_mod_region.bed -wo >$result/$sample_name.overlap
perl src/seq2matrix.pl $result/$sample_name.fa $result/$sample_name.overlap $result/$sample_name.mat
