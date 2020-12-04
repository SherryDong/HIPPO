#!/usr/bin/perl -w

@gene = split ",",$ARGV[0];
$platform = $ARGV[1];
$result=$ARGV[2];
if($platform eq "hg19"){
	$ref_fa = "/cluster/apps/refseq/ANNOVAR/humandb_hg19/hg19_v0_Homo_sapiens_assembly19.fasta";
}
##
foreach $gene (@gene){
	$cmd = "perl src/gene2bed.pl $gene $platform 0 0 >$result/$gene.bed";
	system $cmd;
	$cmd = "seqtk subseq $ref_fa $result/$gene.bed";
	$txt = `$cmd | tail -n 1`; chomp($txt);
	$strand = `cut -f 6 $result/$gene.bed | tail -n 1`;chomp($strand);
	if($strand eq "-"){
		$txt =  reverse($txt);
		$txt =~ s/A/t/g;
		$txt =~ s/T/a/g;
		$txt =~ s/G/c/g;
		$txt =~ s/C/g/g;
		$txt = uc($txt);
	}
	print ">$gene\n$txt\n";
	open O,">$result/$gene.fa";
	print O ">$gene\n$txt\n";
	close O;
}


