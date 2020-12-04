#!/usr/bin/perl -w

$fa = $ARGV[0];
$overlap = $ARGV[1];
###############
open F,$fa;
while(<F>){
	chomp;
	if(/^>(.*)/){
		$name = $1;
		next;
	}
	$seq = $_;
	$seq{$name} = $seq;
	$rev_seq = reverse($seq);
	$rev_seq =~ tr/ATCG/TAGC/;
	$rev_seq{$name} = $rev_seq;
}
close F;
open O,$overlap or die $!;
# 0-based
# chr 7642  7792  ST-E00244:576:HTF53CCXY:4:2106:22242:25323/1  . - 0 150 chr 7636  7643  OID_7636-7643 1
# chr 7642  7792  ST-E00244:576:HTF53CCXY:4:2106:22242:25323/1  . - 0 150 chr 7635  7643  OID_7635-7643 1
while(<O>){
	chomp;
	@a = split "\t";
	$name = $a[3];
	$strand = $a[5];
	$sa = $a[1];$ea = $a[2];$la = $ea-$sa;
	$sb = $a[9];$eb = $a[10];
	$sr = $a[6];$er = $a[7];$lr = $er-$sr;
	if($la > $lr){next;} ## do not consider split reads
	######################	
	$pre = ""; $suf = "";
	$s = $sb - $sa; if($s<0){$count = (0-$s); $s = 0;$pre="N" x $count;}
	$e = $eb - $sa; if($e>$lr){$count = ($e-$lr); $e=$lr; $suf="N" x $count;}
	$s = $s+$sr; $e = $e+$sr; $len = $e-$s; $s = $s+1;
	if($strand eq "-"){
		$txt = substr($rev_seq{$name},$s,$len);	
	}else{
		$txt = substr($seq{$name},$s,$len);	
	}
	$seq = $pre.$txt.$suf;
	if($name =~/(.*)\/[1|2]/){
		$seq_name = $1;
	}
	$region = $a[11];
	#print $seq."\n";
	$final{$seq_name}{$region} = $seq;
	$all_region{$region} = 1;
}
close O;

open O,">$ARGV[2]";
@all_region = sort(keys %all_region);
$n = join "\t",@all_region;
print O "seq_name\t$n\n";
foreach $seq_name (keys %final){
	undef(@n);
	foreach $region (@all_region){
		if($final{$seq_name}{$region}){
			push(@n,$final{$seq_name}{$region});		
		}else{
			push(@n,".");
		}
	}
	$n = join "\t",@n;
	print O "$seq_name\t$n\n";
}
close O;
