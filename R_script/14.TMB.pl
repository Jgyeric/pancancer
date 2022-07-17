use strict;
use warnings;

my %hash=();
#my @samp1e=(localtime(time));

opendir(RD,".") or die $!;
while(my $file=readdir(RD)){
	
	if($file=~/TCGA\-(.+?)\..+\_snv\.tsv$/){
		my $type=$1;
		open(RF,"$file") or die $!;
		while(my $line=<RF>){
			next if($.==1);
			next if($line=~/^\n/);
			next if($line=~/^\#/);
			my @arr=split(/\t/,$line);
			chomp($line);
			${$hash{$type}}{$arr[0]}++;
		}
		close(RF);
    }
}
close(RD);

open(WF,">TMB.txt") or die $!;
print WF "id\tTMB\tCancerType\n";
foreach my $key(keys %hash){
	foreach my $sample(keys %{$hash{$key}}){
		my $tmb=${$hash{$key}}{$sample}/38;
		print WF "$sample\t$tmb\t$key\n";
    }
}
close(WF);


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio
