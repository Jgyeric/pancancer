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


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
