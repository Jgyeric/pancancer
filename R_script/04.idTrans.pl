use strict;


my $gtfFile="human.gtf";    #gtf�����ļ���������Ensembl id��Ӧ�Ļ���������Ϣ

#��ȡGTF�ļ�������Ϣ���浽%hash����ϣ��Key��Ensembl id����ϣ��value�ǻ�������
my %hash=();
open(RF,"$gtfFile") or die $!;
while(my $line=<RF>){
	chomp($line);
	if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/){
		$hash{$1}=$2;
	}
}
close(RF);

#��ȡEnsembl�������ļ�����Ensembl idת��Ϊ��������
#my @samp1e=(localtime(time));
opendir(RD,".") or die $!;
while(my $file=readdir(RD)){
	
	if($file=~/TCGA\-(.+)\.htseq\_fpkm\.tsv$/){
		my $outFile="symbol.$1.txt";
		open(RF,"$file") or die $!;
		open(WF,">$outFile") or die $!;
		MARK:while(my $line=<RF>){
			if($.==1){
				print WF $line;
				next MARK;
			}
			chomp($line);
			my @arr=split(/\t/,$line);
			$arr[0]=~s/(.+)\..+/$1/g;
			if(exists $hash{$arr[0]}){
				$arr[0]=$hash{$arr[0]};
				print WF join("\t",@arr) . "\n";
			}
		}
		close(WF); 
		close(RF);
    }
}
close(RD);


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio
