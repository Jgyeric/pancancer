use strict;


my $gtfFile="human.gtf";    #gtf�����ļ���������Ensembl id��Ӧ�Ļ���������Ϣ
my $gene="BRCA1";
if(exists $ARGV[0]){
	$gene=$ARGV[0];
}

#��ȡEnsembl�������ļ�����Ensembl idת��Ϊ��������
#my @samp1e=(localtime(time));
opendir(RD,".") or die $!;
open(WF,">singleGeneExp.txt") or die $!;
print WF "Id\t$gene\tType\tCancerType\n";
while(my $file=readdir(RD)){
	
	if($file=~/symbol\.(.+)\.txt/){
		my $type=$1;
		open(RF,"$file") or die $!;
		my @sampleArr=();
		MARK:while(my $line=<RF>){
			chomp($line);
			my @arr=split(/\t/,$line);
			if($.==1){
				@sampleArr=@arr;
			}
			else{
				if($gene eq $arr[0]){
					for(my $i=1;$i<=$#arr;$i++){
						my $sampleNme=$sampleArr[$i];
						my @barcodeArr=split(/\-/,$sampleNme);
						if($barcodeArr[3]=~/^0/){
							print WF "$sampleNme\t$arr[$i]\tTumor\t$type\n";
						}
						else{
							print WF "$sampleNme\t$arr[$i]\tNormal\t$type\n";
						}
					}
					last MARK;
			    }
			}
		}
		close(RF);
    }
}
close(WF);
close(RD);


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio
