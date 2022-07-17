use strict;


my $gtfFile="human.gtf";    #gtf配置文件，里面有Ensembl id对应的基因名字信息
my $gene="BRCA1";
if(exists $ARGV[0]){
	$gene=$ARGV[0];
}

#读取Ensembl表达矩阵文件，将Ensembl id转换为基因名字
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


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
