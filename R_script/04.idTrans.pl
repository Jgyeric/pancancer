use strict;


my $gtfFile="human.gtf";    #gtf配置文件，里面有Ensembl id对应的基因名字信息

#读取GTF文件，将信息保存到%hash，哈希的Key是Ensembl id，哈希的value是基因名字
my %hash=();
open(RF,"$gtfFile") or die $!;
while(my $line=<RF>){
	chomp($line);
	if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/){
		$hash{$1}=$2;
	}
}
close(RF);

#读取Ensembl表达矩阵文件，将Ensembl id转换为基因名字
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


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
