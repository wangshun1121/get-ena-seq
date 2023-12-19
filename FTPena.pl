#!/usr/bin/env perl
# Wangshun @ Xiamen
# Version 1.1
# 2020.9.2

# 专为ENA设计的下载脚本，不再考虑NCBI的SRA了
# 输入的文件类型，可以是SRR的编号，SRR的列表，或者直接是project ID
# 下载文件格式有三种：fastq，sra和cram_index
# 下载方式有两种：ftp和aspera
# 只下载文件，文件格式转换另外想办法

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;
use POSIX ":sys_wait_h";
use FindBin qw($Bin);
use Parallel::ForkManager;
use Parallel::Simple qw(prun);

our $ascp='~/.aspera/connect/bin/ascp';
our $KEY='~/.aspera/connect/etc/asperaweb_id_dsa.openssh';


our $work_dir = getcwd;
our $Core = `grep \"process\" /proc/cpuinfo | wc -l `;
chomp($Core); # Core number of server

my $id;
my $prj;
my $list;
my $link;
my $datatype='fastq';
my $ftp; # 指定ftp下载！
my $outdir='.';
my $single=0; # 是否是单端数据
my $threads=4;
my $dryrun;
my $help;

GetOptions(
  'i|id=s' => \$id,
  'p|prj=s' => \$prj,
  'l|list=s' => \$list,
  'd|data=s' => \$datatype,
  'o|outdir=s' => \$outdir,
  't|threads=i' => \$threads,
  'ftp' => \$ftp,
  'dry' => \$dryrun,
  'h|help' => \$help,
);

my $Check=`$ascp -h`;
unless($Check){
  unless($ftp){
    $ftp=1; # 没有安装aspera，则只有ftp才能使用
    print STDERR "Warning: Aspera hasn't been installed, only ftp can use! \n";
    print STDERR "To supress this warning, add -ftp in your $0 command line!\n";
    print STDERR "To install aspera, see https://gitee.com/wangshun1121/get-ena-seq \n";
  }
}

my $tool='aspera';
if($ftp){$tool='ftp';}

if($datatype eq 'cram'){$datatype='cram_index';}

my $usage=<<USAGE;

Usage:
  perl $0 -i SRRXXXXXX -d sra
  perl $0 -l SraAccList.txt -o ./sequences -p 5 -t 8
  perl $0 -p PRJNAXXXXX

  -i|-id	<str> 	   SRA accession ID
  -l|-list 	<file>	   SRA ID list, all IDs should be in one column
  -p|-prj       <str>      SRA project ID, such as PRJNAXXXXXX or PRJEBXXXXXX

  -d|data       <str>      Download data type, can be fastq, sra or cram [$datatype]
  -ftp          <!>        Add -ftp, then use ftp rather than aspera to download data
                          If aspera was not installed, only ftp can choose, and this parameter will be ignored

  -o|-outdir	<dir>	   Output directory [working directory $work_dir]
  -t|threads	<int>	   Threads number used for multi detasets downloading [$Core at most, $threads default]
  -dry          <!>        Typig -dry, then only download commandline will be printed.
  -h|-help                 Show this message

USAGE

if ($help) {
	print $usage;
	exit;
}

unless(($datatype eq 'cram_index') or ($datatype eq 'fastq') or ($datatype eq 'sra')){
  print STDERR "Wrong Data Type: -data can only be fastq or sra or cram, $datatype can't be accepted!\n";
  exit();
}

unless ($id or $list or $link or $prj) {
	print $usage;
	exit;
}

unless(-e $outdir){
  system("mkdir -p $outdir");
}

if($id){
  &download($id,$outdir);
  exit;
}

if($link){
  #帮助文件中暂时把这部分先隐藏吧
  &linkdownload($link,$outdir);
  exit;
}

if($prj){
  $list="$outdir/$prj.list";
  open I,"curl -s \"https://www.ebi.ac.uk/ena/portal/api/links/study?accession=$prj&result=read_run\" -H \"accept: */*\"|cut -f 1 |";
  <I>; # 直接打印SRR的列表即可
  open O,">$list";
  while(<I>){
    print O $_;
  }
  close O;
  close I;
}

if($list){
  my @ids=();
  open I,"<$list" or die("No $list, please check!\n");
  while (<I>) {
    chomp;
    s/\r//g; #Windows换行符考虑兼容
    $_=(split/\t/)[0];
    push(@ids,(split/\t/)[0]);
  }

  if($threads>scalar(@ids)){$threads=scalar(@ids);}

  my $pm=new Parallel::ForkManager($threads);
  foreach my $id (@ids) {
  	my $pid=$pm->start and next;
  	&download($id,$outdir);
  	$pm->finish;
  };
  $pm->wait_all_children;
}

if(-e "$outdir/md5sum"){ # 自动去除重复的md5
  system("sort -k 2 $outdir/md5sum|uniq >$outdir/md5");
  system("rm -f $outdir/md5sum");
}

sub linkdownload{
  #直接给出下载链接，下载SRA数据
  #单纯下载数据，不再直接转换SRA到fastq了
  #这里先暂时用wget来代替。
  my $link=shift;
  my $outdir=shift;
  chdir($outdir);
  my $CMD="wget $link";
  chdir($work_dir);
}

sub download{
  # 下载SRA的核心函数
  my $id=shift;
  my $outdir=shift;

  my $link=();

    # print STDERR "Warning: downloading single end sequences from ENA are not taken into consideration\n";
    # ENA可直接下载fq数据，本版本中暂时仅考虑双端序列的情况
    # $link='era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq';
    #system("curl -s \"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$id&fields=all&format=json&limit=0&result=read_run\"|sed \"s\/\\\",\\\"\/\\\",\\n\\\"\/g\"  >$outdir/$id.json");
    system("curl -s \"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$id&fields=all&limit=0&result=read_run\"|sed \"s\/\\\",\\\"\/\\\",\\n\\\"\/g\"  >$outdir/$id.tsv");
    # 下载记录SRR详细信息的json，后面备用
    # print "curl -s \"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$id&result=read_run&fields=$datatype\_$tool,$datatype\_md5\"\n";
    my $WebInfo=`curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$id&result=read_run&fields=$datatype\_$tool,$datatype\_md5"`;
    unless($WebInfo){
      return 0;
      print "No data for $id in $datatype format!\n";
    }
    # 突然发现ENA的API又变了！ 2023.12.19
    # $WebInfo=(split/\n/,$WebInfo)[1];
    # my ($Links,$md5Line)=(split/\t/,$WebInfo)[-2,-1];
    my %InfoHash = &WebInfo2Hash($WebInfo);
    my $Links = $InfoHash{"$datatype\_$tool"};
    my $md5Line = $InfoHash{"$datatype\_md5"};

    # $WebInfo=(split/\n/,$WebInfo)[1];
    # my ($Links,$md5Line)=(split/\t/,$WebInfo)[-2,-1];
    # 修改ENA Reads获取的API，2020.8.17
      
    my @Links=split/;/,$Links;
    my @md5Values=split/;/,$md5Line;
    my %md5=();
    # note：发现了ENA上居然有数据即保存为单端又保存为双端，例如：SRR6315113
    #       应对这种数据的策略，就是直接把数据判定为双端来下载！
    if(scalar(@Links)>1){
      $link=$Links[-2];
    }else{
      $link=$Links[0];
    }
   
    for(my $i=0;$i<scalar(@Links);$i++){
      if($i<scalar(@Links)-2){next;} # 补上这句话就是应付类似 SRR6315113 的情况
      $md5{$Links[$i]}=$md5Values[$i];
      #print "$Links[$i]\t$md5Values[$i]\n";
    }

    my $AutoSingle=$single;
    if($datatype eq 'fastq'){
      if(scalar(keys %md5)>1 ){$AutoSingle=0;}
      else{ # cran或sra则自动单端吧
        $AutoSingle=1;
      }
    }else{ # cran或sra则自动单端吧
      $AutoSingle=1;
    }
    

    if($AutoSingle){
      # 直接包含了fastq，sra和cram的情况
      my $md5Value=$md5{$link};
      my $fileName=(split/\//,$link)[-1];
        $md5{$fileName}=$md5Value;
      my $CMD="$ascp -QT -l 300m -P33001 -i $KEY -k 1 era-fasp\@$link .";
      system("echo \"$md5Value $fileName\" >> $outdir/md5sum");
      if($ftp){
        $CMD = "wget -c ftp://$link";
      }
      
      if($dryrun){
        print "$CMD\n";
      }else{
        chdir($outdir);
        &GetFile($CMD,$fileName,$md5{$fileName});
        chdir($work_dir);
      }
      
    }else{
      # 只有双端Reads
      my $md5Value=$md5{$link};
      system("echo \"$md5Value $id\_1.fastq.gz\" >> $outdir/md5sum");
      $md5{"$id\_1.fastq.gz"}=$md5Value;
      my $R2=$link;
      $R2=~s/\_1.fastq.gz/\_2.fastq.gz/;
      $md5Value=$md5{$R2};
      system("echo \"$md5Value $id\_2.fastq.gz\" >> $outdir/md5sum");
      $md5{"$id\_2.fastq.gz"}=$md5Value;
      $link=~s/\_1.fastq.gz//;
      my $CMD1="$ascp -QT -l 300m -P33001 -i $KEY -k 1 era-fasp\@$link\_1.fastq.gz .";
      my $CMD2="$ascp -QT -l 300m -P33001 -i $KEY -k 1 era-fasp\@$link\_2.fastq.gz .";
      if($ftp){
        $CMD1="wget -c ftp://$link\_1.fastq.gz";
        $CMD2="wget -c ftp://$link\_2.fastq.gz";
      }
      
      #开双线程下载
      if($dryrun){
        print "$CMD1\n";
        print "$CMD2\n";
      }else{
        chdir($outdir);
        prun(
          sub1=>[\&GetFile,($CMD1,"$id\_1.fastq.gz",$md5{"$id\_1.fastq.gz"})],
          sub2=>[\&GetFile,($CMD2,"$id\_2.fastq.gz",$md5{"$id\_2.fastq.gz"})],
        )or die(Parallel::Simple::errplus());
        chdir($work_dir);
      }
    }

  sub WebInfo2Hash{
    # Web Info 不晓得何时改变了列的顺序，故不能直接用第多少列来代表具体内容
    my $WebInfo = shift; 
    my @a = split/\n/,$WebInfo;
    my @Keys = split/\t/,$a[0];
    my @Infos = split/\t/,$a[1];
    my %Hash = ();
    for(my $i = 0; $i < scalar(@Keys);$i++){
      $Hash{$Keys[$i]} = $Infos[$i];
    }
    return(%Hash);
  }

  sub GetFile{
    # 2019.3.8 添加针对ENA下载数据的MD5校验
    my $CMD=shift; # 文件下载命令
    my $File=shift; # 文件名
    my $md5=shift; # 待校验的md5

    unless(-e $File){&run($CMD);} # 文件没有下载完，则跑
    my $Checked=0;
    my $N=1;
    while($N<5){ # 下载最多尝试5次
        if(-e "$File"){
          $Checked=`md5sum $File`;
          chomp($Checked);
          $Checked=(split/\ /,$Checked)[0];
        }
      if($md5 eq $Checked){
        print STDERR "$File downloaded!\n\n";
        last;
      }
      else{
          $N++;
          if(-e $File){
            system("rm -f $File");
            print STDERR "$File check failed. True md5: $md5; Checked: $Checked\n";
          }
        &run($CMD); # md5校验失败，则删除目标文件，重新下载
        next;
      }
    }
    unless($md5 eq $Checked){
      print STDERR "!!! $File check failed, please check the file manually\n";
    }
  }

  sub run{
  	my $CMD = shift;
    my $log="$CMD\nStart:";
    $log.=localtime();
    $log.="\n";
  	$log.=`$CMD`;
    $log.="End:";
    $log.=localtime();
    $log.="\n";
    print $log;
  }
}
