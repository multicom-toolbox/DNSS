#!/usr/bin/perl
# Script:    PredictSS.pl
# Author:    Matt Spencer
# Date Made:    10/8/2013
# Last Mod:    10/29/2014
# Lastest Modified: 11/2/2016 by Jie Hou
my $GLOBAL_PATH;
BEGIN { 
$GLOBAL_PATH='/space1/jh7x3/DNSS_release/';
}


use strict;
use lib "$GLOBAL_PATH/lib";
use Time qw(formatted_localtime);
use DN_SSpred qw(make_testfeatures predict_SS timeprint make_dnss_file check_err); 
use SeqAlign qw(load_fasta);
use Getopt::Long;

my ($help, $file, @seq, @name, $indir, $outdir,$device);
my $opt = GetOptions(    '-out:s', \$outdir,
            '-help!', \$help,
            '-seq:s@', \@seq,
            '-name:s@', \@name,
            '-device:s', \$device,
            '-indir:s', \$indir,
            '-file!', \$file );

# Required input is either sequence(s) or an input directory
print_help() if ($help || (!@seq && !$indir));

##################################################
##  Setting up directory architecture
##################################################

$outdir = "." unless ($outdir);
$outdir .= "/" unless ($outdir =~ /\/$/);

# These intermediate directories already exist and are
# designed to remain as are. Files within these temporary
# dirs should be automatically deleted between every 
# prediction, though it shouldn't hurt if they aren't.

my $modeldir = "$GLOBAL_PATH/models/";
my $tempdir = "$outdir/temp/";
my $pssmdir = "$outdir/temp/pssm/";
my $probdir = "$outdir/temp/prob/";
my $featdir = "$outdir/temp/feat/";
my $pfeatdir= "$outdir/temp/pfeat/";
my $pprobdir= "$outdir/temp/pprob/";

`mkdir $outdir` unless (-d $outdir);
`mkdir $tempdir` unless (-d $tempdir);
`mkdir $pssmdir` unless (-d $pssmdir);
`mkdir $probdir` unless (-d $probdir);
`mkdir $featdir` unless (-d $featdir);
`mkdir $pfeatdir` unless (-d $pfeatdir);
`mkdir $pprobdir` unless (-d $pprobdir);

# A log file is generated for each prediction. These logs
# can be found at DNSS/logs/ and are labeled according to
# the time at which this program was called. There is 
# currently no mechanism for deleting old log files.

my $time = formatted_localtime();
$time =~ s/[\[\]]//g;
$time =~ s/\//_/g;
$time =~ s/:/_/g;
$time =~ s/\s/-/;
$time =~ s/ //;
my $logfile = "$GLOBAL_PATH/logs/$time.log";

my $ii=1;
while (-f $logfile){
    $logfile = "$GLOBAL_PATH/logs/$time.$ii.log";
    $ii++;
}

my @DN1s = ($modeldir . "DN1.model.dat", $modeldir . "DN2.model.dat");
my @DN1s_util = ($modeldir . "/util_models/DN1.model.util.dat", $modeldir . "/util_models/DN2.model.util.dat");
my $DN2 = $modeldir . "DN3.model.dat";
my $DN2_util = $modeldir . "/util_models/DN3.model.util.dat";

unless (-d $tempdir){
    print "Temp directory was unable to be created. Aborting.\n";
    timeprint($logfile, "Temp directory was unable to be created. Aborting.");
    exit;
}

timeprint($logfile, "Running PredictSS.pl", 2);

##################################################
##  Adjusting for different input options
##################################################

my @fasta;
if ($file) {            # Input specified as individual files
    @fasta = @seq;
    @name = ();
    foreach (@seq) {
        chomp($_);
        my @paths = split('/', $_);
        my @info = split(/\./, $paths[$#paths]);
        my $fullname = join(/./, @info[0..$#info-1]);
        push (@name, $fullname);
    }
}
elsif ($indir){            # Input specified as a directory of files
    $indir .= "/" unless ($indir =~ /\/$/);

    my $list = `ls $indir`;
    my @files = split(/\n/, $list);
    foreach (@files){
        push (@fasta, $indir . $_) if ($_ =~ /fasta$/);
        my @paths = split('/', $_);
        my @info = split(/\./, $paths[$#paths]);
        my $fullname = join(/./, @info[0..$#info-1]);
        push (@name, $fullname);
    }
}
else {                # Input simply a sequence
    for (my $ii=0; $ii<@seq; $ii++){
        my $thisname;
        if ($name[$ii]){
            $thisname = $name[$ii];    
        }
        else {
            $thisname = sprintf("%04s", $ii);
            push (@name, $thisname);
        }
        my $fastfile = $tempdir . "$thisname.fasta";
        push (@fasta, $fastfile);
        open (OUT, ">$fastfile");
        print OUT $seq[$ii];
        close OUT;
    }
}
if($device)
{
    if($device ne 'GPU' and $device ne 'CPU')
    {
        die "The device parameter should be <GPU/CPU>\n";
    }
    timeprint($logfile, "Calling $device to run prediction!!!");
}else{
    $device = 'CPU';
    timeprint($logfile, "Calling CPU to run prediction!!!");
}

# check if GPU can be called correctly
if($device eq 'GPU')
{
    # check GPU is detected
    timeprint($logfile, "Checking GPU environment!!!");
    my $cmd = $GLOBAL_PATH."/scripts/check_GPU.py 1>$GLOBAL_PATH/logs/GPU.log  2>$GLOBAL_PATH/logs/GPU.Err";
    system("python2 $cmd");
    open(FILE,"$GLOBAL_PATH/logs/GPU.log");
    my $errinfo = 'GPU failed to be detected!';
    if (grep{/$errinfo/} <FILE>){
        timeprint($logfile, "!!!Cudamat failed to run on this computer, setting to CPU to run models");
        $device = 'CPU';
    }else{
        timeprint($logfile, "!!!GPU is detected to run cudamat. GPU enabled!");
        $device = 'GPU';
    }
    close FILE;
}
##################################################
##  SS Prediction on input files
##################################################

my $percent = 0;
my $probfile;
my $predfile;
for (my $ii=0; $ii<@fasta; $ii++){
    if ($ii/@fasta*100 >= $percent && $percent < 100){
        timeprint($logfile, "$percent % complete!");
        $percent += 10;
    }
    remove_intermediates();        # Temporary files removed between each prediction
    timeprint($logfile, "Predicting $name[$ii]...");

    unless (-f $fasta[$ii]){
        check_err("err: Fasta file not found: $fasta[$ii]", $logfile);
    }

    ## Generate PSSM file
    my $pssm = $pssmdir . "$name[$ii].pssm";
    my $log = $pssmdir . "$name[$ii].log";
#    if (-f $pssm){
#        goto SKIP;
#    }

    `$GLOBAL_PATH/scripts/generate-pssm.sh $fasta[$ii] $pssm $log`;
    unless (-f $pssm){
        print "PSSM generation failed; using less stringent parameters.\n";
        `$GLOBAL_PATH/scripts/gen-pssm-less-stringent.sh $fasta[$ii] $pssm $log`;
        unless (-f $pssm){
            check_err("err: PSSM generation failed.", $logfile);
            next;
        }
    }

    SKIP:

    ## Generate Feature file
    my $feat = $featdir . "$name[$ii].feat";

    my @files = (0, $pssm);
    my @params = (1,1,0,0,0);
    my @inputs = (\@params, 17, 1, 0);
    my @featlines = make_testfeatures(\@files, \@inputs);
    next if (check_err($featlines[0], $logfile));

    open (FEAT, ">$feat") or die "Open error: $feat\n";
    foreach (@featlines){
        print FEAT "$_";
    }
    close FEAT or die "Close error: $feat\n";

    ## Predict SS using 1st layer DNs
    my @probfiles;
    my @deletefiles;
    my $err=0;
    if($device eq 'GPU')
    {
        for (my $jj=0; $jj<@DN1s; $jj++){
            my $select_model = $DN1s[$jj];
            timeprint ($logfile, "Accessing model $select_model by GPU");
            my ($file1, $file2) = predict_SS( $DN1s[$jj], $feat, $probdir, 'GPU' , $jj);
            $err++ if (check_err($file1, $logfile));
            push (@probfiles, $file1);
            push (@deletefiles, $file2);
        }
        next if ($err);
    }elsif($device eq 'CPU')
    {
        for (my $jj=0; $jj<@DN1s_util; $jj++){
            my $select_model = $DN1s_util[$jj];
            timeprint ($logfile, "Accessing model $select_model by CPU");
            my ($file1, $file2) = predict_SS( $DN1s_util[$jj], $feat, $probdir, 'CPU' , $jj);
            $err++ if (check_err($file1, $logfile));
            push (@probfiles, $file1);
            push (@deletefiles, $file2);
        }
        next if ($err);
    }else{
        die "Error: The device($device) is not correctly assigned! Please set <GPU/CPU>\n";
    }

    ## Generate Prob-Feat file
    # This file will include features including pssm info
    # as well as the predictions from the 1st layer DNs
    my $pfeat = $pfeatdir . "$name[$ii].feat";
    
    my @files = (0, $pssm, @probfiles);
    my @params = (1,0,0,0,0);
    my @inputs = (\@params, 17, 1, 0);
    my @featlines = make_testfeatures(\@files, \@inputs);
    next if (check_err($featlines[0], $logfile));
    open (PFEAT, ">$pfeat") or die "Open error: $pfeat\n";
    foreach (@featlines){
        print PFEAT "$_";
    }
    close PFEAT or die "Close error: $pfeat\n";

    ## Predict SS using L2 DN
    
    if($device eq 'GPU')
    {
        timeprint ($logfile, "Accessing model $DN2 by GPU");
        ($probfile, $predfile) = predict_SS( $DN2, $pfeat, $pprobdir,'GPU');
        next if (check_err($probfile, $logfile));
    }elsif($device eq 'CPU')
    {
        timeprint ($logfile, "Accessing model $DN2_util by CPU");
        ($probfile, $predfile) = predict_SS( $DN2_util, $pfeat, $pprobdir,'CPU');
        next if (check_err($probfile, $logfile));
    }else{
        die "Error: The device($device) is not correctly assigned! Please set <GPU/CPU>\n";
    }
    
    ## Write output DNSS file
    my $dnss = $outdir . "$name[$ii].dnss";
    my $vdnss = $outdir . "$name[$ii].vdnss";
    my $header = ">$name[$ii]";
    my @dirarray = ($probfile);
    my $return = make_dnss_file (\@dirarray, $pssm, $dnss, $vdnss, $header);
    next if (check_err($return, $logfile));
    
    timeprint ($logfile, "$name[$ii] prediction successful! Saved to $dnss and $vdnss");    
}
timeprint ($logfile, "100% complete!");


# This is used to remove temporary files between each prediction
sub remove_intermediates {
#    `rm -f temp/feat/* temp/pfeat/* temp/pprob/* temp/prob/* temp/pssm/*`;
}


sub print_help{


    print "#################################################################################################################\n";
    print "#                                                                                                               #\n";
    print "#   Software     :  DNSS (A Deep Learning Network Approach to ab initio Protein Secondary Structure Prediction) #\n";
    print "#   Release      :  1.1  (October 2016)                                                                         #\n";
    print "#                                                                                                               #\n";
    print "#   Author(s)    :  Matt Spencer, Jesse Eickholt, and Jianlin Cheng                                             #\n";
    print "#   Maintainance :  Jie Hou(jh7x3\@mail.missouri.edu ), Badri Adhikari (bap54\@mail.missouri.edu)                 #\n";
    print "#   Copyright    :  Bioinformatics, Data Mining, and Machine Learning Lab (BDM)                                 #\n";
    print "#                    Department of Computer Science                                                             #\n";
    print "#                    University of Missouri, Columbia                                                           #\n";
    print "#                                                                                                               #\n";
    print "#################################################################################################################\n";

    print "\nRequired input:\n";
    print "\t-seq    : Sequence of interest (can have multiple inputs, e.g. -seq AAA -seq AAB -seq AAC\n";
    print "\nOutput:\n";
    print "\t.dnss    : File giving the sequence and SS prediction horizontally.\n";
    print "\t.vdnss    : Gives confidence levels for each prediction, additionally.\n";
    print "\nOptions:\n";
    print "\t-name    : Name of sequence (in same order). If files are given, file names will be used.\n";
    print "\t-file    : Indicates that fasta inputs are fasta files instead of sequences.\n";
    print "\t-indir : Indicate a directory full of fastas to predict.\n";
    print "\t-out    : Dictate the location of prediction files (default is . )\n";
    print "\t-device: Specify the device(GPU/CPU) to use (default is CPU )\n";
    print "\t-help  : Print this help message.\n\n";
    print "* There are two ways to indicate the protein to predict:\n";
    print "----------------------------------------------------------------------------\n";
    print "(1) Directly give a sequence:\n";
    print "      \n";
    print "   Usage:\n";
    print "   \$  perl PredictSS.pl -seq <AA sequence> -name <Protein Name> -device <GPU/CPU> -out   -out <output folder>\n";
    print "   \n";
    print "   Example:\n";
    print "   a). CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq GNVVIEVDMANGWRGNASGSTSHSGITYSADGVTFAALGDGVGAVFDIARPTTLEDAVIAMVVNVSAEFKASEANLQIFAQLKEDWSKGEWDALAGSSELTADTDLTLTATIDEDDDKFNQTARDVQVGIQAKGTPAGTITIKSVTITLAQEA -name Prot1  -out ./output/\n";
    print "   b). GPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -device GPU -seq GNVVIEVDMANGWRGNASGSTSHSGITYSADGVTFAALGDGVGAVFDIARPTTLEDAVIAMVVNVSAEFKASEANLQIFAQLKEDWSKGEWDALAGSSELTADTDLTLTATIDEDDDKFNQTARDVQVGIQAKGTPAGTITIKSVTITLAQEA -name Prot1  -out ./output/\n";
    print "   \n";
    print "----------------------------------------------------------------------------\n";
    print "(2) Predict from protein file:\n";
    print "   \n";
    print "   Usage:\n";
    print "   \$ perl PredictSS.pl -seq <file name>.fasta -file -device <GPU/CPU>  -out <output folder>\n";
    print "   \n";
    print "   Example:\n";
    print "   a) CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq test/1GNY-A.fasta -file -device CPU -out output/\n";
    print "   b) Gpu: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -seq test/1GNY-A.fasta -file -device GPU -out output/\n";
    print "          \n";
    print "           \n";
    print "----------------------------------------------------------------------------\n";
    print "   \n";
    print "(3) Predicting multiple proteins:\n";
    print "   \n";
    print "   \n";
    print "   Usage:\n";
    print "   \$ perl PredictSS.pl -indir <input directory> -out <output directory> -device <GPU/CPU>\n";
    print "  \n";
    print "   Example:\n";
    print "   a) CPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -indir ./test/ -out ./output/ -device CPU\n";
    print "   b) GPU: perl /home/jh7x3/DNSS/scripts/PredictSS.pl -indir ./test/ -out ./output/ -device GPU\n";
    print "----------------------------------------------------------------------------\n";
    print "\n";
    print " If you have questions or suggestions, please contact:   \n";     
    print " Jie Hou(jh7x3\@mail.missouri.edu), Jianlin Cheng, PhD\n";
    print " Bioinformatics, Data Mining, and Machine Learning Lab (BDM)\n";
    print " Department of Computer Science\n";
    print " University of Missouri, Columbia\n";
    print " Email: chengji\@missouri.edu\n";

    print " ----------------------------------------------------------------------------\n";

    print " Citation: Spencer, Matt, Jesse Eickholt, and Jianlin Cheng. \"A deep learning network approach to ab initio protein secondary structure prediction.\" IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB) 12.1 (2015): 103-112.\n";

    print "\n\n";
    exit;
}
