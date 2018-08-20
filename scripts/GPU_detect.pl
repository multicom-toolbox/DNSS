#!/usr/bin/perl -w
###############################################################################
# Name %n  : GPU_detect.pl
# Desc %d  : A special script to check if cudamat can be run in computer
#
# Author: Jie Hou
# Date: 2016/02/02
###############################################################################

my $GLOBAL_PATH='/space1/jh7x3/DNSS_release/';

$SCRIPT_DIR=$GLOBAL_PATH."scripts/";

# check GPU is detected
print "Checking GPU\n";
my $cmd = $SCRIPT_DIR."check_GPU.py &>$GLOBAL_PATH/logs/GPU.log";
#print("python2 $cmd \n");
system("python2 $cmd");
open(FILE,"$GLOBAL_PATH/logs/GPU.log");
my $errinfo = 'GPU failed to be detected!';
if (grep{/$errinfo/} <FILE>){
	print "Cudamat failed to run on this computer, setting to CPU to run models\n";
	$device = 'CPU';
}else{
	print "GPU is detected to run cudamat. GPU enabled!\n";
	$device = 'GPU';
}
close FILE;


print "DNSS can be run on ".$device."\n";
