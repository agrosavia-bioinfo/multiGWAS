#!/usr/bin/perl -w

use strict;
use File::Basename;

my $dirname = dirname(__FILE__);
my $top = defined($dirname) ? $dirname : '.';

my $libdir = "$top/lib";
opendir (DIR, "$libdir") || die "Could not open $libdir\n";
my @list = readdir(DIR);

my @fl = ();
foreach my $fn(@list){
   if ("$fn" =~ m/\.jar$/){
      push(@fl, "$libdir\/$fn");
   }
}
push(@fl, "$top/sTASSEL.jar");
my $CP = join(":", @fl);
print $CP . "\n";

# Scan @ARGV for Java memory arguments, and put rest in @args
my $java_mem_min_default = "-Xms512m";
my $java_mem_max_default = "-Xmx1536m";
my $java_mem_min = "";
my $java_mem_max = "";
my @args;
for (my $i=0; $i<=$#ARGV; $i++){
   if ($ARGV[$i] =~ m/Xms/) {
      $java_mem_min .= "$ARGV[$i]";
   }
   elsif ($ARGV[$i] =~ m/Xmx/) {
      $java_mem_max .= "$ARGV[$i]";
   }
   else{
      push(@args, $ARGV[$i]);
   }
}

if ($java_mem_min eq "") { $java_mem_min = $java_mem_min_default; }
if ($java_mem_max eq "") { $java_mem_max = $java_mem_max_default; }

print "Memory Settings: $java_mem_min $java_mem_max\n";
print "Tassel Pipeline Arguments: " . "@args\n";

system "java -classpath '$CP' $java_mem_min $java_mem_max net.maizegenetics.pipeline.TasselPipeline @args";
