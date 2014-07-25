#!/usr/bin/perl

use strict;
use warnings;
use Math::Complex;

my ($file,$it,$orbitals,$beta,$matsubaras,$begin,$step) = @ARGV;
defined($step) or die "USAGE: $0 file it orbitals beta matsubaras begin step\n";

my $c0 = int($matsubaras/2);

my $outfile = createFile($file,"gc",$it);
system("perl ../scripts/extract.pl \"#GCKFSC\" $it < $file > $outfile");
commentLines($outfile,1);
addFreq($outfile,"imag");

for (my $i = 0; $i < $orbitals; ++$i) {
	$outfile = createFile($file,"gf0$i",$it);
	system("perl ../scripts/extract.pl \"#GFTMP0 $i\" $it < $file > $outfile");
	addFreq("$outfile","imag");
}


$outfile = createFile($file,"gamma",$it);
system("perl ../scripts/extract.pl \"#gammaFreq\" $it < $file > $outfile");
commentLines("$outfile",1);
addFreq("$outfile","imag");

$outfile = createFile($file,"gf",$it);
system("perl ../scripts/extract.pl \"#gfcluster\" $it < $file > $outfile");
commentLines("$outfile",1);
addFreq("$outfile","real");

$outfile = createFile($file,"G0",$it);
system("perl ../scripts/extract.pl \"#G0\" $it < $file > $outfile");
#commentLines("$outfile",1);
addFreq("$outfile","imag");

$outfile = createFile($file,"oneoverdata2",$it);
system("perl ../scripts/extract.pl \"#ONEOVERDATA2\" $it < $file > $outfile");
deleteColumn("$outfile",0);

$outfile = createFile($file,"sigma",$it);
system("perl ../scripts/extract.pl \"#SIGMA\" $it < $file > $outfile");
commentLines("$outfile",1);
addFreq("$outfile","imag");

$outfile = createFile($file,"delta",$it);
system("perl ../scripts/extract.pl \"#DELTAOMEGA\" $it < $file > $outfile");
deleteColumn("$outfile",0);
#commentLines("$outfile",1);
addFreq("$outfile","imag");

$outfile = createFile($file,"gfMatsubara",$it);
system("perl ../scripts/extract.pl \"#gfMatsubara\" $it < $file > $outfile");
deleteColumn("$outfile",0);
#commentLines("$outfile",1);
addFreq("$outfile","imag");

sub commentLines
{
	my ($file,$l) = @_;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	open(FOUT,"> out.tmp") or die "$0: Cannot open > file out.tmp : $!\n";
	my $c = 0;
	while (<FILE>) {
		if (/^#/) {
			print FOUT;
			next;
		}

		chomp;
		print FOUT "# $_\n";
		$c++;
		last if ($c >= $l);
	}

	while (<FILE>) {
		print FOUT;
	}

	close(FOUT);
	close(FILE);
	system("cp out.tmp $file");
}

sub deleteColumn
{
	my ($file,$l) = @_;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	open(FOUT,"> out.tmp") or die "$0: Cannot open > file out.tmp : $!\n";
	my $c = 0;
	while (<FILE>) {
		if (/^#/) {
			print FOUT;
			next;
		}

		chomp;
		my @temp = split;
		my $n = scalar(@temp);
		for (my $i = 0; $i < $n; ++$i) {
			next if ($i == $l);
			print FOUT "$temp[$i] ";
		}

		print FOUT "\n";
	}

	close(FOUT);
	close(FILE);
	system("cp out.tmp $file");
}

sub addFreq
{
	my ($file,$label) = @_;
	open(FILE,$file) or die "$0: Cannot open file $file : $!\n";
	open(FOUT,"> out.tmp") or die "$0: Cannot open > file out.tmp : $!\n";
	my $c = 0;
	while (<FILE>) {
		if (/^#/) {
			print FOUT;
			next;
		}

		my $freq = getFreq($c,$label);
		print FOUT "$freq ";
		print FOUT;
		$c++;
	}

	close(FOUT);
	close(FILE);
	system("cp out.tmp $file");
}

sub getFreq
{
	my ($c,$label) = @_;
	if ($label eq "real") {
		return $begin + $step * $c;
	}

	my $factor = 2*pi/$beta;
	my $c2 = $c - $c0;
	return ($c2 + 1)*$factor if ($c2 >= 0);
	return $c2*$factor;
}

sub createFile
{
	my ($root,$ext,$it) = @_;
	return "${root}_$it.$ext";
}
