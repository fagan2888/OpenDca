#!/usr/bin/perl -w
use strict;
use lib '../../PsimagLite/scripts';
use Make;

my @drivers = ("openDca");

my $lapack = Make::findLapack();
my ($gslC,$glsL) = Make::findGsl();
Make::backupMakefile();
writeMakefile();

sub writeMakefile
{
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $libs = "$lapack  $glsL  -lm  -lpthread -lpsimaglite";
	my $cxx = "g++";
	my $cppflags = " -O3 -DNDEBUG  -I../Tpem -IEngine -I../../PsimagLite ";
	$cppflags .= " -I../../PsimagLite/src  -I../../dmrgpp/src/Engine $gslC";
	Make::make($fh,\@drivers,"OpenDca","Linux",0,$libs,$cxx,$cppflags,"true",
	"Engine/Version.h","Engine/Version.h gitrev");

	local *FH = $fh;
print FH<<EOF;

gitrev: gitrev.cpp
	\$(CXX) \$(CPPFLAGS) gitrev.cpp -o gitrev

Engine/Version.h: gitrev
	./gitrev > Engine/Version.h

EOF

	close($fh);
	print STDERR "File Makefile has been written\n";
}

