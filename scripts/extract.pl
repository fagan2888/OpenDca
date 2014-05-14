#!/usr/bin/perl

use strict;
use warnings;

my ($label,$counter) = @ARGV;
defined($counter) or die "USAGE: $0 label counter\n";

my $c = 0;

while (<STDIN>) {
	if (/^$label/) {
		last if ($c == $counter);
		$c++;
	}
}

print "$label\n";
while(<STDIN>) {
	last if (/^[a-zA-Z\#]+/);
	s/[\(\)]//g;
	s/,/ /g;
	print;
}
