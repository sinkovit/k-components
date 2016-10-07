#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

my %mdhash;
while(<>) {
    chomp;
    my ($md5sum, $filename) = split;
    $mdhash{$md5sum} = $filename;
}

foreach my $key (keys %mdhash) {
    print "$key $mdhash{$key}\n";
    copy($mdhash{$key}, "UNIQUE");
}
