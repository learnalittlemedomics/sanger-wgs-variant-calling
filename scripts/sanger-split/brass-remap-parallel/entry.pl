#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
exec('python3', "$FindBin::Bin/../brass-remap-parallel/controller.py", @ARGV) or die "Cannot start parallel remapping: $!";
