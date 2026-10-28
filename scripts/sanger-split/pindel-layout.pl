#!/usr/bin/env perl
use strict;
use warnings;
use JSON::PP;
use Sanger::CGP::Pindel::Implement;
for my $method (qw(prepare determine_jobs sanitised_sample_from_bam)) {
 die "Installed cgpPindel lacks required official method: $method\n"
  unless Sanger::CGP::Pindel::Implement->can($method);
}
my $opts = decode_json(do { local $/; <STDIN> });
Sanger::CGP::Pindel::Implement::prepare($opts);
my $result = {version => Sanger::CGP::Pindel::Implement->VERSION,
 tumour_name => $opts->{tumour_name}, normal_name => $opts->{normal_name},
 tumour_dir => Sanger::CGP::Pindel::Implement::sanitised_sample_from_bam($opts->{tumour}),
 normal_dir => Sanger::CGP::Pindel::Implement::sanitised_sample_from_bam($opts->{normal})};
if ($opts->{layout}) {
 my $n = Sanger::CGP::Pindel::Implement::determine_jobs($opts);
 $result->{chromosomes} = [sort keys %{$opts->{seqs}}];
 $result->{jobs} = $n;
}
print encode_json($result);
