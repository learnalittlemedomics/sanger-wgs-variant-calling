#!/usr/bin/env python3
"""Exercise the actual Perl index against a bounded exhaustive pair oracle."""
import itertools
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import unittest

SOURCE = Path(__file__).resolve().parents[2] / 'scripts/sanger-split/get_abs_bkpts_from_clipped_reads.pl'

class ValidatedScope(unittest.TestCase):
    def test_analysis_body_matches_validated_baseline(self):
        # Pin the precise scope of the equivalence argument. Any intentional
        # future analysis change needs its own review and a new baseline.
        text = SOURCE.read_text(encoding='utf-8')
        before, rest = text.split('# BEGIN ordered chromosome-pair index\n', 1)
        _, after = rest.split('# END ordered chromosome-pair index\n', 1)
        original = before + 'for my $i (0..($#lines - 1)) {\n  for my $j (($i+1)..$#lines) {\n' + after
        self.assertEqual(hashlib.sha256(original.encode()).hexdigest(),
                         '26fb9ac5923e66cf9d86bb8be27e19c7462ad6a08b7fc28e5705b80664638aa1')

@unittest.skipUnless(shutil.which('perl'), 'Perl required')
class OrderedIndex(unittest.TestCase):
    def actual(self, rows):
        text = SOURCE.read_text(encoding='utf-8')
        index = text.split('# BEGIN ordered chromosome-pair index\n', 1)[1].split('# END ordered chromosome-pair index', 1)[0]
        # Execute the production loop headers, replacing only the analysis body
        # with its visited (i,j) trace. Input strings are decoded by JSON::PP.
        program = 'use strict; use warnings FATAL => "all"; use JSON::PP; my @lines=@{decode_json(do {local $/; <STDIN>})};\n'
        program += index + 'print "$i,$j\\n"; }}\n'
        result = subprocess.run(['perl', '-e', program], input=json.dumps(rows), text=True, capture_output=True, check=True)
        return [tuple(map(int, line.split(','))) for line in result.stdout.splitlines()]

    def check_case(self, keys):
        rows = [[a, 0, 1, b, 0, 1, str(i), 2, '+', '-'] for i, (a,b) in enumerate(keys)]
        expected = [(i,j) for i in range(len(rows)) for j in range(i+1,len(rows)) if keys[i] == keys[j]]
        self.assertEqual(self.actual(rows), expected)

    def test_edge_cases(self):
        for keys in ([], [('chr1','chr2')], [('chr1','chr2')]*8,
                     [('chr1','chr2'),('chr2','chr1')]*6,
                     [('a|b','c'),('a','b|c'),('a|b','c')],
                     [('chr1','chr1'),('chrX','chrY'),('chr1','chr1')]):
            with self.subTest(keys=keys): self.check_case(keys)

    def test_every_binary_group_assignment_through_length_six(self):
        for n in range(7):
            for labels in itertools.product(range(2), repeat=n):
                self.check_case([('chr1',str(x)) for x in labels])

    def test_shuffled_multiple_groups(self):
        rng=random.Random(91723)
        for n in (10,30,100,300):
            self.check_case([(rng.choice(['chr1','chr2','chrX']), rng.choice(['chr1','chr2','chrY'])) for _ in range(n)])

if __name__ == '__main__': unittest.main()
