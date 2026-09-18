import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from insiphy.alignment import AlignmentBackendError, global_alignment_stats, protein_pair_alignment


class ProteinAlignmentTests(unittest.TestCase):
    def setUp(self):
        availability = patch("insiphy.alignment.shutil.which", return_value="/mock/mafft")
        self.available = availability.start()
        self.addCleanup(availability.stop)
        self.result = SimpleNamespace(returncode=0, stdout=">query\nMK-WX\n>target\nMKQ-X\n", stderr="")
        process = patch("insiphy.alignment.subprocess.run", return_value=self.result)
        self.run = process.start()
        self.addCleanup(process.stop)

    def test_amino_command_preserves_gaps_unknowns_and_query_target_order(self):
        def mafft(command, **kwargs):
            self.assertEqual(Path(command[-1]).read_text(), ">query\nMKWX\n>target\nMKQX\n")
            return self.result

        self.run.side_effect = mafft
        self.result.stdout = ">target description\nMKQ-\nX\n>query\nmk-\nwx\n"
        aligned = protein_pair_alignment("mkwx", "mkqx", threads=4)
        self.assertEqual(aligned, ("MK-WX", "MKQ-X"))
        self.available.assert_called_once_with("mafft")
        self.run.assert_called_once()
        command = self.run.call_args.args[0]
        self.assertEqual(command[:-1], ["/mock/mafft", "--quiet", "--thread", "4", "--amino", "--auto"])

    def test_terminal_alignment_gaps_are_retained(self):
        self.result.stdout = ">query\n--MKWX\n>target\nMKQX--\n"
        self.assertEqual(protein_pair_alignment("MKWX", "MKQX"), ("--MKWX", "MKQX--"))

    def test_nucleotide_mafft_behavior_is_preserved(self):
        self.result.stdout = ">query\nACGT\n>target\nACGT\n"
        stats = global_alignment_stats("ACGT", "ACGT", backend="mafft", threads=2)
        self.assertEqual(stats.backend, "mafft")
        self.assertEqual(stats.identity, 1.0)
        self.assertEqual(stats.coverage, 1.0)
        self.assertEqual(stats.aligned_blocks, [(1, 4, 1, 4)])
        command = self.run.call_args.args[0]
        self.assertEqual(command[:-1], ["/mock/mafft", "--quiet", "--thread", "2", "--nuc", "--auto"])

    def test_missing_mafft_is_explicit(self):
        self.available.return_value = None
        with self.assertRaisesRegex(AlignmentBackendError, "MAFFT was requested but is not available"):
            protein_pair_alignment("MKWX", "MKQX")
        self.run.assert_not_called()

    def test_failed_mafft_is_explicit(self):
        self.result.returncode = 1
        self.result.stderr = "MAFFT could not align the supplied proteins"
        with self.assertRaisesRegex(AlignmentBackendError, "MAFFT could not align"):
            protein_pair_alignment("MKWX", "MKQX")
        self.run.assert_called_once()

    def test_output_must_preserve_both_input_sequences(self):
        outputs = [
            "",
            ">query\nMKWX\n",
            ">query\nMKWX\n>target\nMKQ\n",
            ">query\nMK-W\n>target\nMKQX\n",
            ">query\nMKAX\n>target\nMKQX\n",
            ">query\nMKWX\n>target\nMKQ-\n",
        ]
        for output in outputs:
            with self.subTest(output=output):
                self.result.stdout = output
                with self.assertRaises(AlignmentBackendError):
                    protein_pair_alignment("MKWX", "MKQX")

    def test_caller_owns_stop_and_input_gap_handling(self):
        for query, target in [("", "MKX"), ("MKX", ""), ("MK*", "MKX"), ("M*K", "MKX"), ("MKX", "MK*"), ("MK-X", "MKX")]:
            with self.subTest(query=query, target=target):
                with self.assertRaisesRegex(AlignmentBackendError, "stops normalized by the caller"):
                    protein_pair_alignment(query, target)
        self.run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
