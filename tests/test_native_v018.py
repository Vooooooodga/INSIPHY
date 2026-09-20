"""Real alignment tools recover basic characters from raw synthetic FASTA/GFF3."""
import shutil
import tempfile
import unittest
from pathlib import Path

from intraphy.case import build_case
from intraphy.workflow import run_all
from intraphy.io import read_tsv
from intraphy.evidence.completion import completion_call
from intraphy.verification.native_cases import build_native_example

TOOLS = bool(shutil.which('mafft') and shutil.which('minimap2'))


class EvidenceQualificationTests(unittest.TestCase):
    def test_missing_flanks_do_not_establish_hidden_exon(self):
        row=dict(evidence_status='supports_hidden_segment',inferred_role='CDS')
        self.assertEqual(completion_call(row,0.99,0.55),'ambiguous_evidence')

    def test_explicit_double_flanks_allow_candidate(self):
        row=dict(evidence_status='supports_hidden_segment',inferred_role='CDS',
                 anchor_interval_status='ordered_double_flank_bounded_interval')
        self.assertEqual(completion_call(row,0.99,0.55),'hidden_segment_candidate')

    def test_protein_prediction_keeps_candidate_identity(self):
        row=dict(evidence_status='supports_hidden_segment',inferred_role='predicted_CDS',
                 predicted_role='CDS',inferred_event='protein_cds_projection',
                 anchor_interval_status='ordered_double_flank_bounded_interval')
        self.assertEqual(completion_call(row,0.99,0.55),'predicted_exon_candidate')


@unittest.skipUnless(TOOLS,'Raw integration requires real MAFFT and minimap2')
class NativeIntegrationTests(unittest.TestCase):
    def analyze(self, root, scenario):
        native=root/'native'; case=root/'case'; result=root/'result'
        manifest=build_native_example(native,scenario)
        # The truth sidecar is deliberately removed before analysis.
        (native/'expected.json').unlink()
        build_case(manifest,case,species_tree=native/'species_tree.nwk',threads=1)
        run_all(case,result,threads=1)
        return case,result

    def test_added_intron_from_raw_sequence_is_one_junction_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            case,result=self.analyze(Path(tmp),'splice_difference')
            events=[r for r in read_tsv(result/'branch_structural_events.tsv') if r['event_type']!='no_change']
            required=[r for r in events if r['placement_status']=='required']
            self.assertEqual(len(required),1)
            self.assertEqual(required[0]['event_type'],'intron_gain')
            self.assertEqual(required[0]['child_label'],'Taxon_D')
            self.assertEqual(required[0]['structural_relation'],'exon_split')
            self.assertFalse((result/'compound_structural_events.tsv').exists())
            members=read_tsv(result/'element_correspondence.tsv')
            self.assertFalse(any(r.get('reference_coverage_relation')=='repeated_overlap' for r in members))

    def test_conserved_raw_sequence_has_no_required_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            _,result=self.analyze(Path(tmp),'conserved')
            events=read_tsv(result/'branch_structural_events.tsv')
            self.assertFalse(any(r['event_type']!='no_change' for r in events))

    def test_annotation_dropout_does_not_imply_dna_deletion(self):
        with tempfile.TemporaryDirectory() as tmp:
            _,result=self.analyze(Path(tmp),'annotation_dropout')
            matrix=read_tsv(result/'structural_site_matrix.tsv')
            self.assertFalse(any(r['layer']=='exon_presence' and r['state']=='absent' for r in matrix))
            events=read_tsv(result/'branch_structural_events.tsv')
            self.assertFalse(any(r['event_type']=='sequence_loss' and r['placement_status']=='required' for r in events))
            # Annotation-dependent role and junction differences remain possible.
            # The test does not assert recovery of missing biological RNA evidence.
