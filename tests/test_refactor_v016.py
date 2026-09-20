"""Contract regressions for 0.16. No external alignment executables required."""
import ast
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import typing
import unittest
from unittest.mock import patch

from insiphy.candidate_chain import ChainCandidate, ordered_candidate_chain
from insiphy.coordinates import Interval0
from insiphy.observations.matrix import ObservationMatrix, prepare_observation_matrix
from insiphy.preparation.annotation_index import (
    AnnotationIndex, clear_annotation_cache, iter_annotation, load_annotation_index,
    read_annotation_for_gene,
)
from insiphy.preparation.features import FeatureHierarchy
from insiphy.run_result import RunResult, begin_run, result_model
from insiphy.storage.tabular import iter_tsv, read_tsv, write_tsv


class Refactor016Contracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_empty_intervals_never_overlap(self):
        for empty in (Interval0(0, 0), Interval0(5, 5), Interval0(10, 10)):
            self.assertFalse(empty.overlaps(Interval0(0, 10)))
            self.assertFalse(Interval0(0, 10).overlaps(empty))
        self.assertFalse(Interval0(0, 5).overlaps(Interval0(5, 10)))

    def test_gzip_and_plain_stream_roundtrip(self):
        rows = [{'id': 'x', 'state': 'unknown'}, {'id': 'y', 'state': '0'}]
        for suffix in ('.tsv', '.tsv.gz'):
            path = self.root / ('records' + suffix)
            write_tsv(path, (dict(row) for row in rows), ['id', 'state'])
            iterator = iter_tsv(path)
            self.assertIs(iter(iterator), iterator)
            self.assertEqual(list(iterator), rows)
            self.assertEqual(read_tsv(path), rows)

    def test_tsv_rejects_duplicate_headers(self):
        path = self.root / 'bad.tsv'
        path.write_text('id\tid\nx\ty\n')
        with self.assertRaisesRegex(ValueError, 'duplicate TSV'):
            read_tsv(path)

    def test_tsv_rejects_ragged_rows(self):
        for body in ('id\tstate\nx\n', 'id\tstate\nx\t1\textra\n'):
            path = self.root / 'bad.tsv'
            path.write_text(body)
            with self.assertRaisesRegex(ValueError, 'wrong number'):
                read_tsv(path)

    def test_optional_missing_is_not_required_missing(self):
        path = self.root / 'none.tsv'
        self.assertEqual(read_tsv(path, optional=True), [])
        with self.assertRaises(FileNotFoundError):
            read_tsv(path)

    def test_matrix_snapshot_owns_readonly_rows(self):
        rows = [{'state': 'unknown'}]
        result = ObservationMatrix.snapshot(rows, [], self.root, 'fixture')
        rows[0]['state'] = 'present'
        self.assertEqual(result.rows[0]['state'], 'unknown')
        with self.assertRaises(TypeError):
            result.rows[0]['state'] = 'absent'

    def test_frozen_matrix_does_not_call_builder(self):
        from insiphy.io import write_structural_site_matrix
        src = self.root / 'source'; out = self.root / 'output'
        src.mkdir()
        path = src / 'matrix.tsv'
        rows = [dict(family_id='f', layer='exon_presence', site_id='s', species='A',
                     state='unknown', state_0='absent', state_1='present', evidence='fixture')]
        write_structural_site_matrix(path, rows)
        with patch('insiphy.structural_sites.build_structural_site_matrix', side_effect=AssertionError('must not rebuild')):
            result = prepare_observation_matrix(src, out, path)
        self.assertEqual(result.mode, 'frozen_matrix')
        self.assertEqual(result.rows[0]['state'], 'unknown')
        self.assertTrue((out / 'structural_site_matrix.tsv').exists())

    def test_explicit_run_model_wins_over_stale_files(self):
        (self.root / 'branch_structural_events.tsv').write_text('old\n')
        (self.root / 'structural_changes.tsv').write_text('old\n')
        RunResult('parsimony', 'single-copy', 'species_tree.tsv', ()).write(self.root)
        self.assertEqual(result_model(self.root), 'parsimony')
        RunResult('er-ard', 'single-copy', 'species_tree.tsv', ()).write(self.root)
        self.assertEqual(result_model(self.root), 'er-ard')

    def test_incomplete_result_cannot_render_an_old_success(self):
        RunResult('parsimony', 'single-copy', 'species_tree.tsv', ()).write(self.root)
        begin_run(self.root, 'er-ard', 'single-copy')
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            result_model(self.root)

    def test_conflicting_legacy_results_require_selection(self):
        for name in ('branch_structural_events.tsv', 'structural_changes.tsv'):
            (self.root / name).write_text('header\n')
        with self.assertRaisesRegex(ValueError, 'ambiguous legacy'):
            result_model(self.root)

    def test_visual_status_does_not_interpret_prose(self):
        from insiphy.visualize import visual_status
        row = {'role': 'CDS', 'presence_status': 'present',
               'evidence': 'unknown predicted inferred hidden; these are discussion words'}
        self.assertEqual(visual_status(row, {'membership_call': 'core_member'}, 'exon_like'), 'annotated')
        row['role'] = 'predicted_exon_candidate'
        self.assertEqual(visual_status(row, {}, 'exon_like'), 'predicted')

    def test_probability_display_uses_boolean_field(self):
        from insiphy.visualize import valid_probability_change
        row = {'endpoint_change_probability': '0.8', 'posterior_available': '1',
               'conditioning': 'this wording is not a protocol'}
        self.assertTrue(valid_probability_change(row))
        row['posterior_available'] = '0'
        row['conditioning'] = 'fit_status=success'
        self.assertFalse(valid_probability_change(row))

    def test_gene_names_do_not_assign_history(self):
        from insiphy.case import infer_manifest_copy_role, infer_manifest_source_label
        for symbol in ('jingwei', 'Sdic', 'Adh', 'yellow emperor', 'ordinary'):
            row = {'gene_id': symbol, 'gene_symbol': symbol}
            self.assertEqual(infer_manifest_copy_role(row), 'candidate')
            self.assertEqual(infer_manifest_source_label(row), 'unknown_source')
        self.assertEqual(infer_manifest_copy_role({'copy_role': 'source'}), 'source')

    def test_missing_explicit_tree_is_error(self):
        from insiphy.case import copy_optional_tree
        copy_optional_tree(None, self.root, 'species_tree.tsv')
        with self.assertRaises(FileNotFoundError):
            copy_optional_tree(self.root / 'missing', self.root, 'species_tree.tsv')

    def test_benchmark_retains_distinct_site_identity(self):
        from insiphy.benchmark import benchmark_events
        inp = self.root / 'input'; out = self.root / 'output'; inp.mkdir(); out.mkdir()
        fields = ['family_id', 'event_class', 'site_id']
        truth = [dict(family_id='f', event_class='intron_loss', site_id=s) for s in ('a', 'b')]
        write_tsv(inp / 'truth_events.tsv', truth, fields)
        write_tsv(out / 'candidate_structural_events.tsv', truth[:1], fields)
        result = benchmark_events(inp, out)[0]
        self.assertEqual(result['truth_events'], 2)
        self.assertEqual(result['true_positive'], 1)
        self.assertEqual(result['false_negative'], 1)

    def test_absent_truth_is_not_zero_performance(self):
        from insiphy.benchmark import benchmark_events
        rows = benchmark_events(self.root, self.root / 'out')
        self.assertEqual(rows[0]['status'], 'unavailable')
        self.assertNotIn('precision', rows[0])

    def test_removed_baseline_does_not_fabricate_method_scores(self):
        from insiphy.baseline import evaluate_baselines
        rows = evaluate_baselines(self.root, self.root / 'out')
        self.assertEqual(rows[0]['status'], 'not_a_method_comparison')
        self.assertEqual(rows[0]['score'], 'NA')

    def test_third_near_optimal_path_can_bypass_top_two_shared_node(self):
        # Paths s-a1-x-t, s-a2-x-t and s-y-t have scores 10, 9 and 8.
        def candidate(name, q, t, score):
            return ChainCandidate(name, Interval0(*q), Interval0(*t), score, 'unit')
        candidates = [candidate('s',(0,1),(0,1),0),
                      candidate('a1',(1,2),(1,2),5),
                      candidate('a2',(1,2),(1,2),4),
                      candidate('x',(2,3),(2,3),5),
                      candidate('y',(1,3),(1,3),8),
                      candidate('t',(3,4),(3,4),0)]
        result = ordered_candidate_chain(candidates, score_delta=2, start_ids={'s'}, end_ids={'t'})
        self.assertEqual(result.best_score, 10)
        self.assertEqual(result.retained_ids, frozenset({'s','a1','a2','x','y','t'}))
        self.assertEqual(result.ambiguous_ids, frozenset({'a1','a2','x','y'}))

    def test_disconnected_required_anchors_are_not_a_unique_mapping(self):
        candidates=[ChainCandidate('left',Interval0(0,2),Interval0(4,6),1,'unit'),
                    ChainCandidate('right',Interval0(2,4),Interval0(0,2),1,'unit')]
        result=ordered_candidate_chain(candidates,score_delta=0,start_ids={'left'},end_ids={'right'})
        self.assertEqual(result.ambiguity_status,'no_connecting_chain')
        self.assertEqual(result.ambiguous_ids,frozenset({'left','right'}))

    def _gff(self):
        path = self.root / 'genes.gff'
        path.write_text(
            'chr1\ts\tgene\t100\t400\t.\t+\t.\tID=g1;Name=shared\n'
            'chr1\ts\tmRNA\t100\t400\t.\t+\t.\tID=t1;Parent=g1\n'
            'chr1\ts\tmRNA\t100\t400\t.\t+\t.\tID=t2;Parent=g1\n'
            'chr1\ts\texon\t120\t160\t.\t+\t.\tID=e1;Parent=t1,t2\n'
            'chr1\ts\tCDS\t130\t160\t.\t+\t0\tID=c1;Parent=e1\n'
            'chr1\ts\tgene\t150\t350\t.\t-\t.\tID=g2;Name=shared\n'
            'chr1\ts\texon\t170\t250\t.\t-\t.\tID=e2;Parent=g2\n')
        return path

    def test_annotation_spatial_index_matches_bruteforce(self):
        index = load_annotation_index(self._gff())
        for lo,hi in ((1,99),(120,120),(150,180),(100,400),(400,400),(401,500)):
            expected=[r for r in index.rows if r['seqid']=='chr1' and r['start']<=hi and r['end']>=lo]
            self.assertEqual(index.overlap('chr1',lo,hi),expected)
        self.assertEqual(index.overlap('other',1,100),[])

    def test_one_annotation_parse_for_multiple_gene_queries(self):
        import insiphy.preparation.annotation_index as module
        path=self._gff();clear_annotation_cache()
        with patch.object(module,'iter_annotation',wraps=iter_annotation) as parse:
            read_annotation_for_gene(path,'g1');read_annotation_for_gene(path,'g2');read_annotation_for_gene(path,'g1')
            self.assertEqual(parse.call_count,1)

    def test_annotation_exact_id_wins_and_ambiguous_alias_rejected(self):
        path=self._gff()
        rows,gene,ids,_=read_annotation_for_gene(path,'g1')
        self.assertEqual(gene['id'],'g1')
        self.assertTrue({'g1','t1','t2','e1','c1'}<=ids)
        self.assertIn('g2',{r['id'] for r in rows})  # opposite strand context retained
        with self.assertRaisesRegex(SystemExit,'ambiguously'):
            read_annotation_for_gene(path,'shared')

    def test_annotation_callers_do_not_mutate_cached_rows(self):
        path=self._gff()
        rows,gene,ids,bounds=read_annotation_for_gene(path,'g1')
        gene['attrs']['Name']='modified'
        rows[0]['attrs']['ID']='modified'
        rows2,gene2,_,_=read_annotation_for_gene(path,'g1')
        self.assertEqual(gene2['attrs']['Name'],'shared')
        self.assertEqual(rows2[0]['attrs']['ID'],'g1')

    def test_changed_annotation_file_invalidates_index(self):
        path=self._gff();before=load_annotation_index(path)
        path.write_text(path.read_text()+'chr2\ts\tgene\t1\t5\t.\t+\t.\tID=g3\n')
        after=load_annotation_index(path)
        self.assertIsNot(before,after)
        self.assertIn('g3',after.by_id)

    def test_feature_hierarchy_keeps_multi_parent(self):
        index=FeatureHierarchy(load_annotation_index(self._gff()).rows)
        for transcript in ('t1','t2'):
            ids={r['id'] for r in index.transcript_features(transcript)}
            self.assertTrue({transcript,'e1','c1','g1'}<=ids)

    def test_alignment_type_hints_remain_resolvable(self):
        from insiphy.alignment import AlignmentStats, AlignmentCandidate, AlignmentCandidateSet, AlignmentGap
        for cls in (AlignmentStats,AlignmentCandidate,AlignmentCandidateSet,AlignmentGap):
            self.assertTrue(typing.get_type_hints(cls))

    def test_numeric_kernels_do_not_import_preparation_or_reporting(self):
        package=Path(__file__).resolve().parents[1]/'src/insiphy/inference'
        for name in ('ctmc.py','sankoff.py','posterior.py'):
            tree=ast.parse((package/name).read_text())
            for node in ast.walk(tree):
                if isinstance(node,ast.ImportFrom):
                    self.assertFalse(any(piece in (node.module or '') for piece in ('preparation','annotation','visualize','reporting','structural_sites')),(name,node.module))

    def test_formal_dispatch_import_does_not_load_legacy_numeric_model(self):
        code='import sys; import insiphy.phylogeny; assert "insiphy.experimental.phylogeny" not in sys.modules; assert "insiphy.experimental.tree_model" not in sys.modules'
        completed=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
        self.assertEqual(completed.returncode,0,completed.stderr)

    def test_old_completion_file_not_an_implicit_correspondence_input(self):
        from insiphy.correspondence import infer_correspondence
        inp=self.root/'input';clean=self.root/'clean';dirty=self.root/'dirty'
        for path in (inp,clean,dirty):path.mkdir()
        occurrences=[dict(occurrence_id='a',family_id='f',species='A',gene_copy_id='ga',role='CDS',presence_status='present',start='1',end='6',strand='+',contig='chr1')]
        homology=[dict(homology_id='h',occurrence_id='a',support_type='sequence',confidence='0.9')]
        write_tsv(inp/'segment_occurrences.tsv',occurrences,list(occurrences[0]))
        write_tsv(inp/'segment_homology.tsv',homology,list(homology[0]))
        (inp/'segment_sequences.fasta').write_text('>a\nAAAGGG\n')
        (dirty/'annotation_completion_candidates.tsv').write_text('bad\tbad\nwrong\n')
        infer_correspondence(inp,clean);infer_correspondence(inp,dirty)
        self.assertEqual((clean/'element_correspondence.tsv').read_text(),(dirty/'element_correspondence.tsv').read_text())
        with self.assertRaises(ValueError):
            infer_correspondence(inp,dirty,annotation_completion_path=dirty/'annotation_completion_candidates.tsv')


if __name__=='__main__':
    unittest.main()
