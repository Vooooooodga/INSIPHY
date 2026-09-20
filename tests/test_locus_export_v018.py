"""Coordinate and optional AGAT interface contracts (not real-AGAT tests)."""
from pathlib import Path
from types import SimpleNamespace
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from intraphy.inputs.agat import normalize_annotations
from intraphy.inputs.loci import export_loci
from intraphy.inputs.selection import resolve_inputs
from intraphy.preparation.annotation_index import clear_annotation_cache, load_annotation_index
from intraphy.storage.fasta import read_fasta_interval


class LocusExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.raw=self.root/'raw';self.raw.mkdir()
        (self.raw/'tree.nwk').write_text('(A:1,B:1)root;\n')
        self.dna=('ACGTCGTA'*800)[:5900]
        for sp in ('A','B'):
            (self.raw/f'{sp}.fa').write_text(f'>{sp}_chr\n{self.dna}\n')
            (self.raw/f'{sp}.gff3').write_text(
                '##gff-version 3\n'+
                f'{sp}_chr\ttest\tgene\t5001\t5600\t.\t-\t.\tID={sp}_gene\n'+
                f'{sp}_chr\ttest\tmRNA\t5001\t5600\t.\t-\t.\tID={sp}_tx;Parent={sp}_gene\n'+
                f'{sp}_chr\ttest\texon\t5401\t5600\t.\t-\t.\tID={sp}_ex;Parent={sp}_tx\n'+
                f'{sp}_chr\ttest\tCDS\t5401\t5600\t.\t-\t2\tID={sp}_cds;Parent={sp}_tx\n')
        self.addCleanup(clear_annotation_cache)
        self.args=SimpleNamespace(fasta=[str(self.raw)],gff=[str(self.raw)],orthologs=None,
            family_id='OG1',manifest=None,species_tree=str(self.raw/'tree.nwk'))

    def export(self,flank=1000):
        self.out=self.root/'out'
        return export_loci(resolve_inputs(self.args),self.args.species_tree,str(self.out),flank)

    def test_negative_gene_sequence_stays_forward(self):
        records=self.export();r=records[0]
        sequence=read_fasta_interval(self.out/'OG1/A.fa','A_locus',1,1900)
        self.assertEqual(sequence,self.dna[4000:5900])
        self.assertEqual(r['source_start'],4001)
        self.assertEqual(r['source_end'],5900)
        self.assertEqual(r['sequence_orientation'],'forward_genomic')

    def test_minus_strand_and_cds_phase_preserved(self):
        self.export();features=load_annotation_index(self.out/'OG1/A.gff3').rows
        gene=next(f for f in features if f['type']=='gene')
        cds=next(f for f in features if f['type']=='CDS')
        self.assertEqual((int(gene['start']),int(gene['end']),gene['strand']),(1001,1600,'-'))
        self.assertEqual(cds['phase'],'2')
        self.assertEqual(int(cds['start'])+4000,5401)
        self.assertEqual(int(cds['end'])+4000,5600)

    def test_actual_transcriptional_flanks_and_boundary(self):
        r=self.export()[0]
        self.assertEqual(r['upstream_available_bp'],300)
        self.assertEqual(r['downstream_available_bp'],1000)
        self.assertEqual(r['flank_status'],'truncated_at_input_boundary')

    def test_exported_pair_roundtrips_without_manifest(self):
        self.export();directory=self.out/'OG1'
        selection=resolve_inputs(SimpleNamespace(fasta=[str(directory)],gff=[str(directory)],
            orthologs=None,manifest=None,family_id=None,species_tree=str(directory/'species_tree.nwk')))
        self.assertEqual(len(selection.rows),2)
        self.assertEqual({r['gene_id'] for r in selection.rows},{'A_gene','B_gene'})

    def test_zero_flank_preserves_full_gene_span(self):
        r=self.export(0)[0]
        self.assertEqual((r['source_start'],r['source_end']),(5001,5600))
        self.assertEqual(r['flank_status'],'complete')

    def test_negative_flank_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonnegative'): self.export(-1)

    def test_optional_agat_missing_is_not_internal_repair(self):
        with patch('intraphy.inputs.agat.shutil.which',return_value=None):
            with self.assertRaisesRegex(RuntimeError,'optional'):
                normalize_annotations([str(self.raw/'A.gff3')],str(self.root/'agat'))

    def test_agat_interface_records_command_and_modified_annotation(self):
        def fake_process(command,**kwargs):
            Path(command[command.index('--output')+1]).write_text((self.raw/'A.gff3').read_text())
            self.assertEqual(kwargs['cwd'],(self.root/'agat').resolve())
            return subprocess.CompletedProcess(command,0,'interface fixture','')
        with patch('intraphy.inputs.agat.shutil.which',return_value='/mock/agat'),patch(
            'intraphy.inputs.agat.subprocess.run',side_effect=fake_process):
            records=normalize_annotations([str(self.raw/'A.gff3')],str(self.root/'agat'))
        self.assertTrue(records[0]['annotation_modified'])
        self.assertIn('not_independent',records[0]['interpretation'])
        self.assertTrue((self.root/'agat/agat_commands.json').exists())

    def test_agat_failure_not_accepted(self):
        with patch('intraphy.inputs.agat.shutil.which',return_value='/mock/agat'),patch(
            'intraphy.inputs.agat.subprocess.run',return_value=subprocess.CompletedProcess([],1,'','bad')):
            with self.assertRaisesRegex(RuntimeError,'AGAT failed'):
                normalize_annotations([str(self.raw/'A.gff3')],str(self.root/'agat'))
        self.assertIn('bad',(self.root/'agat/A.agat.log').read_text())

    def test_agat_timeout_preserves_output(self):
        with patch('intraphy.inputs.agat.shutil.which',return_value='/mock/agat'),patch(
            'intraphy.inputs.agat.subprocess.run',side_effect=subprocess.TimeoutExpired('agat',1,b'partial',b'err')):
            with self.assertRaisesRegex(RuntimeError,'timed out'):
                normalize_annotations([str(self.raw/'A.gff3')],str(self.root/'agat'),timeout=1)
        self.assertIn('partial',(self.root/'agat/A.agat.log').read_text())


if __name__=='__main__': unittest.main()
