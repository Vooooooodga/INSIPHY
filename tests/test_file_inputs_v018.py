"""File-input contracts. These tests do not establish biological accuracy."""
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
import gzip
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from intraphy.cli import main
from intraphy.inputs.selection import resolve_inputs
from intraphy.inputs.resources import pair_resources
from intraphy.preparation.annotation_index import clear_annotation_cache
from intraphy.verification.native_cases import build_native_example, SPECIES


class FileInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = self.root / 'raw'
        build_native_example(self.raw)
        self.addCleanup(clear_annotation_cache)

    def args(self, **kwargs):
        defaults = dict(fasta=[str(self.raw)], gff=[str(self.raw)], orthologs=None,
                        family_id=None, manifest=None, species_tree=str(self.raw/'species_tree.nwk'))
        defaults.update(kwargs)
        return SimpleNamespace(**defaults)

    def rewrite(self, path, text):
        path.write_text(text)
        clear_annotation_cache()

    def test_single_locus_needs_no_table(self):
        (self.raw/'manifest.tsv').unlink()
        selection = resolve_inputs(self.args())
        self.assertEqual(len(selection.rows), 4)
        self.assertEqual({r['family_id'] for r in selection.rows}, {'target_gene'})
        self.assertTrue(all(Path(r['genome_fasta']).is_absolute() for r in selection.rows))

    def test_ortholog_fasta_selects_exact_transcripts(self):
        selected = resolve_inputs(self.args(orthologs=[str(self.raw/'orthologs')]))
        self.assertEqual({r['family_id'] for r in selected.rows}, {'example_gene'})
        self.assertTrue(all('transcript' in r['source_member_ids'] for r in selected.rows))

    def test_protein_is_valid_selector_not_genomic_substitute(self):
        members = self.root/'protein.faa'
        members.write_text(''.join(f'>{sp}_transcript\nMKWVL\n' for sp in SPECIES))
        self.assertEqual(len(resolve_inputs(self.args(orthologs=[str(members)])).rows), 4)
        path = self.raw/'Taxon_A.fa'
        text = path.read_text()
        self.rewrite(path, text[:text.index('\n')+1] + 'MKWVL'*250+'\n')
        with self.assertRaisesRegex(ValueError, 'non-DNA'):
            resolve_inputs(self.args())

    def test_agat_gene_metadata_is_explicit_selector(self):
        members = self.root/'agat.fa'
        members.write_text(''.join(f'>external_{sp} gene={sp}_gene seq_id={sp}_chr\nATG\n' for sp in SPECIES))
        self.assertEqual(len(resolve_inputs(self.args(orthologs=[str(members)])).rows), 4)

    def test_species_qualified_identical_ids(self):
        for sp in SPECIES:
            p=self.raw/f'{sp}.gff3'
            self.rewrite(p, p.read_text().replace(sp+'_gene', 'g').replace(sp+'_transcript','tx'))
        members=self.root/'qualified.fa'
        members.write_text(''.join(f'>{sp}|tx\nATG\n' for sp in SPECIES))
        self.assertEqual(len(resolve_inputs(self.args(orthologs=[str(members)])).rows),4)
        members.write_text('>tx\nATG\n')
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            resolve_inputs(self.args(orthologs=[str(members)]))

    def test_missing_family_taxon_not_imputed_absent(self):
        p=self.raw/'orthologs/example_gene.fa'
        self.rewrite(p, ''.join(f'>{sp}_gene\nATG\n' for sp in SPECIES[:-1]))
        with self.assertRaisesRegex(ValueError, 'one gene locus per species'):
            resolve_inputs(self.args(orthologs=[str(p)]))

    def test_duplicate_member_headers_rejected(self):
        p=self.root/'dup.fa';p.write_text('>Taxon_A_gene\nATG\n>Taxon_A_gene\nATG\n')
        with self.assertRaisesRegex(ValueError,'Duplicate FASTA'):
            resolve_inputs(self.args(orthologs=[str(p)]))

    def test_isoforms_of_same_locus_are_not_paralogs(self):
        for sp in SPECIES:
            p=self.raw/f'{sp}.gff3'
            self.rewrite(p,p.read_text()+f'{sp}_chr\tsynthetic\tmRNA\t101\t800\t.\t+\t.\tID={sp}_tx2;Parent={sp}_gene\n')
        p=self.root/'isoforms.fa'
        p.write_text(''.join(f'>{sp}_gene\nATG\n>{sp}_tx2\nATG\n' for sp in SPECIES))
        selection=resolve_inputs(self.args(orthologs=[str(p)]))
        self.assertEqual(len(selection.rows),4)
        self.assertTrue(all(';' in r['source_member_ids'] for r in selection.rows))

    def test_two_loci_in_one_species_rejected(self):
        p=self.raw/'Taxon_A.gff3'
        self.rewrite(p,p.read_text()+'Taxon_A_chr\tsynthetic\tgene\t20\t80\t.\t+\t.\tID=paralog\n')
        m=self.raw/'orthologs/example_gene.fa';m.write_text(m.read_text()+'>paralog\nATG\n')
        with self.assertRaisesRegex(ValueError,'one gene locus per species'):
            resolve_inputs(self.args(orthologs=[str(m)]))

    def test_whole_annotation_requires_selection(self):
        p=self.raw/'Taxon_A.gff3'
        self.rewrite(p,p.read_text()+'Taxon_A_chr\tsynthetic\tgene\t20\t80\t.\t+\t.\tID=neighbor\n')
        with self.assertRaisesRegex(ValueError,'--orthologs'):
            resolve_inputs(self.args())
        self.assertEqual(len(resolve_inputs(self.args(orthologs=[str(self.raw/'orthologs')])).rows),4)

    def test_reused_locus_across_families_rejected(self):
        p=self.raw/'orthologs/second.fa'
        p.write_text((self.raw/'orthologs/example_gene.fa').read_text())
        with self.assertRaisesRegex(ValueError,'more than one input family'):
            resolve_inputs(self.args(orthologs=[str(self.raw/'orthologs')]))

    def test_multiple_families_cannot_take_single_name(self):
        p=self.raw/'orthologs/second.fa';p.write_text('>x\nATG\n')
        with self.assertRaisesRegex(ValueError,'single ortholog FASTA'):
            resolve_inputs(self.args(orthologs=[str(self.raw/'orthologs')],family_id='fixed'))

    def test_combined_genomic_fasta(self):
        p=self.root/'loci.fa';p.write_text(''.join((self.raw/f'{sp}.fa').read_text() for sp in SPECIES))
        selected=resolve_inputs(self.args(fasta=[str(p)]))
        self.assertEqual(len(selected.rows),4)
        self.assertEqual(len({r['genome_fasta'] for r in selected.rows}),1)

    def test_combined_duplicate_record_id_rejected(self):
        p=self.root/'loci.fa';p.write_text('>chr1\nACG\n>chr1\nACG\n')
        with self.assertRaisesRegex(ValueError,'duplicate record ID'):
            resolve_inputs(self.args(fasta=[str(p)]))

    def test_combined_record_cannot_belong_to_two_species(self):
        p=self.root/'combined.fa';p.write_text((self.raw/'Taxon_A.fa').read_text().replace('Taxon_A_chr','chr1'))
        for sp in SPECIES:
            g=self.raw/f'{sp}.gff3';self.rewrite(g,g.read_text().replace(sp+'_chr','chr1'))
        with self.assertRaisesRegex(ValueError,'distinct sequence IDs'):
            resolve_inputs(self.args(fasta=[str(p)]))

    def test_gzip_files(self):
        for p in list(self.raw.glob('*.fa'))+list(self.raw.glob('*.gff3')):
            with gzip.open(str(p)+'.gz','wt') as handle: handle.write(p.read_text())
            p.unlink()
        self.assertEqual(len(resolve_inputs(self.args()).rows),4)

    def test_tree_required(self):
        with self.assertRaisesRegex(ValueError,'rooted --species-tree'):
            resolve_inputs(self.args(species_tree=None))

    def test_tree_tip_names_must_match(self):
        p=self.raw/'species_tree.nwk';p.write_text(p.read_text().replace('Taxon_A','Other'))
        with self.assertRaisesRegex(ValueError,'species filename stems'):
            resolve_inputs(self.args())

    def test_fasta_gff_stem_mismatch(self):
        (self.raw/'Taxon_A.gff3').rename(self.raw/'Wrong.gff3')
        with self.assertRaisesRegex(ValueError,'same species stem'):
            resolve_inputs(self.args())

    def test_duplicate_resource_path_rejected(self):
        with self.assertRaisesRegex(ValueError,'more than once'):
            pair_resources([str(self.raw),str(self.raw/'Taxon_A.fa')],[str(self.raw)])

    def test_spliced_input_with_original_coords_rejected(self):
        p=self.raw/'Taxon_A.fa';self.rewrite(p,'>Taxon_A_chr\nATG\n')
        with self.assertRaisesRegex(ValueError,'exceed FASTA length'):
            resolve_inputs(self.args())

    def test_gene_only_annotation_not_exon_inferred(self):
        p=self.raw/'Taxon_A.gff3'
        self.rewrite(p,'\n'.join(line for line in p.read_text().splitlines() if '\tgene\t' in line)+'\n')
        with self.assertRaisesRegex(ValueError,'no extractable'):
            resolve_inputs(self.args())

    def test_annotation_without_gene_has_actionable_error(self):
        p=self.raw/'Taxon_A.gff3'
        self.rewrite(p,'\n'.join(line for line in p.read_text().splitlines() if '\tgene\t' not in line)+'\n')
        with self.assertRaisesRegex(ValueError,'normalize the hierarchy'):
            resolve_inputs(self.args())

    def test_family_path_traversal_rejected(self):
        with self.assertRaisesRegex(ValueError,'Invalid family_id'):
            resolve_inputs(self.args(family_id='../escape'))

    def test_check_does_not_require_aligners(self):
        with patch('shutil.which',return_value=None),redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            result=main(['check','--fasta',str(self.raw),'--gff',str(self.raw),
                         '--species-tree',str(self.raw/'species_tree.nwk')])
        self.assertEqual(result,0)

    def test_legacy_table_optional_but_not_mixed(self):
        selection=resolve_inputs(self.args(fasta=None,gff=None,manifest=str(self.raw/'manifest.tsv')))
        self.assertEqual(selection.source_mode,'manifest')
        with self.assertRaisesRegex(ValueError,'not both'):
            resolve_inputs(self.args(manifest=str(self.raw/'manifest.tsv')))


if __name__=='__main__':
    unittest.main()
