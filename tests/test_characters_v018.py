"""Structural characters, coordinate evidence and dependence are separate facts."""
import json
import tempfile
import unittest
from pathlib import Path

from helpers_v018 import observation
from intraphy.observations.characters import (
    add_interval_dependencies, character_id, coordinate_evidence,
    validate_sequence_applicability, write_character_catalogue)
from intraphy.storage.tabular import read_tsv, write_tsv


class CharacterTests(unittest.TestCase):
    def test_identifier_escapes_delimiters(self):
        self.assertNotEqual(character_id(("a/b","c","d")),character_id(("a","b/c","d")))

    def test_catalogue_does_not_merge_equal_patterns(self):
        rows=[observation(site,sp,state) for site in ("x","y") for sp,state in zip("ABC",(0,1,1))]
        with tempfile.TemporaryDirectory() as tmp:
            catalogue=write_character_catalogue(tmp,rows,[])
            self.assertEqual(len(catalogue),2)
            self.assertEqual(len({r['count_unit_id'] for r in catalogue}),2)

    def test_plotting_fragments_do_not_create_characters(self):
        rows=[observation("x",sp,1,layer="exon_presence") for sp in "ABC"]
        coordinates=[dict(family_id="family",layer="exon_presence",site_id="x",
                         species="A",gene_copy_id="g",occurrence_id="e",contig="chr",
                         start=i,end=i+4,strand="+",element_id="x",
                         coordinate_type="genomic_aligned_interval_1_based_closed") for i in (1,6,11)]
        with tempfile.TemporaryDirectory() as tmp:
            first=write_character_catalogue(tmp,rows,coordinates)
            second=write_character_catalogue(tmp,rows,coordinates+coordinates)
            self.assertEqual(first,second)
            self.assertEqual(len(first),1)

    def test_shared_transcript_alone_does_not_group(self):
        rows=[dict(observation(s,"A",1),parent_transcript_ids="same") for s in ("x","y")]
        self.assertTrue(all(r['linked_group_id']=="NA" for r in add_interval_dependencies(rows,[])))

    def test_overlapping_intervals_are_linked_within_layer(self):
        rows=[observation(s,"A",1,layer="exon_presence") for s in ("x","y")]
        coords=[dict(family_id="family",layer="exon_presence",site_id=s,species="A",
                     gene_copy_id="g",contig="chr",strand="+",start=start,end=end,
                     coordinate_type="genomic_aligned_interval_1_based_closed")
                for s,start,end in (("x",1,10),("y",8,20))]
        linked=add_interval_dependencies(rows,coords)
        self.assertEqual(linked[0]['linked_group_id'],linked[1]['linked_group_id'])
        self.assertTrue(linked[0]['linked_group_id'].startswith('overlap:'))

    def test_separate_genomes_do_not_create_overlap(self):
        rows=[observation(s,sp,1,layer="exon_presence") for s,sp in (("x","A"),("y","B"))]
        coords=[dict(r,contig="chr",gene_copy_id="g",strand="+",start=1,end=10,
                     coordinate_type="genomic_aligned_interval_1_based_closed") for r in rows]
        self.assertTrue(all(r['linked_group_id']=="NA" for r in add_interval_dependencies(rows,coords)))

    def test_missing_coordinates_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            cat=write_character_catalogue(tmp,[observation("x","A",1)],[])
            self.assertEqual(cat[0]['coordinates_available'],0)

    def test_frozen_matrix_does_not_recover_current_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp,'element_correspondence.tsv').write_text('invalid_current_annotation\n')
            self.assertEqual(coordinate_evidence([observation("x","A",1)],tmp,frozen=True),[])

    def test_frozen_coordinate_sidecar_is_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            row=dict(observation("x","A",1),start="7",end="12",coordinate_type="test")
            write_tsv(Path(tmp,'character_coordinates.tsv'),[row],list(row))
            self.assertEqual(coordinate_evidence([row],tmp,frozen=True)[0]['start'],'7')

    def test_role_observation_requires_present_dna(self):
        for role in (0,1):
            with self.subTest(role=role), self.assertRaises(ValueError):
                validate_sequence_applicability([observation("x","A",0,layer="exon_presence"),
                    observation("x","A",role,layer="exon_role")])

    def test_unknown_role_is_allowed_when_dna_absent(self):
        validate_sequence_applicability([observation("x","A",0,layer="exon_presence"),
                    observation("x","A","unknown",layer="exon_role")])
