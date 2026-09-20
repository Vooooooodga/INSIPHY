"""Native-input preservation for twelve synthetic boundary cases, not event accuracy."""
import sys
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from build_raw_cases_v017 import build_raw_cases,CASES
from intraphy.preprocess import extract_gene
from intraphy.storage.tabular import read_tsv

class NativeInputsV017(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.root=build_raw_cases(Path(cls.temp.name)/'raw')
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def run_case(self,name):
        results=[]
        for sp in ('A','B'):
            out=self.root/name/('extracted_'+sp)
            extract_gene(self.root/name/f'{sp}.fa',self.root/name/f'{sp}.gff3','g'+sp,name,sp,'g'+sp,out)
            occurrences=read_tsv(out/'segment_occurrences.tsv');raw=read_tsv(out/'raw_gene_features.tsv')
            self.assertTrue(occurrences);self.assertTrue(raw)
            self.assertTrue(all(int(r['end'])>=int(r['start']) for r in occurrences))
            results.append((occurrences,raw,read_tsv(out/'transcript_paths.tsv')))
        return results


def make_case(name):
    def test(self):
        results=self.run_case(name)
        if name=='short_complete_alternative':
            for _,_,paths in results:
                self.assertTrue(any('alt' in p.get('transcript_id','') for p in paths))
                self.assertTrue(all(p.get('partial_start')!='1' and p.get('partial_end')!='1' for p in paths))
        if name=='microexon9':
            self.assertTrue(any(int(r['end'])-int(r['start'])+1==9 and r.get('role')!='intron' for r in results[1][0]))
        if name=='nested_antisense':
            self.assertTrue(any('other' in r.get('id','') for r in results[1][1]))
        if name=='retained_intron_paths':
            self.assertGreaterEqual(len({r['transcript_id'] for r in results[1][2]}),2)
        if name=='missing_annotation_split_DNA':
            self.assertFalse(any(r.get('type')=='exon' and int(r['start'])==401 for r in results[1][1]))
        if name=='annotated_split':
            self.assertEqual(len([r for r in results[1][1] if r.get('type')=='exon']),4)
    test.__name__='test_raw_'+name
    return test

for case in CASES:setattr(NativeInputsV017,'test_raw_'+case,make_case(case))

if __name__=='__main__':unittest.main()
