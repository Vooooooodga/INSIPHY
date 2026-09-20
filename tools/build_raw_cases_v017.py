#!/usr/bin/env python3
"""Synthetic native FASTA/GFF boundary cases for extraction tests.

Fixtures are deliberately not a real-species or statistically calibrated benchmark.
Expected policies describe which evidence exists, not hard-coded inferred events.
"""
from pathlib import Path
import argparse
import json
from random import Random
from insiphy.storage.tabular import write_tsv

CASES={
 'conserved_annotations': 'Paired native exon ranges are preserved.',
 'missing_annotation_split_DNA': 'B has no middle-exon annotation; only DNA/context candidates may be recovered.',
 'annotated_split': 'One A exon versus two B exons with a 55bp interval; no automatic event polarity.',
 'short_complete_alternative': 'Internal alternative start/end must not imply partial.',
 'assembly_gap': 'N-rich middle interval is unknown, not a deletion.',
 'microexon9': '9bp annotated exon requires coding/flank context for homology.',
 'retained_intron_paths': 'Both annotated paths exist; repertoire is not a conflict.',
 'nested_antisense': 'Other-owner exon is context, not focal-gene exonic usage.',
 'internal_repeat': 'Two copies of a local segment are not one-to-two complementary coverage.',
 'UTR_CDS_role': 'A CDS boundary change does not erase exon identity.',
 'translation_exception': 'CDS exception is not automatic structure-gain evidence.',
 'rearranged_candidate': 'Permuted DNA pieces are outside a single monotone event model.'}


def build_raw_cases(root):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    rng=Random(17017);base=''.join(rng.choice('ACGT') for _ in range(1500)); records=[]
    for case,scope in CASES.items():
        dest=root/case;dest.mkdir(exist_ok=True)
        for sp in ('A','B'):
            seq=list(base);exons=[(101,220),(401,520),(901,1080)];extra=[]
            if sp=='B':
                if case=='missing_annotation_split_DNA':
                    seq[400:460]=list(base[400:460]);seq[515:575]=list(base[460:520]);exons=[(101,220),(901,1080)]
                if case=='annotated_split':exons=[(101,220),(401,460),(516,575),(901,1080)]
                if case=='assembly_gap':seq[400:520]=['N']*120;exons=[(101,220),(901,1080)]
                if case=='microexon9':exons=[(101,220),(401,409),(901,1080)]
                if case=='internal_repeat':seq[600:720]=list(base[400:520]);exons=[(101,220),(401,520),(601,720),(901,1080)]
                if case=='rearranged_candidate':seq[400:520]=list(base[900:1020]);seq[900:1020]=list(base[400:520])
            strand='-' if case=='nested_antisense' and sp=='B' else '+'
            lines=['##gff-version 3',f'chr{sp}\tfixture\tgene\t101\t1400\t.\t{strand}\t.\tID=g{sp}',
                   f'chr{sp}\tfixture\tmRNA\t101\t1080\t.\t{strand}\t.\tID=t{sp};Parent=g{sp}']
            for i,(a,b) in enumerate(exons,1):
                lines.append(f'chr{sp}\tfixture\texon\t{a}\t{b}\t.\t{strand}\t.\tID=e{sp}{i};Parent=t{sp}')
            if case=='short_complete_alternative':
                lines += [f'chr{sp}\tfixture\tmRNA\t151\t1050\t.\t+\t.\tID=t{sp}alt;Parent=g{sp}']
                for i,(a,b) in enumerate([(151,220),(401,520),(901,1050)],1):
                    lines.append(f'chr{sp}\tfixture\texon\t{a}\t{b}\t.\t+\t.\tID=alt{sp}{i};Parent=t{sp}alt')
            if case=='retained_intron_paths':
                lines += [f'chr{sp}\tfixture\tmRNA\t101\t1080\t.\t+\t.\tID=t{sp}ret;Parent=g{sp}',
                          f'chr{sp}\tfixture\texon\t101\t520\t.\t+\t.\tID=ret{sp};Parent=t{sp}ret',
                          f'chr{sp}\tfixture\texon\t901\t1080\t.\t+\t.\tID=retend{sp};Parent=t{sp}ret']
            if case=='nested_antisense':
                anti='+' if strand=='-' else '-'
                lines += [f'chr{sp}\tfixture\tgene\t450\t490\t.\t{anti}\t.\tID=other{sp}',
                          f'chr{sp}\tfixture\tncRNA\t450\t490\t.\t{anti}\t.\tID=otherRNA{sp};Parent=other{sp}',
                          f'chr{sp}\tfixture\texon\t450\t490\t.\t{anti}\t.\tID=otherex{sp};Parent=otherRNA{sp}']
            if case in ('UTR_CDS_role','translation_exception'):
                for i,(a,b) in enumerate(exons,1):
                    start=a+(30 if i==1 else 0)+(3 if case=='UTR_CDS_role' and sp=='B' and i==1 else 0)
                    attrs=f'ID=cds{sp}{i};Parent=t{sp}'
                    if case=='translation_exception':attrs+=';exception=ribosomal%20slippage'
                    lines.append(f'chr{sp}\tfixture\tCDS\t{start}\t{b}\t.\t+\t0\t{attrs}')
            (dest/f'{sp}.fa').write_text(f'>chr{sp}\n'+''.join(seq)+'\n')
            (dest/f'{sp}.gff3').write_text('\n'.join(lines)+'\n')
        (dest/'scope.txt').write_text(scope+'\nSynthetic native input only; not a certified evolutionary reconstruction.\n')
        records.append({'case':case,'files':'A.fa;A.gff3;B.fa;B.gff3','scope':scope})
    write_tsv(root/'cases.tsv',records,['case','files','scope'])
    return root

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    print(build_raw_cases(p.parse_args().output))
