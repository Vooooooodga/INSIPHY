#!/usr/bin/env python3
"""Render the editable v18 methods diagram; SVG requires only the standard library."""
from pathlib import Path
from xml.etree.ElementTree import Element, SubElement, ElementTree
import argparse

NS = "http://www.w3.org/2000/svg"
INK, MUTED, LINE = "#172B3A", "#4D6373", "#B7C8D2"
TEAL, BLUE, GOLD = "#007F86", "#315A9B", "#906315"


def draw(path):
    svg=Element('svg',xmlns=NS,width='1500',height='1120',viewBox='0 0 1500 1120',
                role='img',**{'aria-labelledby':'title description'})
    SubElement(svg,'title',id='title').text='IntraPhy: phylogenetic analysis of local gene structure'
    SubElement(svg,'desc',id='description').text=(
        'Genome sequence and supplied annotation support local correspondence. '
        'A frozen matrix separates DNA presence, exon identity and splice junctions. '
        'Parsimony reports elementary changes; known dependence blocks independent CTMC fitting.')
    defs=SubElement(svg,'defs')
    marker=SubElement(defs,'marker',id='arrow',viewBox='0 0 10 10',refX='9',refY='5',
                      markerWidth='7',markerHeight='7',orient='auto-start-reverse')
    SubElement(marker,'path',d='M 0 0 L 10 5 L 0 10 z',fill=MUTED)
    SubElement(svg,'rect',width='1500',height='1120',fill='white')
    def rect(x,y,w,h,fill='white',stroke=LINE,r=12):
        return SubElement(svg,'rect',x=str(x),y=str(y),width=str(w),height=str(h),
                          rx=str(r),fill=fill,stroke=stroke,**{'stroke-width':'1.5'})
    def text(x,y,value,size=21,weight='normal',fill=INK,anchor='start'):
        node=SubElement(svg,'text',x=str(x),y=str(y),fill=fill,
                        **{'font-family':'DejaVu Sans, Arial, sans-serif','font-size':str(size),
                           'font-weight':weight,'text-anchor':anchor})
        node.text=value
    def lines(x,y,values,size=20,step=29,fill=MUTED):
        for i,value in enumerate(values): text(x,y+step*i,value,size,fill=fill)
    def arrow(points,dashed=False):
        SubElement(svg,'polyline',points=' '.join(f'{x},{y}' for x,y in points),fill='none',
                   stroke=MUTED,**{'stroke-width':'2','marker-end':'url(#arrow)',
                                   'stroke-dasharray':'6 5' if dashed else 'none'})
    def panel(x,w,letter,title,color):
        rect(x,125,w,405,fill='#F8FBFD')
        rect(x,125,w,59,fill=color,stroke=color)
        text(x+23,164,letter,25,'bold','white')
        text(x+60,163,title,23,'bold','white')

    text(52,58,'IntraPhy',37,'bold')
    text(253,58,'Phylogenetic analysis of local gene structure',28)
    text(53,94,'Sequence evidence → structural characters → conditional reconstruction',22,fill=MUTED)
    panel(50,420,'A','Input and annotation',TEAL)
    panel(540,420,'B','Local correspondence',BLUE)
    panel(1030,420,'C','Structural characters',TEAL)
    text(73,222,'Single-copy orthologous genes',22,'bold')
    lines(73,258,['Genome FASTA + GFF3 / GTF','Supplied rooted species tree'],21)
    SubElement(svg,'line',x1='75',y1='315',x2='444',y2='315',stroke=LINE)
    text(73,350,'Retain native structure',22,'bold')
    lines(73,386,['Exons, CDS / UTR and introns','Strand and genomic coordinates','All supplied transcript paths'],21)
    text(73,501,'Annotation remains observation-dependent.',16,fill=MUTED)
    text(563,222,'Complementary evidence',22,'bold')
    lines(563,258,['Protein alignment → codon projection','DNA alignment in bounded intervals','Unique flanks + local order'],20)
    text(563,376,'Resolve actual aligned blocks',22,'bold')
    lines(563,412,['Distinguish fragments from repeats','Retain alternative correspondences','Unresolved evidence stays unknown'],20)
    text(563,501,'A reference sequence is not an ancestor.',16,fill=MUTED)
    arrow([(474,328),(535,328)])
    arrow([(964,328),(1025,328)])
    labels=[('DNA presence','present / absent'),('Exon identity','exonic / non-exonic'),
            ('Splice junction','present / absent')]
    for i,(label,states) in enumerate(labels):
        y=215+i*72
        text(1053,y,label,23,'bold')
        text(1053,y+28,states,20,fill=MUTED)
    text(1053,445,'Unknown ≠ absent',22,'bold',TEAL)
    lines(1053,477,['DNA absence makes exon role inapplicable.','Flanking exons do not prove whole-intron homology.'],15,step=24)

    rect(50,565,1400,89,fill='#EAF4F5',stroke='#7FAFB3')
    text(750,600,'Character catalogue + frozen observation matrix',26,'bold',TEAL,'middle')
    text(750,632,'Actual intervals · state definitions · dependence · observation masks · complete coordinates retained',20,fill=MUTED,anchor='middle')
    arrow([(1240,535),(1240,560)])
    arrow([(397,655),(397,700)])
    arrow([(1103,655),(1103,700)])
    rect(50,706,685,262,fill='#F8FBFD')
    rect(765,706,685,262,fill='#FBFAF6')
    text(75,748,'D  Maximum parsimony',27,'bold',BLUE)
    text(75,782,'Default · fixed rooted tree · equal transition costs',20,fill=MUTED)
    lines(75,827,['All globally optimal states and branch placements',
                   'Required / possible elementary character changes',
                   'One compatible minimum-change reconstruction'],21,step=32)
    text(75,942,'Alternative placements are not added as events.',19,'bold',BLUE)
    text(790,748,'E  Conditional CTMC analysis',27,'bold',GOLD)
    text(790,782,'Optional · same matrix and tree · ER / ARD / foreground',19,fill=MUTED)
    lines(790,827,['Known linked included characters → no independent fit',
                   'Valid fits only: rates, model tests and conditional states',
                   'Finite-sample calibration remains unassessed'],20,step=32)
    text(790,942,'Model-test significance does not validate a named event.',18,'bold',GOLD)
    arrow([(397,970),(397,1003)])
    arrow([(1103,970),(1103,1003)])
    rect(50,1008,1400,79,fill=INK,stroke=INK)
    text(750,1040,'Local evidence, elementary changes and target-specific structural histories',24,'bold','white','middle')
    text(750,1070,'No inferred mutation-event total · no compound events · no complete ancestral transcript or mechanism',19,fill='#D7E6EE',anchor='middle')
    ElementTree(svg).write(path,encoding='utf-8',xml_declaration=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[1]/'docs/figures')
    parser.add_argument('--png',action='store_true',help='Also render a 300-dpi PNG; requires CairoSVG.')
    args=parser.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    target=args.output_dir/'architecture.svg'; draw(target)
    if args.png:
        import cairosvg
        cairosvg.svg2png(url=str(target),write_to=str(target.with_suffix('.png')),
                        output_width=3000,output_height=2240)
    print(target)

if __name__=='__main__':
    main()
