"""Three-step biological overview; detail belongs in the two companion figures."""
from .drawing import Canvas, INK, MUTED, TEAL, PURPLE, AMBER, GREY, PALE, LINE
from .example import TREE as TOY_TREE, CHARACTERS
TOY_STATES = CHARACTERS['intron_presence']

def draw(path):
    p=Canvas('IntraPhy | Reconstructing changes in gene structure',
        'Genome sequences + gene annotations + a rooted species tree, for a set of single-copy orthologs.',
        width=1800,height=820)
    p.section(48,142,'1','Identify homologous regions','Sequence alignment and genomic context')
    p.section(648,142,'2','Define structural characters','One biological question for each aligned region')
    p.section(1248,142,'3','Reconstruct their evolution','Ancestral states, branch changes and rates')
    p.line(615,132,615,700,stroke=LINE,width=1)
    p.line(1215,132,1215,700,stroke=LINE,width=1)
    xy=p.tree(62,258,width=112,step=86)
    # Ribbons are drawn before genes and connect only a confirmed DNA interval.
    p.ribbon(230+104*.96,230+168*.96,269,230+104*.96,230+168*.96,333)
    for s in 'ABCD':
        p.gene(230,xy[s][1]-11,width=336,segment=s!='D',exonic=None if s=='C' else True,
               intron=s!='B')
    p.lines(65,617,['Solid boxes: annotated exons',
                    'Ribbons: aligned homologous sequence',
                    'Dashed box: candidate, not a confirmed exon'],16,step=28)
    for yy,kind,title,left,right in [
        (250,'sequence','Presence of homologous DNA','Present','Deletion supported'),
        (385,'exonic','Exonic status of homologous DNA','Exonic','Intronic*'),
        (520,'intron','Intron at a homologous position','Intron present','Continuous exon')]:
        p.text(660,yy,title,20,bold=True)
        p.local(678,yy+26,kind,1,w=160)
        p.local(988,yy+26,kind,0,w=160)
        p.line(856,yy+37,969,yy+37,stroke=LINE)
        p.text(758,yy+78,left,16,anchor='middle')
        p.text(1068,yy+78,right,16,anchor='middle')
    p.lines(660,647,['Unknown is not absence.',
        '* Intronic only in an informative annotated transcript.'],15,step=25,fill=MUTED)
    xr=p.tree(1265,258,width=106,step=86)
    for s in 'ABCD':
        p.gene(1421,xr[s][1]-11,width=305,segment=s!='D',exonic=None if s=='C' else True,intron=s!='B')
    # Branch annotations refer to elementary changes, not physical mutation counts.
    p.dot(1348,xr['B'][1],PURPLE,5)
    p.text(1408,xr['B'][1]-25,'Intron loss',16,bold=True,fill=PURPLE)
    p.dot(1348,xr['D'][1],TEAL,5)
    p.text(1408,xr['D'][1]-25,'Sequence loss',16,bold=True,fill=TEAL)
    p.lines(1259,617,['Infer where each character changed.',
        'Estimate ancestral-state probabilities and gain/loss rates*.',
        '* Optional likelihood analysis, when estimable.'],15,step=28)
    p.rect(48,722,1704,49,fill=PALE,stroke='none',r=6)
    p.text(67,753,'Local characters, not complete ancestral transcripts. A gene may contain multiple changes; mutation mechanisms are not inferred.',18)
    return p.write(path)
