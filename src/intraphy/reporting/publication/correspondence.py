"""Sequence correspondence and observation qualification schematic."""
from .drawing import box, line, text, track, document, RULE, TEAL, GOLD, PURPLE, RED, MUTED


def draw():
    b = [text(42, 46, "Sequence correspondence and structural observations", 27, bold=True),
         text(42, 76, "A shared coordinate system is conditional on qualified sequence evidence", 15, color=MUTED)]
    # Genome and coding evidence converge on ordered interval correspondence.
    b += [text(55, 131, "A   Genomic locus", 19, bold=True), text(55, 160, "annotated exons", 14, color=MUTED)]
    b += track(62, 217, 500, 100, [(4,24),(31,45),(69,80),(87,98)], height=24)
    b += [text(55, 190, "genomic DNA", 13, color=MUTED),
          text(55, 387, "coding alignment → genome", 14, color=MUTED)]
    b += track(62, 346, 500, 100, [(4,24),(31,45),(69,80),(87,98)], color=PURPLE, height=15)
    b += [line(182, 217, 182, 345, GOLD, 2), line(302,217,302,345,GOLD,2),
          line(412,217,412,345,GOLD,2), line(507,217,507,345,GOLD,2),
          text(61, 425, "collinear sequence blocks", 13, color=MUTED),
          line(580, 280, 632, 280, arrow=True)]
    # Symbolic candidate DAG: no scores or unsupported probability labels.
    b += [text(650, 131, "B   Ordered candidate chain", 19, bold=True),
          text(650, 159, "compatible order and strand", 14, color=MUTED)]
    nodes = ((670,220,"Flank"),(834,220,"Path A"),(834,302,"Path B"),(1000,260,"Flank"))
    for x,y,label in nodes:
        b += [box(x,y,126,42,"white",TEAL,5), text(x+63,y+29,label,13,anchor="middle")]
    b += [line(796,241,834,241,arrow=True), line(796,241,834,323,TEAL,1.5,arrow=True),
          line(960,241,1000,281,arrow=True), line(960,323,1000,281,arrow=True),
          text(650, 379, "Alternative compatible chains", 13, color=MUTED),
          line(1138, 280, 1180, 280, arrow=True)]
    # Three distinct evidence outcomes.
    b += [text(1190, 131, "C   Tip observations", 19, bold=True)]
    rows = ((171, ("DNA present", "exon status unknown"), TEAL, "?"),
            (268, ("DNA absent", "paired flanks qualify"), GOLD, "∅"),
            (365, ("copy / strand / alignment", "conflict unresolved"), RED, "!"))
    for y,labels,color,mark in rows:
        b += [box(1190,y,542,88,"white",RULE,5), text(1210,y+56,mark,23,color=color,bold=True),
              text(1250,y+39,labels[0],14), text(1250,y+75,labels[1],14)]
    b += [line(42, 473, 1732, 473, RULE, 1),
          text(48, 514, "Compatible alternatives retain complete exon configurations.", 16),
          text(48, 544, "Correspondence is conditional on sequence qualification.", 16, color=MUTED)]
    return document(1780, 590, "Sequence correspondence and structural observations", b,
                    "Genomic and coding-projection evidence support ordered correspondence hypotheses. DNA present with unknown exon status, supported DNA absence, and unresolved conflicts are shown as distinct outcomes.")
