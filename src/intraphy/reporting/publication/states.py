"""Schematic definitions of local configurations, repertoires and tip compatibility."""
from .drawing import box, line, text, track, document, PALE, RULE, TEAL, PURPLE, GOLD, MUTED


def draw():
    b = [text(42, 46, "Local structural states and observation constraints", 27, bold=True),
         text(42, 77, "Ordered local exon structures on shared genomic material", 15, color=MUTED)]
    # Configurations represent complete ordered structures on material m.
    b += [box(45, 120, 500, 280, PALE), text(68, 158, "A   Exon configuration C", 19, bold=True),
          text(75, 217, "C₁", 16, bold=True)]
    b += track(126, 211, 340, 100, [(4,31),(49,72),(78,96)], height=27)
    b += [text(75, 282, "C₂", 16, bold=True)]
    b += track(126, 276, 340, 100, [(4,31),(49,96)], color=PURPLE, height=27)
    b += [text(68, 338, "[a b]  ↔  [a]—[b]", 18), text(68, 365, "Ordered exon intervals", 13, color=MUTED)]
    # Coexistence versus a set of alternatives.
    b += [box(582, 120, 520, 280, "white", RULE), text(605, 158, "B   Local repertoire R", 19, bold=True),
          text(607, 196, "R = {C₁, C₂} on m", 14, color=MUTED),
          text(608, 248, "AND", 17, color=TEAL, bold=True), text(695, 248, "coexisting structures", 15),
          text(695, 280, "Distinct members", 13, color=MUTED),
          text(608, 316, "OR", 17, color=GOLD, bold=True), text(695, 316, "complete alternatives", 13),
          text(608, 351, "Observation sets a lower bound", 13, color=MUTED),
          text(608, 384, "additional members may occur", 13, color=MUTED)]
    # Unknown tip is a compatibility vector of ones, not a uniform prior.
    b += [box(1140, 120, 602, 280, PALE), text(1163, 158, "C   Tip compatibility", 19, bold=True),
          text(1163, 198, "Allowed states", 14, color=MUTED)]
    for j, (x,label) in enumerate(((1190,"C₁"),(1315,"C₂"),(1440,"C₃"),(1565,"…"))):
        b += [f'<circle cx="{x}" cy="252" r="25" fill="white" stroke="{TEAL}" stroke-width="2"/>',
              text(x,258,label,16,anchor="middle",bold=True), text(x,309,"1",15,anchor="middle",color=TEAL,bold=True)]
    b += [text(1163, 342, "Fully unknown → all weights 1", 13),
          text(1163, 378, "Known DNA presence retained", 13, color=MUTED)]
    # Catalogue closure emphasizes legal intermediate states.
    b += [line(42, 411, 1742, 411, RULE, 1), text(48, 453, "D   Candidate exon boundaries", 18, bold=True),
          text(48, 482, "observed + supported candidates", 14, color=MUTED),
          line(545, 465, 640, 465, arrow=True), text(670,453,"admissible closure",18,bold=True),
          text(670,482,"reachable states",14,color=MUTED)]
    for x,y,label in ((1015,450,"C₁"),(1155,450,"C₂"),(1085,520,"C₃"),(1235,520,"…")):
        b += [f'<circle cx="{x}" cy="{y}" r="20" fill="{PALE}" stroke="{PURPLE}" stroke-width="2"/>',
              text(x,y+6,label,14,anchor="middle",bold=True)]
    b += [line(1035,450,1135,450,arrow=True),line(1168,466,1100,504,arrow=True),
          line(1105,520,1215,520,arrow=True),
          text(1320, 468, "Finite state space", 15, bold=True),
          text(1320, 506, "Conditional on catalogue", 15, bold=True)]
    return document(1790, 570, "Local structural states and observation constraints", b,
                    "Panel A defines an ordered exon configuration on genomic material. Panel B defines a local structural repertoire: coexisting configurations map injectively to distinct members, while tip observations constrain a lower bound and allow additional members. Panel C shows compatibility weights; fully unknown observations have weight one for every admissible state, with known material-presence constraints retained. Panel D closes the finite state space over candidate exon boundaries and admissible edits.")
