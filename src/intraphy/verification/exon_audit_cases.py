"""Audit regressions: observation damage is distinct from a structural mutation."""
from __future__ import annotations
from pathlib import Path
import json
from .exon_cases import write_exon_example
from ..storage.fasta import parse_fasta

AUDIT_SCENARIOS = ("mis_split_unchanged_DNA", "mis_fusion_unchanged_DNA",
    "one_boundary_unchanged_DNA", "both_boundaries_unchanged_DNA",
    "first_exon_dropout", "last_exon_dropout", "split_with_changed_signals",
    "shift_with_changed_signal", "exact_exon_deletion")


def write_audit_example(output_dir, name):
    if name not in AUDIT_SCENARIOS: raise ValueError("Unknown structural audit scenario")
    template = {"mis_split_unchanged_DNA": "intronization", "one_boundary_unchanged_DNA": "donor_shift",
                "split_with_changed_signals": "intronization", "shift_with_changed_signal": "donor_shift",
                "mis_fusion_unchanged_DNA": "utr"}.get(name, "conserved")
    raw = write_exon_example(output_dir, template)
    file = raw/"Species_D.gff3"; lines = file.read_text().splitlines()
    fasta = raw/"Species_D.fa"; records = parse_fasta(fasta); seq = next(iter(records.values())); contig = next(iter(records))
    if name in {"first_exon_dropout", "last_exon_dropout"}:
        n = 1 if name.startswith("first") else 4
        lines = [l for l in lines if f"_e{n};" not in l and f"_c{n};" not in l]
    elif name == "both_boundaries_unchanged_DNA":
        changed = []
        for line in lines:
            if "_e2;" in line or "_c2;" in line:
                x = line.split("\t");x[3]=str(int(x[3])+3);x[4]=str(int(x[4])-3);line="\t".join(x)
            changed.append(line)
        lines = changed
    elif name == "mis_fusion_unchanged_DNA":
        changed=[]
        for line in lines:
            if "_e2;" in line: continue
            if "_e1;" in line:
                x=line.split("\t");x[4]="530";line="\t".join(x)
            changed.append(line)
        lines=changed
    elif name == "split_with_changed_signals":
        seq = seq[:410]+"GT"+seq[412:438]+"AG"+seq[440:]
    elif name == "shift_with_changed_signal":
        seq = seq[:527]+"GT"+seq[529:]
    elif name == "exact_exon_deletion":
        left,right=350,530;seq=seq[:left]+seq[right:];changed=[]
        for line in lines:
            if line.startswith("#"): changed.append(line);continue
            x=line.split("\t");a,b=int(x[3])-1,int(x[4])
            if a>=left and b<=right and x[2] in {"exon","CDS"}: continue
            if a>=right: a-=right-left
            elif a>=left: a=left
            if b>=right: b-=right-left
            elif b>left: b=left
            if a>=b: continue
            x[3],x[4]=str(a+1),str(b);changed.append("\t".join(x))
        lines=changed
    file.write_text("\n".join(lines)+"\n");fasta.write_text(f">{contig}\n{seq}\n")
    truth = json.loads((raw/"truth.json").read_text())
    truth.update(audit_scenario=name, annotation_only=(name.startswith("mis_") or "unchanged_DNA" in name or "dropout" in name))
    (raw/"truth.json").write_text(json.dumps(truth,indent=2)+"\n")
    return raw


def check_audit_result(name, summary, details, result_dir):
    total=sum(r.get("minimum_structural_edits") or 0 for r in summary)
    annotation=sum(min((h.get("minimum_structural_edits") for h in u.get("views",{}).get("annotation",{}).get("histories",())
                       if h.get("minimum_structural_edits") is not None),default=0) for u in details["units"])
    events=[e for u in details["units"] for e in u.get("views",{}).get("evidence",{}).get("events",()) if e["support"]=="required"]
    genomic_mutation=name in {"split_with_changed_signals","shift_with_changed_signal","exact_exon_deletion"}
    if genomic_mutation:
        assert total == 1 and len(events)==1,(name,total,events)
    else:
        assert total==0 and not events,(name,total,events)
        assert annotation>0,(name,"original annotation history was lost")
    if name.endswith("dropout"):
        rows=parse_fasta(Path(result_dir)/"alignment_evidence/family_00001/genomic_loci.fa")
        assert {len(v) for v in rows.values()} == {1140},(name,"available locus DNA was cropped")
    return {"minimum_evidence_edits":total,"minimum_annotation_edits":annotation,
            "required_evidence_edits":len(events),"annotation_only":not genomic_mutation}
