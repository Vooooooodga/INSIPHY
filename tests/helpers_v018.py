"""Small, explicit structural matrices for v18 unit tests."""
from pathlib import Path
from intraphy.observations.schema import write_structural_site_matrix
from intraphy.topology import SpeciesTree


def observation(site, species, state, *, layer="splice_junction", group="NA"):
    zero, one = (("not_exonic", "exonic") if layer == "exon_role" else ("absent", "present"))
    value = "unknown" if state == "unknown" else (one if state else zero)
    return dict(family_id="family", layer=layer, site_id=site, species=species,
                state=value, state_0=zero, state_1=one, annotation_view="repertoire",
                applicability="applicable", observation_mask="missing" if state == "unknown" else "observed",
                observation_reason="test_missing" if state == "unknown" else "test_observation",
                evidence="explicit_test_input", site_kind="within_element_junction" if layer == "splice_junction" else "sequence_unit",
                linked_group_id=group, discovery_rule="independent_catalogue")


def write_case(directory, rows):
    root=Path(directory); source=root/"input"; source.mkdir(parents=True)
    source.joinpath("species_tree.tsv").write_text(
        "node_id\tparent_id\tlabel\tbranch_length\nroot\t\troot\t0\n"
        "a\troot\tA\t1\nb\troot\tB\t1\nc\troot\tC\t1\n")
    matrix=source/"matrix.tsv"; write_structural_site_matrix(matrix, rows)
    return source, root/"output", matrix


def small_tree():
    return SpeciesTree([
        dict(node_id="root",parent_id="",label="root",branch_length=0),
        dict(node_id="a",parent_id="root",label="A",branch_length=1),
        dict(node_id="inner",parent_id="root",label="inner",branch_length=1),
        dict(node_id="b",parent_id="inner",label="B",branch_length=1),
        dict(node_id="c",parent_id="inner",label="C",branch_length=1)])
