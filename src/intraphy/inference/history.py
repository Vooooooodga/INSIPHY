"""A reproducible, globally compatible minimum-change reconstruction."""
from ..observations.characters import character_id


def optimal_history(tree, observations, minimum, node_states, branch_pairs):
    """Select one witness of the optimum; ties do not define probabilities.

    Endpoint pairs come from the inside/outside Sankoff recursion. Conditional
    choices on a tree are compatible; the final cost is independently checked.
    """
    states = {tree.root: min(node_states[tree.root])}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            allowed = {b for a, b in branch_pairs[(parent, child)] if a == states[parent]}
            if not allowed:
                raise ValueError("No compatible child state in optimal reconstruction")
            states[child] = states[parent] if states[parent] in allowed else min(allowed)
    cost = sum(states[parent] != states[child] for parent, child in tree.edges())
    if cost != minimum:
        raise ValueError("Reconstructed history does not attain the global minimum")
    for label, node in tree.leaf_by_label.items():
        observed = observations.get(label, "unknown")
        if observed in {0, 1} and states[node] != observed:
            raise ValueError("Reconstructed history conflicts with an observed tip")
    return states


def history_rows(key, tree, observations, minimum, node_states, branch_pairs):
    if all(value == "unknown" for value in observations.values()):
        return []
    states = optimal_history(tree, observations, minimum, node_states, branch_pairs)
    return [{"family_id": key[0], "layer": key[1], "site_id": key[2],
             "character_id": character_id(key), "node_id": node,
             "node_label": tree.label[node], "state_index": states[node],
             "history_id": "one_minimum_change_reconstruction",
             "selection_rule": "lowest_root_state_then_parent_state_preferred",
             "probability": "NA", "minimum_changes": int(minimum)}
            for node in tree.preorder()]
