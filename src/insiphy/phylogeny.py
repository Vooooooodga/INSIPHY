"""Public inference dispatch; formal engines do not load legacy models."""

def infer_phylogeny(
    input_dir,
    output_dir,
    bootstrap_replicates=0,
    stochastic_maps=0,
    seed=7,
    foreground_branches=None,
    analysis_scope="single-copy",
    model="parsimony",
    branch_length_mode="supplied",
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
    structural_site_matrix_path=None,
    annotation_view="repertoire",
    analysis_range="all",
    min_callable_fraction=0.70,
):
    if analysis_scope == "single-copy":
        if model not in {"parsimony", "er-ard", "foreground"}:
            raise ValueError(f"unsupported single-copy model: {model!r}")
        if bootstrap_replicates or stochastic_maps:
            raise SystemExit(
                "single-copy analysis uses parsimony or analytic likelihood; "
                "--bootstrap-replicates and --stochastic-maps must be 0"
            )
        if model != "foreground" and foreground_branches:
            raise SystemExit("--foreground-branches requires --model foreground")
        if int(threads) < 1:
            raise SystemExit("--threads must be at least 1")
        from insiphy.run_result import begin_run
        begin_run(output_dir, model, analysis_scope)
        from insiphy.observations.matrix import prepare_observation_matrix
        matrix = prepare_observation_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                                            analysis_range=analysis_range, min_callable_fraction=min_callable_fraction)
        if model == "parsimony":
            from insiphy.parsimony import infer_single_copy_parsimony

            return infer_single_copy_parsimony(
                input_dir,
                output_dir,
                threads=threads,
                structural_site_matrix_path=structural_site_matrix_path,
                annotation_view=annotation_view,
                observation_matrix=matrix,
            )
        from insiphy.structural_phylogeny import infer_single_copy_phylogeny

        return infer_single_copy_phylogeny(
            input_dir,
            output_dir,
            model=model,
            foreground_branches=foreground_branches,
            branch_length_mode=branch_length_mode,
            ascertainment=ascertainment,
            threads=threads,
            root_frequency=root_frequency,
            root_presence=root_presence,
            structural_site_matrix_path=structural_site_matrix_path,
            annotation_view=annotation_view,
            observation_matrix=matrix,
        )
    if analysis_scope != "experimental-multicopy":
        raise ValueError(f"unsupported analysis_scope: {analysis_scope!r}")
    if analysis_range != "all" or min_callable_fraction != 0.70:
        raise ValueError("callable-scope policy is a formal single-copy feature")
    from insiphy.experimental.phylogeny import infer_experimental_phylogeny
    from insiphy.run_result import record_run_result
    result = infer_experimental_phylogeny(
        input_dir, output_dir, bootstrap_replicates=bootstrap_replicates,
        stochastic_maps=stochastic_maps, seed=seed, foreground_branches=foreground_branches,
        analysis_scope=analysis_scope, model=model, branch_length_mode=branch_length_mode,
        ascertainment=ascertainment, threads=threads, root_frequency=root_frequency,
        root_presence=root_presence, structural_site_matrix_path=structural_site_matrix_path,
        annotation_view=annotation_view,
    )
    record_run_result(output_dir, input_dir, "experimental-multicopy", analysis_scope)
    return result

_LEGACY_NAMES = frozenset({'call_scope_for_pattern', 'read_structural_tree', 'interpretation_hints', 'interpretation_caveat_for_scope', 'interpretation_hint_for_pattern', 'best_branch_support', 'tree_depths', 'mrca_node', 'make_event', 'element_phylogenetic_coverage', 'classify_branch_event', 'tree_tip_label_for_occ', 'role_from_occurrences', 'state_from_occurrences', 'add_q_values', 'copy_tip_label', 'event_support_summary', 'ancestor_set', 'copy_tip_label_from_adjacency', 'add_character', 'informative_source_labels', 'internal_homology_phylogenetic_coverage', 'completion_event_class', 'read_foreground_edges', 'fallback_element_rows', 'has_multicopy_species', 'tree_tip_label_for_adjacency', 'phylogenetic_coverage'})

def __getattr__(name):
    if name in _LEGACY_NAMES:
        from insiphy.experimental import phylogeny
        return getattr(phylogeny, name)
    raise AttributeError(name)
