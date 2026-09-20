"""Explicit adapters for older exported results, never used by numeric kernels."""

def posterior_row(row):
    result = dict(row)
    if "posterior_available" not in result:
        # Compatibility with pre-0.15 exports; current rows have a typed flag.
        tokens = set(str(result.get("conditioning", "")).split(";"))
        result["posterior_available"] = (
            "fit_status=success" in tokens and
            result.get("structural_change_type") != "posterior_not_reported"
        )
    return result
