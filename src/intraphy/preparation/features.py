"""Indexes over annotation records. These indexes never decide biological states."""
from collections import defaultdict


class FeatureHierarchy:
    """Preserve all records for each ID and Parent; do not collapse isoforms."""

    def __init__(self, rows):
        self.rows = tuple(rows)
        self.by_id = defaultdict(list)
        self.children = defaultdict(list)
        for row in self.rows:
            identifier = row.get("id")
            if identifier not in {None, "", "NA", ".", "unknown"}:
                self.by_id[identifier].append(row)
            for parent in self.parents(row):
                self.children[parent].append(row)

    @staticmethod
    def parents(row):
        # Exported raw rows use parent; parsed annotations also provide parents.
        values = row.get("parents")
        if not isinstance(values, (list, tuple, set, frozenset)):
            values = str(row.get("parent") or "").split(",")
        return tuple(dict.fromkeys(str(x).strip() for x in values
                                   if str(x).strip() not in {"", "NA", ".", "unknown"}))

    def transcript_features(self, transcript_id):
        """Return descendants and ancestor attributes in the historical traversal order."""
        selected, visited = [], set()
        pending = [transcript_id]
        while pending:
            identifier = pending.pop()
            if identifier in {None, "", "NA", ".", "unknown"} or identifier in visited:
                continue
            visited.add(identifier)
            selected.extend(self.by_id.get(identifier, ()))
            for row in self.children.get(identifier, ()):
                child_id = row.get("id")
                if child_id not in {None, "", "NA", ".", "unknown"}:
                    pending.append(child_id)
                else:
                    selected.append(row)
        pending = list(self.by_id.get(transcript_id, ()))
        while pending:
            row = pending.pop()
            for parent in self.parents(row):
                if parent not in visited:
                    visited.add(parent)
                    ancestors = self.by_id.get(parent, ())
                    selected.extend(ancestors)
                    pending.extend(ancestors)
        return selected
