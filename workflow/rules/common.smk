# Helpers for summarizing BUSTED-S-MH JSON output.

def _legacy_branch_count(branches):
    """Count branches in older HyPhy list-based branch sections."""
    if isinstance(branches, dict):
        counts = [_legacy_branch_count(value) for value in branches.values()]
        known_counts = [count for count in counts if count is not None]
        return sum(known_counts) if known_counts else None
    if isinstance(branches, list):
        return len(branches)
    if isinstance(branches, int) and not isinstance(branches, bool):
        return branches
    return None


def _labeled_branch_counts(tested):
    """Read HyPhy branch-name -> branch-set mappings under ``tested``."""
    if not isinstance(tested, dict):
        return None

    labels = []
    for partition in tested.values():
        if not isinstance(partition, dict):
            continue
        partition_labels = list(partition.values())
        if partition_labels and all(isinstance(label, str) for label in partition_labels):
            normalized = [label.strip().lower() for label in partition_labels]
            if all(label in {"test", "background"} for label in normalized):
                labels.extend(normalized)

    if not labels:
        return None

    return labels.count("test"), labels.count("background")


def count_busteds_mh_branches(data):
    """Return ``(tested, background)`` branch counts across HyPhy formats."""
    tested = data.get("tested")

    labeled_counts = _labeled_branch_counts(tested)
    if labeled_counts is not None:
        return labeled_counts

    return (
        _legacy_branch_count(tested),
        _legacy_branch_count(data.get("background")),
    )