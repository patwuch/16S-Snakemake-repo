"""Name a tree node by the rank its classifier actually reached.

write_res.R fills an empty rank with the nearest named ancestor plus "_sp" (for
species) or "_unresolved" (any other rank), so that each unresolved lineage
keys its own node: "g__COE1.s__" becomes "g__COE1.s__COE1_sp". Printed as-is,
that tip reads as a species. `resolve(lineage)` walks the lineage from the root
with the same rule, and returns the deepest real token plus the ranks that were
filled in, so a plotter can draw "g__COE1 (no species call)" instead.
"""

RANK = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class", "o": "order",
        "f": "family", "g": "genus", "s": "species"}


def resolve(lineage):
    """(deepest real token, missing rank words) for a dot-separated lineage."""
    tokens = [t for t in lineage.split(".") if t and t != "root"]
    last_named, real, missing = None, None, []
    for tok in tokens:
        prefix, name = tok[:1], tok[3:] if tok[1:3] == "__" else tok
        tag = "sp" if prefix == "s" else "unresolved"
        if last_named is not None and name == f"{last_named}_{tag}":
            missing.append(RANK.get(prefix, prefix))
        else:
            last_named, real, missing = name, tok, []
    return (real or (tokens[-1] if tokens else lineage)), missing


def suffix(missing):
    """' (no genus or species call)' for the ranks resolve() reports missing."""
    if not missing:
        return ""
    words = missing[0] if len(missing) == 1 else ", ".join(missing[:-1]) + " or " + missing[-1]
    return f" (no {words} call)"
