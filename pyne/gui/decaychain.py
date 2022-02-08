"""Tools for generating decay chain figures for a nuclide. This requires the
python-graphviz library.
"""

from graphviz import Digraph

from pyne import data
from pyne import nucname
from pyne import rxname


PRETTY_RX = {
    "a": "α",
    "alpha": "α",
    "bplus": "β+",
    "bplus_p": "β+p",
    "bplus_n": "β+n",
    "bminus": "β-",
    "bminus_p": "β-p",
    "bminus_n": "β-n",
}


def graph(nuc):
    """Returns a graphviz Digraph object for the decay chain starting from nuc."""
    i = nucname.name(nuc)
    name = nucname.name(nuc)
    dot = Digraph(comment="Decay chain for " + nuc)
    dot.node(name)
    nodes_seen = {name}
    kids = data.decay_children(i)
    from_to = {(i, j) for j in kids}
    edges_seen = set()
    while len(from_to) != 0:
        new_from_to = set()
        for ft in from_to:
            i, j = ft
            jname = nucname.name(j)
            if jname not in nodes_seen:
                dot.node(jname)
                nodes_seen.add(jname)
            if ft not in edges_seen:
                iname = nucname.name(i)
                try:
                    label = rxname.name(i, j, "decay")
                except RuntimeError:
                    label = "sf"
                label = PRETTY_RX.get(label, label)
                br = data.branch_ratio(i, j)
                if br < 1.0:
                    label += ", br={0:.3}".format(br)
                dot.edge(iname, jname, label=label)
                edges_seen.add(ft)
            kids = data.decay_children(j)
            new_from_to |= {(j, k) for k in kids}
        from_to = new_from_to
    return dot
