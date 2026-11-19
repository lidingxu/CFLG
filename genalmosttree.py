"""Generate almost-tree networks (tree_A: Pruefer-random labeled tree, tree_B: random BFS tree) plus extra edges.

Every instance uses its own random.Random(seed) for the tree, the extra edges and the edge lengths, so the
output is reproducible. (The earlier tree_A instances drew their extra edges and edge lengths from Python's
unseeded global random stream and could not be regenerated; they were replaced by this set.)

Default: the benchmark tree sets (one instance per generator x n x p, instance seed 1), written to
benchmarks/tree_{A,B}/<prefix>.<n>.<extra>.s<seed>
    python genalmosttree.py
    python genalmosttree.py --ns 100 200 300 --ps 0.3 0.4 0.5 --seeds 1 2 3 --outdir other_dir
"""
import argparse
import os
import random

import networkx as nx


def add_extra_edges(tree, extra_edges, rng):
    # add extra edges (creating cycles) until extra_edges were added or too many attempts failed
    nodes = list(tree.nodes())
    added = 0
    attempts = 0
    max_attempts = 10 * extra_edges  # avoid infinite loops
    while added < extra_edges and attempts < max_attempts:
        u, v = rng.sample(nodes, 2)
        if not tree.has_edge(u, v):
            tree.add_edge(u, v)
            added += 1
        attempts += 1
    return tree


def generate_almost_tree(n, extra_edges, rng):
    # random labeled spanning tree (uniform via Pruefer sequences) plus extra edges
    tree = nx.random_labeled_tree(n, seed=rng)
    return add_extra_edges(tree, extra_edges, rng)


def random_bfs_tree(n, extra_edges, rng):
    # BFS tree where each parent gets 1-3 children, plus extra edges
    available = list(range(n))
    root = available.pop(rng.randrange(len(available)))
    tree = nx.Graph()
    tree.add_node(root)
    queue = [root]
    while available and queue:
        parent = queue.pop(0)
        num_children = rng.randint(1, min(3, len(available)))  # small branching
        children = rng.sample(available, num_children)
        for child in children:
            tree.add_edge(parent, child)
            queue.append(child)
            available.remove(child)
    return add_extra_edges(tree, extra_edges, rng)


def graph_to_string(G, rng):
    lines = [f"{G.number_of_nodes()} {G.number_of_edges()}"]
    for u, v in G.edges():
        length = rng.uniform(1, 2)  # edge lengths uniform in [1, 2]
        lines.append(f"{u + 1} {v + 1} {length}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ns", type=int, nargs="+", default=[100, 150, 200, 250, 300, 350])
    parser.add_argument("--ps", type=float, nargs="+", default=[0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--outdir", default="benchmarks")
    args = parser.parse_args()

    count = 0
    for prefix, method in (("treeA", generate_almost_tree), ("treeB", random_bfs_tree)):
        outsubdir = os.path.join(args.outdir, "tree_" + prefix[-1])
        os.makedirs(outsubdir, exist_ok=True)
        for n in args.ns:
            for p in args.ps:
                extra = int(n * p)
                for s in args.seeds:
                    rng = random.Random(f"{prefix}-{n}-{p}-{s}")  # one deterministic stream per instance
                    G = method(n, extra, rng)
                    assert nx.is_connected(G)
                    with open(os.path.join(outsubdir, f"{prefix}.{n}.{extra}.s{s}"), "w") as f:
                        f.write(graph_to_string(G, rng))
                    count += 1
    print(f"wrote {count} instances to {args.outdir}")


if __name__ == "__main__":
    main()
