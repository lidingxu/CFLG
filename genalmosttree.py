import networkx as nx
import matplotlib.pyplot as plt
import random
import os

def generate_almost_tree(n, extra_edges, seed=None):
    # Step 1: Generate a random spanning tree with n nodes
    tree = nx.random_labeled_tree(n, seed=seed)

    # Step 2: Add a few random extra edges to make it "almost a tree"
    nodes = list(tree.nodes())
    added = 0
    attempts = 0
    max_attempts = 10 * extra_edges  # avoid infinite loops

    while added < extra_edges and attempts < max_attempts:
        u, v = random.sample(nodes, 2)
        if not tree.has_edge(u, v) and not nx.has_path(tree, u, v):
            # avoid adding an edge that makes it a tree again
            tree.add_edge(u, v)
            added += 1
        elif not tree.has_edge(u, v):
            # we allow adding edges that create cycles
            tree.add_edge(u, v)
            added += 1
        attempts += 1

    return tree

def random_bfs_tree(n, extra_edges, seed=None):
    random.seed(seed)
    G = nx.empty_graph(n)
    available = list(range(n))
    root = available.pop(random.randrange(len(available)))
    tree = nx.Graph()
    tree.add_node(root)
    queue = [root]
    while available and queue:
        parent = queue.pop(0)
        num_children = random.randint(1, min(3, len(available)))  # small branching
        children = random.sample(available, num_children)
        for child in children:
            tree.add_edge(parent, child)
            queue.append(child)
            available.remove(child)

    # Step 2: Add a few random extra edges to make it "almost a tree"
    nodes = list(tree.nodes())
    added = 0
    attempts = 0
    max_attempts = 10 * extra_edges  # avoid infinite loops

    while added < extra_edges and attempts < max_attempts:
        u, v = random.sample(nodes, 2)
        if not tree.has_edge(u, v) and not nx.has_path(tree, u, v):
            # avoid adding an edge that makes it a tree again
            tree.add_edge(u, v)
            added += 1
        elif not tree.has_edge(u, v):
            # we allow adding edges that create cycles
            tree.add_edge(u, v)
            added += 1
        attempts += 1

    return tree

def graph_to_string(G):
    lines = []
    lines.append(f"{G.number_of_nodes()} {G.number_of_edges()}")
    for u, v in G.edges():
        # If you want to assign a random length, for example between 1 and 10:
        length = random.uniform(1, 2)
        lines.append(f"{u + 1} {v + 1} {length}")
    return "\n".join(lines)
ns = [100, 150, 200, 250, 300, 350]
ps = [0.1,0.15,0.2]
# Example usage
outdir = "benchmarks"
os.makedirs(outdir, exist_ok=True)
i = 0
for method in [generate_almost_tree, random_bfs_tree]:
    for n in ns:
        for p in ps:
            extra = int(n * p)
            G = method(n, extra, seed=42)
            gstring = graph_to_string(G)
            subdir = f"tree_{'A' if i == 0 else 'B'}"
            subdir_ = f"tree{'A' if i == 0 else 'B'}"
            outsubdir = os.path.join(outdir, subdir)
            os.makedirs(outsubdir, exist_ok=True)
            filename = f"{subdir_}.{n}.{extra}"
            filepath = os.path.join(outsubdir, filename)
            with open(filepath, "w") as f:
                f.write(gstring)

    i += 1
