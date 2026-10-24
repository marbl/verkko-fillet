#!/usr/bin/env python3
"""
Simplified reimplementation of neighborhood segment extraction from https://github.com/snurk/gfacpp.
Extract segment names within -r hops of seed nodes in a GFA graph.
Prints matching segment names to stdout, one per line.
"""

import argparse
import sys
from collections import defaultdict, deque


def parse_gfa(path):
    """
    Read a GFA file and return:
      segments: set of segment names
      adj: dict mapping segment -> set of neighbor segments (undirected)
    Only S and L lines are consulted; orientations are collapsed since
    neighborhood expansion treats the graph as undirected.
    """
    segments = set()
    adj = defaultdict(set)

    with open(path) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            fields = line.rstrip('\n').split('\t')
            rec = fields[0]
            if rec == 'S':
                # S <name> <sequence> [tags...]
                if len(fields) >= 2:
                    segments.add(fields[1])
            elif rec == 'L':
                # L <from> <from_orient> <to> <to_orient> <overlap> [tags...]
                if len(fields) >= 5:
                    a, b = fields[1], fields[3]
                    adj[a].add(b)
                    adj[b].add(a)
    return segments, adj


def read_seeds(path):
    seeds = set()
    with open(path) as f:
        for line in f:
            name = line.strip()
            if name and not name.startswith('#'):
                seeds.add(name)
    return seeds


def collect_neighborhood(adj, seeds, radius):
    """BFS out to `radius` hops from any seed. Returns set of segment names."""
    visited = {}  # name -> shortest hop distance
    queue = deque()
    for s in seeds:
        visited[s] = 0
        queue.append((s, 0))

    while queue:
        node, depth = queue.popleft()
        if depth == radius:
            continue
        for nb in adj.get(node, ()):
            nd = depth + 1
            if nb not in visited or visited[nb] > nd:
                visited[nb] = nd
                queue.append((nb, nd))

    return visited


def main():
    ap = argparse.ArgumentParser(
        description="Print GFA segments within -r hops of seed nodes."
    )
    ap.add_argument('gfa', help='input GFA file')
    ap.add_argument('-n', '--nodes', required=True,
                    help='file listing seed segment names (one per line)')
    ap.add_argument('-r', '--radius', type=int, default=10,
                    help='neighborhood radius in hops (default: 10)')
    args = ap.parse_args()

    segments, adj = parse_gfa(args.gfa)
    seeds = read_seeds(args.nodes)

    missing = seeds - segments
    if missing:
        print(f"warning: {len(missing)} seed(s) not found in graph, e.g. "
              f"{sorted(missing)[:3]}", file=sys.stderr)

    neighborhood = collect_neighborhood(adj, seeds, args.radius)

    for name in sorted(neighborhood):
        print(name)


if __name__ == '__main__':
    main()