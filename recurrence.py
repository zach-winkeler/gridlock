import itertools
from copy import copy
from functools import cache
from itertools import combinations


# represents a polynomial of one variable called "k"
class Poly:
    def __init__(self, cs):
        trailhead = len(cs)
        for i in reversed(list(range(len(cs)))):
            if cs[i] == 0:
                trailhead = i
            else:
                break
        self.cs = cs[:trailhead]

    def __add__(self, other):
        return Poly([a + b for a, b in itertools.zip_longest(self.cs, other.cs, fillvalue=0)])

    def __sub__(self, other):
        return Poly([a - b for a, b in itertools.zip_longest(self.cs, other.cs, fillvalue=0)])

    def __mul__(self, other):
        if isinstance(other, int):
            return Poly([other * c for c in self.cs])
        else:
            d = len(self.cs) - 1 + len(other.cs) - 1
            zipped = list(itertools.zip_longest(self.cs, other.cs, [0] * (d + 1), fillvalue=0))
            return Poly([sum(zipped[i][0] * zipped[n - i][1] for i in range(n + 1)) for n in range(d + 1)])

    def __rmul__(self, other):
        return self * other

    def __repr__(self):
        return " + ".join([f"{c}k^{i}" for i, c in enumerate(self.cs)])

    def __hash__(self):
        return hash(self.cs)

    def __eq__(self, other):
        return self.cs == other.cs

    def evaluate(self, v):
        return sum(c * (v ** i) for i, c in enumerate(self.cs))

    def degree(self):
        return max(i for i, c in enumerate(self.cs) if c != 0)


zero = Poly([])


def k(power):
    return Poly([0, ] * power + [1])


def initialize_from_set(elems):
    return Partition({e: i for i, e in enumerate(elems)},
                     {i: e for i, e in enumerate(elems)},
                     [i for i in range(len(elems))])


def group_by_key(elems, key):
    groups = {}
    for e in elems:
        k = key(e)
        if k in groups:
            groups[k].append(e)
        else:
            groups[k] = [e]
    return sorted(groups.values(), key=len, reverse=True)


# represents a set partition
# elements of the set are first indexed, and then their indices are stored as parts
class Partition:
    def __init__(self, to_num, from_num, parts):
        self.to_num = to_num
        self.from_num = from_num
        self.parts = parts

    def copy(self):
        return Partition(copy(self.to_num), copy(self.from_num), copy(self.parts))

    def join_two(self, e1, e2):
        out = self.copy()
        i1, i2 = self.to_num[e1], self.to_num[e2]
        p1, p2 = self.parts[i1], self.parts[i2]
        if p1 < p2:
            for i, p in enumerate(self.parts):
                if p == p2:
                    out.parts[i] = p1
        elif p2 < p1:
            for i, p in enumerate(self.parts):
                if p == p1:
                    out.parts[i] = p2

        return out

    def join_many(self, es):
        out = self.copy()
        for e in es[1:]:
            out = out.join_two(es[0], e)
        return out

    def find(self, e):
        return self.parts[self.to_num[e]]

    def __len__(self):
        return len(set(self.parts))

    def __hash__(self):
        return sum(p * (29 ** i) for i, p in enumerate(self.parts))

    def __eq__(self, other):
        return all(p == other.parts[i] for i, p in enumerate(self.parts))


def lo(g, strategy=lambda g: sorted(g.nodes, key=g.degree, reverse=True)):
    # constructs initial vertex partition
    # finds leaves and bivalent vertices to save time
    def initial_values():
        partition = initialize_from_set(list(g.nodes))

        # bivalent vertices share the same color as their neighbors
        bivalent_vertices = [v for v in g.nodes if g.degree(v) == 2]
        for v in bivalent_vertices:
            e1, e2 = g.edges(v)
            partition = partition.join_two(e1[0], e1[1])
            partition = partition.join_two(e2[0], e2[1])

        # leaves share the same color as their neighbor
        leaves = [v for v in g.nodes if g.degree(v) == 1]
        for v in leaves:
            (e,) = g.edges(v)
            partition = partition.join_two(e[0], e[1])

        return partition, {v for v in g.nodes
                           if g.degree(v) == 0
                           or sum(partition.find(v) == partition.find(w) for w in g.neighbors(v)) > g.degree(v) / 2}

    # computes the LO polynomial corresponding to the union of the graphs with at least one of the identities added
    def union_lo(partition, finished, identities):
        out = zero
        for i in range(1, len(identities) + 1):
            for choices in combinations(identities, i):
                new_partition = partition.copy()
                for choice in choices:
                    new_partition = new_partition.join_many(choice)
                out += (1 if i % 2 == 1 else -1) * lo_helper(new_partition, finished)

        return out

    # tries to infer a coarser partition and a larger set of finished vertices
    # before calling the cached helper function
    def lo_helper(partition, finished):
        updated = True
        while updated:
            updated = False
            for v in g.nodes:
                if v not in finished:
                    # does v share a color with a strict majority of neighbors?
                    if sum(partition.find(v) == partition.find(w) for w in g.neighbors(v)) > (g.degree(v) / 2):
                        finished = finished.union({v})
                        updated = True
                    else:
                        # do a weak majority of the non-same-color neighbors of v share a color?
                        non_blue_neighbors = [w for w in g.neighbors(v) if partition.find(w) != partition.find(v)]
                        w_groups = group_by_key(non_blue_neighbors, lambda w: partition.find(w))
                        if len(w_groups) >= 1 and len(w_groups[0]) >= g.degree(v) / 2:
                            partition = partition.join_two(v, w_groups[0][0])
                            updated = True
                            if len(w_groups[0]) > g.degree(v) / 2:
                                finished = finished.union({v})
                            if len(w_groups) >= 2 and len(w_groups[1]) >= g.degree(v) / 2:
                                partition = partition.join_two(v, next(w for w in non_blue_neighbors
                                                                       if partition.find(w) == w_groups[1][0]))

        return lo_helper_cached(partition, frozenset(finished))

    # this recursive helper function does most of the computation
    # caches results to save time on repeated calls
    # requires inputs to be frozen
    @cache
    def lo_helper_cached(partition, finished):
        finished = set(finished)
        if len(finished) == len(g.nodes):
            return k(len(partition))
        else:
            v = next(v for v in strategy(g) if v not in finished)
            neighbors = list(g.neighbors(v))
            bs = {b for b in neighbors if partition.find(v) == partition.find(b)}
            ws = {w for w in neighbors if partition.find(v) != partition.find(w)}
            w_groups = group_by_key(ws, lambda w: partition.find(w))
            w_reps = [group[0] for group in w_groups]
            if len(bs) == 0:  # at least two neighbors should share a color with v
                return sum(((1 if j % 2 == 0 else -1) * (j - 1) *
                            sum((lo_helper(partition.join_many((v,) + new_same_neighbors), finished)
                                 for new_same_neighbors in combinations(ws, j)), start=zero)
                            for j in range(2, len(ws) + 1)), start=zero)
            elif len(bs) == 1:  # at least one more neighbor should share a color with v
                return sum((-(1 if j % 2 == 0 else -1) *
                            sum((lo_helper(partition.join_many((v,) + new_same_neighbors), finished)
                                 for new_same_neighbors in combinations(w_reps, j)), start=zero)
                            for j in range(1, len(w_reps) + 1)), start=zero)
            else:  # we can't assume any more neighbors share a color with v
                return lo_helper(partition, finished.union({v})) \
                    - union_lo(partition,
                               finished.union({v}),
                               set(combinations(ws, len(bs))).union(
                                   {(v, w) for w in w_reps})) \
                    + sum((-(1 if j % 2 == 0 else -1) *
                           sum((lo_helper(partition.join_many((v,) + new_same_neighbors), finished)
                                for new_same_neighbors in combinations(w_reps, j)), start=zero)
                           for j in range(1, len(w_reps) + 1)), start=zero)

    return lo_helper(*initial_values())
