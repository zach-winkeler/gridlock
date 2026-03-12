import unittest
from random import randrange
from unittest import TestCase

from networkx import Graph, grid_graph, random_unlabeled_tree, petersen_graph
from numpy import array
from scipy.special import stirling2

from recurrence import lo, Poly, k


# computes the falling factorial power of k
def ff(power):
    if power == 1:
        return k(1)
    else:
        return Poly([1 - power, 1]) * ff(power - 1)


def to_falling_basis(poly):
    n = poly.degree()
    mat = array([[stirling2(b, a) for b in range(n + 1)] for a in range(n + 1)])
    return mat @ poly.coef


def from_adjacency_list(adj_list):
    g = Graph()
    for v, adj in enumerate(adj_list):
        g.add_edges_from([(v, w) for w in adj])
    return g


class TestLO(TestCase):

    def test1(self):
        g = Graph([(1, 2), (1, 3), (2, 3), (4, 5), (4, 6), (5, 6), (7, 8), (7, 9), (8, 9),
                   (10, 11), (10, 12), (11, 12),
                   (1, 4), (1, 7), (1, 10), (4, 10), (7, 10)])
        self.assertEqual(ff(4) + 2 * ff(2) + ff(1), lo(g))

    def test2(self):
        g = Graph([(1, 2), (1, 3), (2, 3), (4, 5), (4, 6), (5, 6), (7, 8), (7, 9), (8, 9),
                   (10, 11), (10, 12), (11, 12),
                   (1, 4), (1, 7), (1, 10), (4, 7), (4, 10), (7, 10)])
        self.assertEqual(ff(4) + 3 * ff(2) + ff(1), lo(g))

    def test3(self):
        g = Graph([(1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7), (7, 8), (7, 9)])
        self.assertEqual(ff(3) + 2 * ff(2) + ff(1), lo(g))

    def test4(self):
        g = Graph([(1, 2), (1, 3), (2, 3), (4, 5), (4, 6), (5, 6), (7, 8), (7, 9), (8, 9),
                   (10, 11), (10, 12), (11, 12), (13, 14), (13, 15), (14, 15), (16, 17), (16, 18), (17, 18),
                   (1, 5), (4, 8), (7, 11), (10, 14), (13, 17), (16, 2),
                   (3, 12), (6, 15), (9, 18)])
        self.assertEqual(k(6), lo(g))

    def test5(self):
        g = Graph([(1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                   (2, 3), (3, 4), (4, 5), (5, 6), (6, 2),
                   (2, 5), (3, 7), (4, 7), (3, 8), (4, 8)])
        self.assertEqual(k(2), lo(g))

    def test_grid(self):
        g = grid_graph(dim=(4, 4))
        self.assertEqual(k(4) - 2 * k(3) + k(2) + k(1), lo(g))

    def test_path(self):
        g = Graph([(1, 2), (2, 3), (3, 4), (4, 5),
                   (1, 6), (1, 7), (1, 8),
                   (2, 9), (2, 10),
                   (3, 11), (3, 12),
                   (4, 13), (4, 14), (4, 15)])
        self.assertEqual(k(4) - 2 * k(3) + 3 * k(2) - k(1), lo(g))

    def test_tree(self):
        g = from_adjacency_list([[1], [0, 2, 4], [1, 5], [6], [1, 7], [2], [3, 7], [4, 6, 10], [11], [12], [7, 11, 13],
                                 [8, 10], [9, 13], [10, 12, 14, 16], [13], [16], [13, 15, 17], [16]])

        self.assertEqual(2 * k(3) - 2 * k(2) + k(1), lo(g))

    def test_petersen(self):
        g = petersen_graph()
        self.assertEqual(6 * ff(2) + ff(1), lo(g))

    def test_random_trees(self):
        for _ in range(10):
            n = randrange(3, 50)
            g = random_unlabeled_tree(n)
            p = lo(g)
            self.assertEqual(0, p.evaluate(0))
            self.assertEqual(1, p.evaluate(1))
            self.assertLess(p.degree(), (n // 3) + 1)

    def test_matt(self):
        g = from_adjacency_list(
            [[1, 2, 3, 6], [0, 2, 4, 7], [0, 1, 5, 8], [4, 5, 0, 9], [3, 5, 1, 10], [3, 4, 2, 11], [7, 8, 0, 9],
             [6, 8, 1, 10], [6, 7, 2, 11], [10, 11, 3, 6], [9, 11, 4, 7], [9, 10, 5, 8]])

        self.assertEqual(k(4) - k(3) - 2 * k(2) + 3 * k(1), lo(g))

    def test_matt2(self):
        g = from_adjacency_list(
            [[1, 2, 4, 7], [0, 2, 6, 7], [0, 1, 3, 4], [4, 5, 2, 11], [3, 5, 0, 2], [3, 4, 9, 11], [7, 8, 1, 10],
             [6, 8, 0, 1], [6, 7, 9, 10], [10, 11, 5, 8], [9, 11, 6, 8], [9, 10, 3, 5]])

        self.assertEqual(k(1), lo(g))

    # ladder-shaped graph, non-alternating LO-polynomial
    def test_nonalternating(self):
        g = Graph(
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
             (6, 7), (7, 8), (8, 9), (9, 10), (10, 11),
             (0,6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)])

        self.assertEqual(ff(3) + 4 * ff(2) + ff(1), lo(g))


if __name__ == '__main__':
    unittest.main()
