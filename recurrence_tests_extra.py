import unittest
from random import randrange
from unittest import TestCase

from networkx import Graph, grid_graph, random_unlabeled_tree, petersen_graph
from numpy import array
from scipy.special import stirling2

from recurrence import lo, Poly, k
from recurrence_tests import from_adjacency_list


class TestLO(TestCase):

    def test_matt3(self):
        g = from_adjacency_list(
            [[1, 2, 3, 4, 10], [0, 2, 3, 5, 11], [0, 1, 3, 6, 8], [0, 1, 2, 7, 9], [5, 6, 7, 0, 14], [4, 6, 7, 1, 15],
             [4, 5, 7, 2, 12], [4, 5, 6, 3, 13], [9, 10, 11, 2, 12], [8, 10, 11, 3, 13], [8, 9, 11, 0, 14],
             [8, 9, 10, 1, 15], [13, 14, 15, 6, 8], [12, 14, 15, 7, 9], [12, 13, 15, 4, 10], [12, 13, 14, 5, 11]])

        self.assertEqual(2 * k(4) - 6 * k(3) + 14 * k(2) - 9 * k(1), lo(g))

    def test_matt4(self):
        g = from_adjacency_list(
            [[1, 2, 3, 4, 5], [0, 2, 3, 8, 9], [0, 1, 3, 4, 5], [0, 1, 2, 8, 9], [5, 6, 7, 0, 2], [4, 6, 7, 0, 2],
             [4, 5, 7, 12, 14], [4, 5, 6, 12, 14], [9, 10, 11, 1, 3], [8, 10, 11, 1, 3], [8, 9, 11, 13, 15],
             [8, 9, 10, 13, 15], [13, 14, 15, 6, 7], [12, 14, 15, 10, 11], [12, 13, 15, 6, 7], [12, 13, 14, 10, 11]])

        self.assertEqual(2 * k(4) + 8 * k(3) - 16 * k(2) + 7 * k(1), lo(g))


if __name__ == '__main__':
    unittest.main()
