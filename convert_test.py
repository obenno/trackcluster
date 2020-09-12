#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/2018 6:52 PM
# @Author  : Runsheng     
# @File    : convert_test.py

import unittest
from convert import *


class GffConvertTest(unittest.TestCase):
    def setUp(self):
        self.gff = GFF("./test/genes/unc52/unc52.gff")

    def test_convert(self):
        bigg_list= gff_to_bigGenePred(self.gff)

        with open("./test/genes/unc52/unc52_gff.bed", "w") as fw:
            for bigg in bigg_list:
                fw.write(bigg.to_str())
                fw.write("\n")

    def test_cigar_count(self):
        cigar_tuple = [(4, 95), (0, 28), (2, 1), (0, 14), (2, 2), (0, 5), (2, 1),
                       (0, 17), (2, 2), (0, 16), (1, 1), (0, 4), (1, 1), (0, 10),
                       (2, 1), (0, 18), (2, 1), (0, 4), (2, 1), (0, 43), (2, 2),
                       (0, 14), (1, 1), (0, 13), (2, 1), (0, 6), (1, 1), (0, 21),
                       (1, 1), (0, 1), (1, 1), (0, 3), (2, 1), (0, 6), (2, 2),
                       (0, 1), (2, 1), (0, 4), (3, 836), (0, 20), (1, 1),
                       (0, 14), (1, 2), (0, 1), (1, 1), (0, 12), (2, 1), (0, 1),
                       (2, 1), (0, 14), (2, 1), (0, 19), (1, 1), (0, 18), (3, 740),
                       (0, 12), (2, 2), (0, 31), (2, 1), (0, 14), (1, 1), (0, 33),
                       (2, 1), (0, 2), (2, 1), (0, 18), (1, 3), (0, 33), (2, 2),
                       (0, 11), (1, 3), (0, 7), (2, 1), (0, 5), (2, 3), (0, 28),
                       (1, 1), (0, 4), (2, 2), (0, 11), (1, 1), (0, 2), (2, 1),
                       (0, 32), (1, 1), (0, 2), (2, 2), (0, 1), (2, 2), (0, 10),
                       (1, 1), (0, 12), (2, 1), (0, 4), (3, 3399), (0, 1), (2, 2),
                       (0, 35), (2, 2), (0, 7), (1, 2), (0, 17), (1, 2), (0, 18),
                       (1, 2), (0, 1), (1, 2), (0, 14), (1, 1), (0, 15), (2, 3),
                       (0, 5), (3, 248), (0, 17), (1, 3), (0, 12), (1, 1), (0, 4),
                       (2, 1), (0, 19), (2, 1), (0, 9), (2, 2), (0, 13), (1, 5), (0, 8),
                       (1, 1), (0, 10), (2, 2), (0, 7), (1, 3), (0, 2), (2, 1), (0, 3),
                       (2, 1), (0, 6), (1, 1), (0, 6), (1, 1), (0, 7), (4, 15)]

        print(cigar_count(cigar_tuple))


    def tearDown(self):
        self.bigg = None