from collections import namedtuple

RnaStructureInfo = namedtuple("RnaStructureInfo", ["seq", "structure",
                                                   "length", "mfe"])

FoldResult = namedtuple("FoldResult", ["is_optimal", "seq", "structure", "mfe",
                                       "tree_dist", "n_tree_dist", "bp_dist",
                                       "n_bp_dist"])