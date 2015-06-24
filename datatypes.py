from collections import namedtuple

RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                  "real_structure", "mfe"])

FoldResult = namedtuple("FoldResult", ["is_optimal", "sequence", "structure",
                        "mfe", "distance"])