from collections import namedtuple

RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                  "real_structure"])

FoldResult = namedtuple("FoldResult", ["optimality", "sequence", "structure",
                        "mfe", "distance"])