from collections import namedtuple

RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                  "real_structure"])

RNAfoldResult = namedtuple("RNAfoldResult", ["sequence", "structure",
                                "mfe", "distance"])

MfoldResult = namedtuple("MfoldResult", ["structure_name", "sequence",
                              "structure", "mfe", "distance"])

MfoldBestResult = namedtuple("MfoldBestResult", ["sequence",
                                   "structure", "mfe", "distance"])


