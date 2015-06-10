import subprocess as sp
from collections import namedtuple
from working_with_RNAdistance import run_rnadistance
import logging


def run_rnafold_for_one_seq(name, seq, real_structure):
    """Run RNAfold and RNAdistance for one sequence
    Takes RNA nucleotide sequence with its name and real structure as input.
    Returns named tuple containing RNAfold predicted structure in dot-bracket
    form, its mfe and distance to real structure as a value,
    """

    module_logger = logging.getLogger('mRNA_structure_project.RNAfold')

    module_logger.info("Running RNAfold for %s...", name)

    rnafold_proc = sp.Popen("RNAfold", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)

    rnafold_stdout = rnafold_proc.communicate(seq.encode())
    rnafold_result = rnafold_stdout[0].decode().splitlines()[1].split()
    #Rnafold returns input nucleotide sequence (1st line) and its structure
    #with mfe (2nd line).

    module_logger.info("RNAfold finished its work")

    RNAfold_result = namedtuple("RNAfold_result", ["structure", "mfe",
                                "distance"])

    folding_string = rnafold_result[0]
    #RNAfold returns mfe as '(mfe_value)'
    mfe = rnafold_result[1][1:][:-1]

    dist = run_rnadistance(folding_string, real_structure)

    rnafold_tuple = RNAfold_result(folding_string, mfe, dist)

    return rnafold_tuple


def run_rnafold_for_collection(input_collection):
    """Run RNAfold for collection of sequences
    Returns dictionary, containing name of the structure as a key and
    named tuple RNAfold predicted structure in dot-bracket form, its
    mfe and distance to real structure as a value.
    """

    rnafold_collection = {}

    for el in input_collection:
        structure_name = el.name

        rnafold_collection[structure_name] = run_rnafold_for_one_seq(el.name,
                                                                     el.seq, el.real_structure)
    return rnafold_collection
