import subprocess as sp
from rnadistance import get_distance
import logging
import datatypes as dt


def run_for_one_seq(mRNA_tuple):
    """Run RNAfold and RNAdistance for one sequence
    Takes RNA nucleotide sequence with its name and real structure as input.
    Returns named tuple containing RNAfold predicted structure in dot-bracket
    form, its nucleotide sequence, mfe and distance to real structure,
    """

    module_logger = logging.getLogger('mRNA_structure_project.RNAfold')

    module_logger.info("Running RNAfold for %s...", mRNA_tuple.name)

    rnafold_proc = sp.Popen("RNAfold", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE, stderr=sp.PIPE)

    rnafold_stdout, rnafold_stderr = rnafold_proc.communicate(mRNA_tuple.seq.encode())
    rnafold_result = rnafold_stdout.decode().splitlines()[1].split()
    #Rnafold returns input nucleotide sequence (1st line) and its structure
    #with mfe (2nd line).

    rnafold_error = rnafold_stderr.decode()
    if rnafold_error != '':
        module_logger.warning("RNAfold worked with errors:\n%s...", rnafold_error)

    module_logger.info("RNAfold finished its work")

    #RNAfold returns only one optimal structure
    optimality = "optimal"

    folding_string = rnafold_result[0]
    #RNAfold returns mfe as '(mfe_value)'
    mfe = float(rnafold_result[1][1:][:-1])

    dist = get_distance(folding_string, mRNA_tuple.real_structure)

    rnafold_tuple = dt.FoldResult(optimality, mRNA_tuple.seq, folding_string,
                                  mfe, dist)

    return rnafold_tuple


def run_for_collection(input_collection):
    """Run RNAfold for collection of sequences
    Returns dictionary, containing name of the structure as a key and
    named tuple of RNAfold predicted structure in dot-bracket form, its
    nucleotide sequence, mfe and distance to real structure as a value.
    """

    return {el.name: run_for_one_seq(el) for el in input_collection}