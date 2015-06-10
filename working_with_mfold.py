import os
from collections import namedtuple
import subprocess as sp
import tempfile
import shutil
from working_with_RNAdistance import run_rnadistance
import logging


def run_mfold_for_one_seq(name, seq, real_structure):
    """Run mfold and RNAdistance for one sequence
    Takes RNA nucleotide sequence with its name and real structure as input.
    mfold returns several structures for each sequences (optimal and
    suboptimal).
    As a result this function returns list containing named tuples of
    counting number of mfold predicted structure,
    predicted structure in dot-bracket form, its mfe and distance to
    real structure.
    """

    mfold_result = namedtuple("Mfold_result", ["structure_name",
                              "structure", "mfe", "distance"])

    #Temporary file in fasta format containing RNA sequence should be created
    #This file will be passed to mfold as an input

    # create temporary folder for mfold input and intermediary results
    mfold_inter_out = tempfile.mkdtemp()

    # remember path to the main working directory
    main_folder = os.getcwd()

    # mfold writes its results only in working directory!
    # enter the temporary folder for intermediary results
    os.chdir(mfold_inter_out)

    with open('./{}.fasta'.format(name), 'w') as mfold_input:
        mfold_input.write('>')
        mfold_input.write(name)
        mfold_input.write('\n')
        mfold_input.write(seq)
        mfold_input.write('\n')

    module_logger = logging.getLogger('mRNA_structure_project.mfold')

    module_logger.info("Running mfold for %s...", name)

    mfold_1_proc = sp.Popen("mfold SEQ='./{}.fasta'".format(name), shell=True)
    mfold_1_proc.wait()

    # script Ct2B.pl should be located in the main working directory

    # mfold intermediary files contain the name of input file,

    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                            name), shell=True)
    mfold_2_proc.wait()

    module_logger.info("mfold finished its work")

    mfold_list = []

    with open("{}/{}.b".format(mfold_inter_out, name)) as mfold_output:
        #First line of mfold output is the RNA sructure nucleotide sequence
        mfold_output_lines = mfold_output.readlines()[1:]
        for index, string in enumerate(mfold_output_lines):
            # dot-bracket structure and the value of mfe are
            # divided by tabulation
            raw_result = string.strip().split("\t")
            structure_name = str(index + 1)
            folding_string = raw_result[0]
            mfe = raw_result[1][1:][:-1]
            dist = run_rnadistance(folding_string, real_structure)
            mfold_tuple = mfold_result(structure_name, folding_string, mfe,
                                       dist)
            mfold_list.append(mfold_tuple)

    # remove temporary folder
    shutil.rmtree(mfold_inter_out)
    # return to the main working directory
    os.chdir(main_folder)

    return mfold_list


def run_mfold_for_collection(input_collection):
    """Run mfold for collection of sequences
    Returns dictionary, containing name of the structure as a key and list of
    named tuples of mfold predicted structure counting number, predicted
    structure in dot-bracket form, its mfe and distance to real structure as
    a value.
    """

    mfold_dict = {}

    for el in input_collection:
        mfold_dict[el.name] = run_mfold_for_one_seq(el.name, el.seq,
                                                    el.real_structure)

    return mfold_dict
