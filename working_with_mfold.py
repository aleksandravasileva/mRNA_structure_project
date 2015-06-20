import os
import subprocess as sp
import tempfile
import shutil
from working_with_RNAdistance import run_rnadistance
import logging
import datatypes as dt


def run_mfold_for_one_seq(mRNA_tuple):
    """Run mfold and RNAdistance for one sequence
    Takes named tuple of RNA nucleotide sequence with its name and real
    structure as input.
    mfold returns several structures for each sequence (optimal and
    suboptimal).
    As a result this function returns list containing named tuples of
    counting number of mfold predicted structure, predicted structure
    in dot-bracket form, its nucleotide sequence, mfe and distance to
    real structure.
    """

    #Temporary file in fasta format containing RNA sequence should be created
    #This file will be passed to mfold as an input

    #Create temporary folder for mfold input and intermediary results
    mfold_inter_out = tempfile.mkdtemp()

    #Remember path to the main working directory
    main_folder = os.getcwd()

    #mfold writes its results only in working directory!
    #Enter the temporary folder for intermediary results
    os.chdir(mfold_inter_out)

    with open('./{}.fasta'.format(mRNA_tuple.name), 'w') as mfold_input:
        mfold_input.write('>')
        mfold_input.write(mRNA_tuple.name)
        mfold_input.write('\n')
        mfold_input.write(mRNA_tuple.seq)
        mfold_input.write('\n')

    module_logger = logging.getLogger('mRNA_structure_project.mfold')

    module_logger.info("Running mfold for %s...", mRNA_tuple.name)

    mfold_1_proc = sp.Popen("mfold SEQ='./{}.fasta'".format(mRNA_tuple.name),
                            shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    mfold_stdout, mfold_stderr = mfold_1_proc.communicate()
    mfold_error = mfold_stderr.decode()
    if mfold_error != '':
        module_logger.warning("mfold worked with errors:\n%s...", mfold_error)

    #Script Ct2B.pl should be located in the main working directory

    #mfold intermediary files contain the name of input file

    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                            mRNA_tuple.name), shell=True)
    mfold_2_proc.wait()

    module_logger.info("mfold finished its work")

    mfold_list = []

    with open("{}/{}.b".format(mfold_inter_out, mRNA_tuple.name)) as mfold_output:
        #First line of mfold output is the RNA sructure nucleotide sequence
        mfold_output_lines = mfold_output.readlines()[1:]
        for index, string in enumerate(mfold_output_lines):
            #dot-bracket structure and the value of mfe are
            #divided by tabulation
            raw_result = string.strip().split("\t")
            structure_name = str(index + 1)
            sequence = mRNA_tuple.seq
            folding_string = raw_result[0]
            mfe = float(raw_result[1][1:][:-1])
            dist = run_rnadistance(folding_string, mRNA_tuple.real_structure)
            mfold_tuple = dt.MfoldResult(structure_name, sequence, folding_string, mfe,
                                       dist)
            mfold_list.append(mfold_tuple)

    #Remove temporary folder
    shutil.rmtree(mfold_inter_out)
    #Return to the main working directory
    os.chdir(main_folder)

    return mfold_list


def run_mfold_for_collection(input_collection):
    """Run mfold for collection of sequences.
    Returns dictionary, containing name of the structure as a key and list of
    named tuples of mfold predicted structure counting number, predicted
    structure in dot-bracket form, its nucleotide sequence, mfe and distance to
    real structure as a value.
    """
    return {el.name: run_mfold_for_one_seq(el) for el in input_collection}


def find_best_mfold_structure(mfold_list):
    """mfold predicts several optimal and suboptimal foldings.
    This function determines the best predicted structure (the closest
    to the real structure according RNAdistance value).
    Returns named tuple containing the best mfold predicted structure
    in dot-bracket form, its nucleotide sequence, mfe and distance to
    real structure.
    """

    min_dist = mfold_list[0].distance
    result_sequence = mfold_list[0].sequence
    result_structure = mfold_list[0].structure
    result_mfe = mfold_list[0].mfe
    for el in mfold_list:
        if el.distance < min_dist:
            min_dist = el.distance
            result_structure = el.structure
            result_mfe = el.mfe

    mfold_best_result_tuple = dt.MfoldBestResult(result_sequence,
                                                result_structure, result_mfe, min_dist)

    return mfold_best_result_tuple


def create_mfold_best_structures_collection(input_collection):
    """Returns dictionary of the best mfold predicted structures for
    the collection of sequences.
    Name of the structure is a key and named tuple of the best mfold
    predicted structure in dot-bracket form, its nucleotide sequence,
    mfe and distance to real structure is a value.
    """
    return {key: find_best_mfold_structure(value) for
            key, value in input_collection.items()}