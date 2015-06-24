import os
import subprocess as sp
import tempfile
import shutil
import rnadistance as rnadist
import logging
import datatypes as dt


def get_structs_for_one_seq(name, mRNA_tuple, name_of_logger):
    """Run mfold and RNAdistance for one sequence
    Takes named tuple of RNA nucleotide sequence with its name and real
    structure as input.
    mfold returns several structures for each sequence (optimal and
    suboptimal).
    As a result this function returns list containing named tuples of
    counting number of mfold predicted structure, predicted structure
    in dot-bracket form, its nucleotide sequence, mfe and distances to
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

    with open('./{}.fasta'.format(name), 'w') as mfold_input:
        mfold_input.write('>')
        mfold_input.write(name)
        mfold_input.write('\n')
        mfold_input.write(mRNA_tuple.seq)
        mfold_input.write('\n')

    global module_logger
    module_logger = logging.getLogger(name_of_logger + '.mfold')

    module_logger.info("Running mfold for %s...", name)

    mfold_1_proc = sp.Popen("mfold SEQ='./{}.fasta'".format(name),
                            shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    mfold_stdout, mfold_stderr = mfold_1_proc.communicate()
    mfold_error = mfold_stderr.decode()
    if mfold_error != '':
        module_logger.warning("mfold worked with errors:\n%s...", mfold_error)

    #Script Ct2B.pl should be located in the main working directory

    #mfold intermediary files contain the name of input file

    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                            name), shell=True)
    mfold_2_proc.wait()

    module_logger.info("mfold finished its work")

    mfold_list = []

    with open("{}/{}.b".format(mfold_inter_out, name)) as mfold_output:
        #First line of mfold output is the RNA sructure nucleotide sequence
        mfold_output_lines = mfold_output.readlines()[1:]
        for index, string in enumerate(mfold_output_lines):
            #dot-bracket structure and the value of mfe are
            #divided by tabulation
            raw_result = string.strip().split("\t")

            #mfold returns several structures: the first structure (index=0)
            #is optimal, all others are suboptimal
            if index == 0:
                is_optimal = True
            else:
                is_optimal = False
            seq = mRNA_tuple.seq
            folding_string = raw_result[0]
            mfe = float(raw_result[1][1:][:-1])
            tree_dist = rnadist.get_tree_dist(folding_string,
                                              mRNA_tuple.structure)
            n_tree_dist = rnadist.get_n_tree_dist(folding_string,
                                                  mRNA_tuple.structure)
            bp_dist = rnadist.get_bp_dist(folding_string,
                                              mRNA_tuple.structure)
            n_bp_dist = rnadist.get_n_bp_dist(folding_string,
                                                  mRNA_tuple.structure)

            mfold_tuple = dt.FoldResult(is_optimal, seq, folding_string,
                                        mfe, tree_dist, n_tree_dist,
                                        bp_dist, n_bp_dist)
            mfold_list.append(mfold_tuple)

    #Remove temporary folder
    shutil.rmtree(mfold_inter_out)
    #Return to the main working directory
    os.chdir(main_folder)

    return mfold_list


def get_structs_for_collection(input_collection, name_of_logger):
    """Run mfold for collection of sequences.
    Returns dictionary, containing name of the structure as a key and list of
    named tuples of mfold predicted structure counting number, predicted
    structure in dot-bracket form, its nucleotide sequence, mfe and distances to
    real structure as a value.
    """
    return {key: get_structs_for_one_seq(key, value, name_of_logger)
            for key, value in input_collection.items()}


def find_best_structure(mfold_list):
    """mfold predicts several optimal and suboptimal foldings.
    This function determines the best predicted structure (the closest
    to the real structure according RNAdistance value).
    Returns named tuple containing the best mfold predicted structure
    in dot-bracket form, its nucleotide sequence, mfe and distances to
    real structure.
    """

    min_tree_dist = mfold_list[0].tree_dist
    result_optimality = mfold_list[0].is_optimal
    result_seq = mfold_list[0].seq
    result_structure = mfold_list[0].structure
    result_mfe = mfold_list[0].mfe
    result_n_tree_dist = mfold_list[0].n_tree_dist
    result_bp_dist = mfold_list[0].bp_dist
    result_n_bp_dist = mfold_list[0].n_bp_dist

    for el in mfold_list:
        if el.tree_dist < min_tree_dist:
            min_tree_dist = el.tree_dist
            result_optimality = el.is_optimal
            result_structure = el.structure
            result_mfe = el.mfe
            result_n_tree_dist = el.n_tree_dist

    mfold_best_result_tuple = dt.FoldResult(result_optimality, result_seq,
                                            result_structure, result_mfe,
                                            min_tree_dist, result_n_tree_dist,
                                            result_bp_dist, result_n_bp_dist)

    return mfold_best_result_tuple


def get_best_structs_collection(input_collection, name_of_logger):
    """Returns dictionary of the best mfold predicted structures for
    the collection of sequences.
    Name of the structure is a key and named tuple of the best mfold
    predicted structure in dot-bracket form, its nucleotide sequence,
    mfe and distances to real structure is a value.
    """

    mfold_collection = get_structs_for_collection(input_collection,
                                                  name_of_logger)

    module_logger.info('Created mfold predicted structures collection')
    module_logger.debug('RNAfold predicted structures collection:\n%s',
                 mfold_collection)

    return {key: find_best_structure(value)
            for key, value in mfold_collection.items()}