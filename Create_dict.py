# Имя файла Create_dict теперь не очень подходит. Потом исправлю

import sys
import os
import subprocess as sp
from collections import namedtuple
import tempfile
import shutil

def main():
    input_filename, output_folder = parse_args()

    # create output folder
    os.makedirs(output_folder, exist_ok=True)

    rna_collection = create_mRNA_structure_collection(input_filename)
    print(rna_collection)

    rnafold_collection = run_rnafold_for_collection(rna_collection)
    print(rnafold_collection)

    mfold_collection = run_mfold_for_collection(rna_collection)
    print(mfold_collection)

def parse_args():
    """Get paths to the input file and output directory
    and convert them to absolute paths
    """

    input_filename = sys.argv[1]

    output_folder="./mRNA_structure_project_out"
    if len(sys.argv) == 3:
        output_folder = sys.argv[2]
    elif len(sys.argv) != 2:
        sys.exit("ERROR! Wrong number of arguments!")

    abs_in_path = os.path.abspath(input_filename)
    abs_out_path = os.path.abspath(output_folder)

    return abs_in_path, abs_out_path


def create_mRNA_structure_collection(input_filename):
    """Takes input file and creates list, containing named tuples of
    experimentally proved RNA structures names, their sequences and
    structures in dot-bracket form.
    """

    mrna_collection = []

    RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                                       "real_structure"])

    with open(input_filename, "r") as input_file:
        for line in input_file:
            try:
                if not line.strip():
                    continue
                if line[0] == ">":
                    name_of_seq = line[1:].strip()
                    line = next(input_file)

                    first_sym = line[0]
                    while first_sym.upper() < "A" or first_sym.upper() > "Z":
                        line = next(input_file)

                    seq = line.strip()
                    line = next(input_file)

                    first_sym = line[0]
                    while first_sym != "(" and first_sym != ".":
                        line = next(input_file)

                    structure = line.strip()

                    mRNA_tuple = RnaStructureInfo(name_of_seq, seq, structure)
                    mrna_collection.append(mRNA_tuple)
            except StopIteration:
                break

    return mrna_collection

def run_rnadistance(structure_name, pred_structure, real_structure):
    """Compare predicted structure with real structure using
     RNAdistance tool.
    """

    # create temporary folder for RNAdistance intermediary files
    dist_inter= tempfile.mkdtemp()

    # remember path to the main working directory
    main_folder = os.getcwd()

    # enter the temporary folder for intermediary files
    os.chdir(dist_inter)

    path_to_file = "./{}_vs_real_structure.txt".format(structure_name)

    with open(path_to_file, "w") as file_for_dist:
        file_for_dist.write(pred_structure)
        file_for_dist.write('\n')
        file_for_dist.write(real_structure)

    rnadist_proc = sp.Popen("RNAdistance < {}".format(path_to_file),
                                                 shell=True, stdout=sp.PIPE)
    rnadist_result = rnadist_proc.communicate()[0].decode().strip()

    # remove temporary folder
    shutil.rmtree(dist_inter)

    # return to the main working directory
    os.chdir(main_folder)

    return rnadist_result

def run_rnafold_for_one_seq(name, seq, real_structure):
    """Run RNAfold and RNAdistance for one sequence
    Takes RNA nucleotide sequence with its name and real structure as input.
    Returns tuple containing RNAfold predicted structure, its mfe and distance
    to real structure,
    """

    print("\nRunning RNAfold for {}...\n".format(name))

    rnafold_proc = sp.Popen("RNAfold", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)

    rnafold_raw_result = rnafold_proc.communicate(seq.encode())
    rnafold_result = rnafold_raw_result[0].decode().splitlines()[1].split()

    print("\nRNAfold worked successfully!\n")

    draft_folding_string = rnafold_result[0]
    mfe = rnafold_result[1][1:][:-1]

    dist = run_rnadistance(name, draft_folding_string, real_structure)

    return draft_folding_string, mfe, dist

def run_rnafold_for_collection(input_collection):
    """Run RNAfold for collection of sequences
    Returns dictionary, containing name of the structure as a key and
    named tuple RNAfold predicted structure in dot-bracket form, its
    mfe and distance to real structure as a value.
    """

    rnafold_collection = dict()

    RNAfold_result = namedtuple("RNAfold_result", ["structure", "mfe",
                                                   "distance"])

    for el in input_collection:
        structure_name = el.name

        draft_folding_string, mfe, dist = run_rnafold_for_one_seq(el.name,
                                                el.seq, el.real_structure)

        rnafold_collection[structure_name] = RNAfold_result(draft_folding_string,
                                                            mfe, dist)
    return rnafold_collection

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

    mfold_result = namedtuple("mfold_result", ["structure_name",
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

    mfold_1_proc = sp.Popen("mfold SEQ='./{}.fasta'".format(name), shell=True)
    mfold_1_proc.wait()

    # script Ct2B.pl should be located in the main working directory

    # mfold intermediary files contain the name of input file,

    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                                                     name), shell=True)
    mfold_2_proc.wait()

    print("\nmfold worked successfully!\n")

    mfold_list = []

    with open("{}/{}.b".format(mfold_inter_out, name)) as mfold_output:
       i = 0
       for line in mfold_output:
           if i != 0:
               # dot-bracket structure and the value of mfe are
               # divided by tabulation
               raw_result = line.strip().split("\t")
               structure_name = str(i)
               draft_folding_string = raw_result[0]
               mfe = raw_result[1][1:][:-1]
               dist = run_rnadistance(name, draft_folding_string, real_structure)
               mfold_tuple = mfold_result(structure_name, draft_folding_string, mfe, dist)
               mfold_list.append(mfold_tuple)
           i += 1

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

    mfold_dict = dict()

    for el in input_collection:
        mfold_dict[el.name] = run_mfold_for_one_seq(el.name, el.seq, el.real_structure)

    return mfold_dict



if __name__ == '__main__':
    main()