#!/usr/bin/env python3

import sys
import subprocess as sp
import os
import os.path
import tempfile
import shutil

def main():
    abs_in_path, abs_out_path = parse_args()

    # create output folder
    os.makedirs("{}".format(abs_out_path), exist_ok=True)

    result_folding_rnafold = run_proc_rnafold(abs_in_path)

    result_folding_mfold = run_proc_mfold(abs_in_path)

    run_rnadistance(abs_in_path, abs_out_path, result_folding_rnafold)

    run_rnadistance(abs_in_path, abs_out_path, result_folding_mfold)

def parse_args():
    """Get paths to the input file and output directory
    and convert them to absolute paths
    """

    input_path = sys.argv[1]
    fasta_extensions = ['.fasta', '.fas', '.fsa', '.fa']
    if not any(input_path.lower().endswith(extension)
                  for extension in fasta_extensions):
        sys.exit("ERROR! Wrong input file format!")

    output_path="./mRNA_structure_project_out"
    if len(sys.argv) == 3:
        output_path = sys.argv[2]
    elif len(sys.argv) != 2:
        sys.exit("ERROR! Wrong number of arguments!")

    abs_in_path = os.path.abspath(input_path)
    abs_out_path = os.path.abspath(output_path)

    return abs_in_path, abs_out_path

def get_file_name(in_path):
    """Return name of the file"""

    file_name_with_format = os.path.basename(in_path)
    file_name = os.path.splitext(file_name_with_format)[0]

    return file_name

def get_seq_name(in_path):
    """Return name of the sequence"""

    with open (in_path, "r") as input_file:
        seq_name = input_file.readline()[1:].strip()

    return seq_name

def get_real_structure(in_path):
    """Get real sequence structure in the dot-bracket form
     from the input file
    """

    with open(in_path, "r") as input_fasta:
        for line in input_fasta:
            if line[0] == "." or line[0] == "(":
                return line.strip()

def run_proc_rnafold(in_path):
    """Run RNAfold
    Returns list of lists, containing predicted
    structure in dot-bracket form and mfe.
    """

    rnafold_proc = sp.Popen("RNAfold < {}".format(in_path), shell=True,
                            stdout=sp.PIPE)
    rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()
    print("\nRNAfold worked successfully!\n")

    draft_folding_string = rnafold_result[2].split(" ")[0]
    mfe = rnafold_result[2].split(" ")[1]

    seq_name = get_seq_name(in_path)

    result_folding_list = ["{}_RNAfold_output".format(seq_name)]

    result_folding_list.append([draft_folding_string, mfe])

    return result_folding_list

def run_proc_mfold(in_path):
    """Run mfold
    Returns list of lists, containing predicted
    structure in dot-bracket form and mfe.
    """

    # create temporary folder for mfold intermediary results
    mfold_inter_out = tempfile.mkdtemp()

    # remember path to the main working directory
    main_folder = os.getcwd()

    # mfold writes its results only in working directory!
    # enter the temporary folder for intermediary results
    os.chdir(mfold_inter_out)

    mfold_1_proc = sp.Popen("mfold SEQ='{}'".format(in_path), shell=True)
    mfold_1_proc.wait()

    # script Ct2B.pl should be located in the main working directory

    # mfold intermediary files contain the name of input file,
    # that's why file_name variable is used
    file_name = get_file_name(in_path)
    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                                                     file_name), shell=True)
    mfold_2_proc.wait()
    print("\nmfold worked successfully!\n")

    seq_name = get_seq_name(in_path)

    result_folding_list = ["{}_mfold_output".format(seq_name)]

    with open("{}/{}.b".format(mfold_inter_out, file_name)) as mfold_result:
        i = 0
        for line in mfold_result:
            if i != 0:
                # dot-bracket structure and the value of mfe are
                # divided by tabulation
                draft_folding_string = line.split("\t")[0]
                mfe = line.split("\t")[1].strip()
                result_folding_list.append([draft_folding_string, mfe])
            i += 1

    return result_folding_list

    # remove temporary folder
    shutil.rmtree(mfold_inter_out)

    # return to the main working directory
    os.chdir(main_folder)

def run_rnadistance(in_path, out_path, result_folding_tool):
    """Compare all predicted structures with real structure using
     RNAdistance tool.
    """

    print("Comparing {} with real structure...".format(result_folding_tool[0]))

    # create folder for RNAdistance main result
    os.makedirs("{}/RNAdistance_out".format(out_path), exist_ok=True)

    # create temporary folder for RNAdistance intermediary files
    dist_inter_out = tempfile.mkdtemp()

    # for each predicted structure create file containing
    # this structure and the real structure
    # this file is an input for the RNAdistance
    real_structure = get_real_structure(in_path)
    i = 0
    for structure in result_folding_tool:
        if i != 0:
            with open("{}/{}_{}".format(dist_inter_out, result_folding_tool[0],
                                                     i), "w") as file_for_dist:
                file_for_dist.write(structure[0])
                file_for_dist.write("\n")
                file_for_dist.write(real_structure)
        i += 1

    with open ("{}/RNAdistance_out/{}_distances.txt".format(out_path,
                      result_folding_tool[0]), "w") as dist_out_file:
        for file in os.listdir(dist_inter_out):
            path_to_file = "{}/{}".format(dist_inter_out, file)
            rnadist_proc = sp.Popen("RNAdistance < {}".format(path_to_file),
                                                 shell=True, stdout=sp.PIPE)
            rnadist_result = rnadist_proc.communicate()[0].decode()

            dist_out_file.write(file)
            dist_out_file.write("\n")
            dist_out_file.write(rnadist_result.strip())
            dist_out_file.write("\n")
            dist_out_file.write("\n")

    print("\nRNAdistance worked successfully!\n")

    # remove temporary folder
    shutil.rmtree(dist_inter_out)

if __name__ == '__main__':
    main()