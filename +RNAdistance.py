#!/usr/bin/env python3

import sys
import subprocess as sp
import os
import os.path

def main():
    abs_in_path, abs_out_path = parse_paths(*parse_args())

    run_proc_rnafold(abs_in_path, abs_out_path)

    running_mfold(abs_in_path, abs_out_path)
    processing_mfold_output(abs_in_path, abs_out_path)

    run_rnadistance(abs_in_path, abs_out_path)

def parse_args():
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

    return input_path, output_path

def parse_paths(in_path, out_path):
    """Return absolute paths to the input file and output directory"""

    abs_in_path = os.path.abspath(in_path)
    abs_out_path = os.path.abspath(out_path)

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

def get_real_str(in_path, out_path):
    """Get real sequence structure in the dot-bracket form
     from the input file
    """
    with open(in_path, "r") as input_fasta:
        for line in input_fasta:
            if line[0] == "." or line[0] == "(":
                return line.strip()

def run_proc_rnafold(in_path, out_path):
    """Run RNAfold and get dot-bracket structure with mfe from its output"""

    # create folder for RNAfold main results (folder "Predictors_out"
    # inside output folder)
    os.makedirs("{}/Predictors_out".format(out_path), exist_ok=True)

    rnafold_proc = sp.Popen("RNAfold < {}".format(in_path), shell=True,
                            stdout=sp.PIPE)
    rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()
    print("\nRNAfold worked successfully!\n")

    draft_folding_string = rnafold_result[2].split(" ")[0]
    mfe = rnafold_result[2].split(" ")[1]

    # remember path to the main working directory
    main_folder = os.getcwd()

    # enter the main RNAfold results folder
    os.chdir("{}/Predictors_out".format(out_path))

    seq_name = get_seq_name(in_path)

    with open("./{}_RNAfold_output.fasta".format(seq_name),
              "w") as rnafold_out:
        rnafold_out.write(draft_folding_string)
        rnafold_out.write("\n")
        rnafold_out.write(mfe)

    # return to the main working directory
    os.chdir(main_folder)

def running_mfold(in_path, out_path):
    """Run mfold"""

    # create folders for mfold main (folder "Predictors_out")
    # and intermediary results (folder "mfold_inter_out")
    os.makedirs("{}/Predictors_out/mfold_inter_out".format(out_path),
                exist_ok=True)

    # remember path to the main working directory
    main_folder = os.getcwd()

    # mfold writes its results only in working directory!
    # enter the folder for intermediary results
    os.chdir("{}/Predictors_out/mfold_inter_out".format(out_path))

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

    # return to the main working directory
    os.chdir(main_folder)

def processing_mfold_output(in_path, out_path):
    """Get dot-bracket structure with mfe from mfold output"""

    # remember path to the main working directory
    main_folder = os.getcwd()

    # enter the main mfold results folder
    os.chdir("{}/Predictors_out".format(out_path))

    file_name = get_file_name(in_path)
    seq_name = get_seq_name(in_path)
    with open("./mfold_inter_out/{}.b".format(file_name)) as mfold_result:
        i = 0
        for line in mfold_result:
            if i != 0:
                # main mfold output files should contain the name
                # of RNA sequence,
                # that's why seq_name variable is used
                with open("./{}_mfold_output_{}.fasta".format(seq_name, i),
                          "w") as mfold_out:
                    # dot-bracket structure and the value of mfe are
                    # divided by tabulation
                    draft_folding_string = line.split("\t")[0]
                    mfe = line.split("\t")[1]

                    mfold_out.write(draft_folding_string)
                    mfold_out.write("\n")
                    mfold_out.write(mfe)
            i += 1

    # return to the main working directory
    os.chdir(main_folder)

def run_rnadistance(in_path, out_path):
    """Compare all predicted structures (from folder "Predictors_out")
    with real structure using RNAdistance tool.
    """

    # create folders for RNAdistance main results (folder "RNAdistance_out")
    # and intermediary files (folder "RNAdistance_inter")
    os.makedirs("{}/RNAdistance_out/RNAdistance_inter".format(out_path),
                exist_ok=True)

    # create list containing names of all files with
    # predicted structures (all files from "Predictors_out" folder)
    all_files = []
    pred_path = "{}/Predictors_out".format(out_path)
    for file in os.listdir(pred_path):
        if os.path.isfile(os.path.join(pred_path, file)):
            all_files.append(file)

    # remember path to the main working directory
    main_folder = os.getcwd()

    # enter the main RNAdistance results folder
    os.chdir("{}/RNAdistance_out".format(out_path))

    # for each predicted structure create new file
    # containing this structure and the real structure

    # this file is an input for the RNAdistance
    real_str = get_real_str(in_path, out_path)
    for file in all_files:
        with open("{}/Predictors_out/{}".format(out_path, file), "r")\
                as pred_str_file:
            file_name = get_file_name("{}/Predictors_out/{}".format(out_path,
                                                                    file))

            with open("./RNAdistance_inter/{}_for_dist.fasta".format(file_name),
                "w") as file_for_dist:
                file_for_dist.write(pred_str_file.readline().strip())
                file_for_dist.write("\n")
                file_for_dist.write(real_str)

    with open ("./All_distances.txt", "w") as dist_out_file:
        for file in os.listdir("./RNAdistance_inter"):
            path_to_file = "./RNAdistance_inter/{}".format(file)
            rnadist_proc = sp.Popen("RNAdistance < {}".format(path_to_file),
                                    shell=True, stdout=sp.PIPE)
            rnadist_result = rnadist_proc.communicate()[0].decode()

            dist_out_file.write(file)
            dist_out_file.write("\n")
            dist_out_file.write(rnadist_result.strip())
            dist_out_file.write("\n")
            dist_out_file.write("\n")

    print("\nRNAdistance worked successfully!\n")

    # return to the main working directory
    os.chdir(main_folder)

if __name__ == '__main__':
    main()