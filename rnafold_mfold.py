#!/usr/bin/env python3

import sys
import subprocess as sp
import os

def main():
    abs_in_path, abs_out_path = parse_paths(*parse_args())

    file_name, seq_name = parse_names(abs_in_path, abs_out_path)

    processing_rnafold_output(abs_in_path, abs_out_path, seq_name,
                              *running_rnafold(abs_in_path, abs_out_path))

    running_mfold(abs_in_path, abs_out_path, file_name)
    processing_mfold_output(abs_in_path, abs_out_path, seq_name)

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
    abs_in_path = os.path.abspath(in_path)
    abs_out_path = os.path.abspath(out_path)
    return abs_in_path, abs_out_path

def parse_names(in_path, out_path):
    # read name of the input file
    file_name_with_format = os.path.basename(parse_paths(in_path, out_path)[0])
    file_name = os.path.splitext(file_name_with_format)[0]
    # read name of RNA sequence
    with open (in_path, "r") as input_file:
        seq_name = input_file.readline()[1:]
    return file_name, seq_name

def running_rnafold(in_path, out_path):
    os.makedirs(out_path, exist_ok=True)
    rnafold_proc = sp.Popen("RNAfold < {}".format(in_path),
                               shell=True, stdout=sp.PIPE)
    rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()
    print("RNAfold worked successfully!")
    draft_folding_string = rnafold_result[2].split(" ")[0]
    mfe = rnafold_result[2].split(" ")[1]
    return draft_folding_string, mfe

def processing_rnafold_output(in_path, out_path, seq_name,
                              draft_folding_string, mfe):
    with open("{}/{}_RNAfold_output.fasta".format(out_path,
                            seq_name), "w") as rnafold_out:
        rnafold_out.write(draft_folding_string)
        rnafold_out.write("\n")
        rnafold_out.write(mfe)

def running_mfold(in_path, out_path, file_name):
    # create folder for mfold intermediary results
    os.makedirs("{}/mfold_inter_out".format(out_path), exist_ok=True)
    # remember path to the main working directory
    main_folder = os.getcwd()
    # mfold writes its results only in working directory,
    # that's why enter the folder for intermediary results
    os.chdir("{}/mfold_inter_out".format(out_path))
    mfold_1_proc = sp.Popen("mfold SEQ='{}'".format(in_path), shell=True,
                                          stdout=sp.PIPE, stderr=sp.PIPE)
    mfold_1_proc.wait()
    # script Ct2B.pl should be located in the main working directory
    # mfold intermediary files contain the name of input file,
    # that's why file_name variable is used
    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                     file_name), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    mfold_2_proc.wait()
    print("mfold worked successfully!")
    # return to the main working directory
    os.chdir(main_folder)

def processing_mfold_output(in_path, out_path, seq_name):
    with open("{}/mfold_inter_out/{}.b".format(out_path,
        parse_names(in_path, out_path)[0])) as mfold_result:
        i = 0
        for line in mfold_result:
            if i != 0:
                # main mfold output files should contain the name
                # of RNA sequence,
                # that's why seq_name variable is used
                with open("{}/{}_mfold_output_{}.fasta".format(out_path,
                         seq_name, i), "w") as mfold_out:
                    # dot-bracket structure and the value of mfe are
                    # divided by tabulation
                    draft_folding_string = line.split("\t")[0]
                    mfe = line.split("\t")[1]
                    mfold_out.write(draft_folding_string)
                    mfold_out.write("\n")
                    mfold_out.write(mfe)
            i += 1

if __name__ == '__main__':
    main()
