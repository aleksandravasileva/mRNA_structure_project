#!/usr/bin/env python3

import sys
import subprocess as sp
import os

def main():
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

    running_rnafold(parse_paths(input_path, output_path)[0],
                    parse_paths(input_path, output_path)[1])

    running_mfold(parse_paths(input_path, output_path)[0],
                  parse_paths(input_path, output_path)[1])
    processing_mfold_output(parse_paths(input_path, output_path)[0],
                            parse_paths(input_path, output_path)[1])

def parse_paths(in_path, out_path):
    full_in_path = os.path.abspath(in_path)
    full_out_path = os.path.abspath(out_path)
    return full_in_path, full_out_path

def parse_names(in_path, out_path):
    file_name_with_format = os.path.basename(parse_paths(in_path, out_path)[0])
    file_name = os.path.splitext(file_name_with_format)[0]
    with open (in_path, "r") as input_file:
        seq_name = input_file.readline()[1:]
        # reading name of RNA sequence
    return file_name, seq_name

def running_rnafold(in_path, out_path):
    os.makedirs(out_path, exist_ok=True)
    rnafold_proc = sp.Popen("RNAfold < {}".format(in_path),
                               shell=True, stdout=sp.PIPE)
    rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()
    print("RNAfold worked successfully!")
    draft_folding_string = rnafold_result[2].split(" ")[0]
    mfe = rnafold_result[2].split(" ")[1]
    with open("{}/{}_RNAfold_output.fasta".format(out_path,
                         parse_names(in_path, out_path)[1]), "w") as rnafold_out:
        rnafold_out.write(draft_folding_string)
        rnafold_out.write("\n")
        rnafold_out.write(mfe)

def running_mfold(in_path, out_path):
    os.makedirs("{}/mfold_inter_out".format(out_path), exist_ok=True)
        # creating folder for mfold intermediary results
    main_folder = os.getcwd()
        # remembering path to the main working directory
    os.chdir("{}/mfold_inter_out".format(out_path))
        # entering the folder for intermediary results
        # mfold writes its results only in working directory
    mfold_1_proc = sp.Popen("mfold SEQ='{}'".format(in_path), shell=True)
    mfold_1_proc.wait()
    mfold_2_proc = sp.Popen("{0}/Ct2B.pl {1}.ct > {1}.b".format(main_folder,
                            parse_names(in_path, out_path)[0]), shell=True)
        # script Ct2B.pl should be located in the main working directory
        # mfold intermediary files contain the name of input file,
        # that's why file_name variable is used
    mfold_2_proc.wait()
    print("mfold worked successfully!")
    os.chdir(main_folder)
        # returning to the main working directory

def processing_mfold_output(in_path, out_path):
    with open("{}/mfold_inter_out/{}.b".format(out_path,
        parse_names(in_path, out_path)[0])) as mfold_result:
        i = 0
        for line in mfold_result:
            if i != 0:
                with open("{}/{}_mfold_output_{}.fasta".format(out_path,
                         parse_names(in_path, out_path)[1], i), "w") as mfold_out:
                      # main mfold output files will contain the name
                      # of RNA sequence
                    draft_folding_string = line.split("\t")[0]
                    mfe = line.split("\t")[1]
                      # dot-bracket structure and the value of mfe are
                      # divided by tabulation
                    mfold_out.write(draft_folding_string)
                    mfold_out.write("\n")
                    mfold_out.write(mfe)
            i += 1

if __name__ == '__main__':
    main()
