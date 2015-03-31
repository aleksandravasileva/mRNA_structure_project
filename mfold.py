import sys
import subprocess as sp
import os

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
                            file_name), shell=True)
        # script Ct2B.pl should be located in the main working directory
        # mfold intermediary files contain the name of input file,
        # that's why file_name variable is used
    mfold_2_proc.wait()
    os.chdir(main_folder)
        # returning to the main working directory

def processing_mfold_output(out_path):
    with open("{}/mfold_inter_out/{}.b".format(out_path, file_name))\
            as mfold_result:
        i = 0
        for line in mfold_result:
            if i != 0:
                with open("{}/{}_mfold_output_{}.fasta".format(out_path,
                                             seq_name, i), "w") as mfold_out:
                      # main mfold output files will contain the name
                      # of RNA sequence
                    mfold_out.write(line.split("\t")[0])
                      # dot-bracket structure and the value of mfe are
                      # divided by tabulation
                    mfold_out.write("\n")
                    mfold_out.write(line.split("\t")[1])
            i += 1

output_path="./mRNA_structure_project_out"

if len(sys.argv) == 3:
    output_path = sys.argv[2]
elif len(sys.argv) != 2:
    sys.exit("ERROR! Wrong number of arguments!")

input_path = sys.argv[1]

full_in_path = os.path.abspath(input_path)
full_out_path = os.path.abspath(output_path)

file_name = os.path.splitext(os.path.basename(full_in_path))[0]

with open (full_in_path, "r") as input_file:
    seq_name = input_file.readline()[1:]
    # reading name of RNA sequence

running_mfold(full_in_path, full_out_path)
processing_mfold_output(full_out_path)

