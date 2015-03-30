import sys
import subprocess as sp
import os

out_path='./mRNA_structure_project_out'

if len(sys.argv) == 3:
    out_path = sys.argv[2]
elif len(sys.argv) != 2:
    sys.exit("Wrong number of arguments!")

if out_path[0] == ".":
    full_out_path = os.getcwd() + out_path[1:]
else:
    full_out_path = out_path

in_path = sys.argv[1]

if in_path[0] == ".":
    full_in_path = os.getcwd() + in_path[1:]
else:
    full_in_path = in_path

with open (in_path, 'r') as input_file:
    _ = 0
    for line in input_file:
        if _ == 0:
            seq_name = line[1:]
        _ += 1

main_folder = os.getcwd()

file_name = (sys.argv[1].split("/")[len(sys.argv[1].split("/")) - 1]).split(".")[0]

os.makedirs(out_path, exist_ok=True)

os.mkdir("{}/mfold_inter_out".format(out_path))

os.chdir("{}/mfold_inter_out".format(out_path))

mfold_1_proc = sp.Popen("mfold SEQ='{}'".format(full_in_path),
                        shell=True, stdout=sp.PIPE)

mfold_1_proc.wait()

mfold_2_proc = sp.Popen("{}/Ct2B.pl {}.ct > {}.b".format(main_folder,
                        file_name, file_name), shell=True, stdout=sp.PIPE)

mfold_2_proc.wait()

with open("./{}.b".format(file_name)) as mfold_result:
    _ = 0
    for line in mfold_result:
        if _ != 0:
            with open("{}/{}_mfold_output_{}.fasta".format(full_out_path,
                      seq_name, _), 'w') as mfold_out:
                mfold_out.write(line.split("\t")[0])
        _ += 1