import sys
import subprocess as sp
import os

path = sys.argv[1]
out_path = sys.argv[2]

#out_path путь к папке, которая должна быть создана
# и в которую будут записываться все промежуточные файлы

rnafold_proc = sp.Popen("RNAfold < {}".format(path), shell=True, stdout=sp.PIPE)

rnafold_proc.wait()
rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()

seq_name = rnafold_result[0].split(">")[1]

draft_folding_string = rnafold_result[2].split(" ")[0]

try:
    os.mkdir("{}".format(out_path))
except FileExistsError:
    print("Don't create directory next time")
finally:
    with open('{}/{}_RNAfold_output.fasta'.format(out_path, seq_name),
              'w') as rnafold_new_file:
        rnafold_new_file.write(draft_folding_string)

