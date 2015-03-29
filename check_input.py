import sys
import os

if len(sys.argv) == 3:
    out_path = sys.argv[2]
elif len(sys.argv) == 2:
    out_path = "./mRNA_structure_project_out"
else:
    sys.exit("Wrong number of arguments!")

if ".fasta" in sys.argv[1] or ".FASTA" in sys.argv[1] or \
   ".fas"   in sys.argv[1] or ".fsa" in sys.argv[1] or \
                              ".fa"  in sys.argv[1]:
    in_path = sys.argv[1]
else:
    sys.exit("Wrong input file format!")

try:
    os.mkdir(out_path)
except FileExistsError:
    pass




