# Имя файла Create_dict теперь не очень подходит. Потом исправлю

import sys
import os
import subprocess as sp
import collections

def main():
    abs_in_path, abs_out_path = parse_args()

    # create output folder
    os.makedirs("{}".format(abs_out_path), exist_ok=True)

    rna_collection = create_mRNA_structure_collection(abs_in_path)
    print(rna_collection)

    rnafold_collection = run_proc_rnafold(rna_collection)
    print(rnafold_collection)

def parse_args():
    """Get paths to the input file and output directory
    and convert them to absolute paths
    """

    input_path = sys.argv[1]

    output_path="./mRNA_structure_project_out"
    if len(sys.argv) == 3:
        output_path = sys.argv[2]
    elif len(sys.argv) != 2:
        sys.exit("ERROR! Wrong number of arguments!")

    abs_in_path = os.path.abspath(input_path)
    abs_out_path = os.path.abspath(output_path)

    return abs_in_path, abs_out_path

def create_mRNA_structure_collection(path_to_input):
    """Takes input file and creates list, containing named tuples of
    experimentally proved RNA structures names, their sequences and
    structures in dot-bracket form.
    """

    mrna_collection = []

    with open(path_to_input, "r") as input_file:
        lines = input_file.readlines()
        lines_iter = iter(lines)

    i = 0

    for i in range(len(lines)):
        try:
            line = next(lines_iter)

            if not line.strip():
                continue
            if line[0] == ">":
                name_of_seq = line[1:].strip()
                line = next(lines_iter)

                first_sym = line[0]
                while first_sym < "A" or first_sym > "Z":
                    line = next(lines_iter)

                seq = line.strip()
                line = next(lines_iter)

                first_sym = line[0]
                while first_sym != "(" and first_sym != ".":
                    line = next(lines_iter)

                structure = line.strip()

                mRNA_tuple = collections.namedtuple("mRNA_tuple", ["name_of_seq",
                                                    "seq", "structure"])
                m_t = mRNA_tuple(name_of_seq, seq, structure)
                mrna_collection.append(m_t)
            i += 1
        except StopIteration:
            break

    return mrna_collection

def run_proc_rnafold(input_collection):
    """Run RNAfold
    Takes collection of RNA structures as input.
    Returns list, containing named tuples of RNAfold predicted structures names,
    predicted structures in dot-bracket form and their mfe.
    """

    rnafold_collection = []

    for el in input_collection:
        print("\nRunning RNAfold for {}...\n".format(el.name_of_seq))

        rnafold_proc = sp.Popen("RNAfold", shell=True, stdin=sp.PIPE,
                                stdout=sp.PIPE)

        rnafold_raw_result = rnafold_proc.communicate(el.seq.encode())
        rnafold_result = rnafold_raw_result[0].decode().splitlines()

        print("\nRNAfold worked successfully!\n")

        draft_folding_string = rnafold_result[1].split(" ")[0]
        mfe = rnafold_result[1].split(" ")[1][1:][:-1]
        structure_name = "{}_RNAfold_output".format(el.name_of_seq)

        RNAfold_tuple = collections.namedtuple("RNAfold_tuple",
                                        ["structure_name", "structure", "mfe"])
        r_t = RNAfold_tuple(structure_name, draft_folding_string, mfe)
        rnafold_collection.append(r_t)

    print(rnafold_collection)
    return rnafold_collection


if __name__ == '__main__':
    main()