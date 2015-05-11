# Имя файла Create_dict теперь не очень подходит. Потом исправлю

import sys
import os
import subprocess as sp
import collections

def main():
    input_filename, output_folder = parse_args()

    # create output folder
    os.makedirs(output_folder, exist_ok=True)

    rna_collection = create_mRNA_structure_collection(input_filename)
    print(rna_collection)

    rnafold_collection = run_rnafold_for_collection(rna_collection)
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

    RnaStructureInfo = collections.namedtuple("mRNA_tuple", ["name", "seq",
                                                       "structure"])

    with open(path_to_input, "r") as input_file:
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

                    m_t = RnaStructureInfo(name_of_seq, seq, structure)
                    mrna_collection.append(m_t)
            except StopIteration:
                break

    return mrna_collection

def run_rnafold_for_one_seq(name, seq):
    """Run RNAfold for one sequence
    Takes collection of RNA structures as input.
    Returns tuple containing RNAfold predicted structure and its mfe,
    """

    print("\nRunning RNAfold for {}...\n".format(name))

    rnafold_proc = sp.Popen("RNAfold", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)

    rnafold_raw_result = rnafold_proc.communicate(seq.encode())
    rnafold_result = rnafold_raw_result[0].decode().splitlines()

    print("\nRNAfold worked successfully!\n")

    draft_folding_string = rnafold_result[1].split(" ")[0]
    mfe = rnafold_result[1].split(" ")[1][1:][:-1]

    return draft_folding_string, mfe

def run_rnafold_for_collection(input_collection):
    """Run RNAfold for collection of sequences
    Returns list, containing named tuples of RNAfold predicted structures names,
    predicted structures in dot-bracket form and their mfe.
    """

    rnafold_collection = []

    RNAfold_result = collections.namedtuple("RNAfold_tuple", ["structure_name",
                                            "structure", "mfe"])

    for el in input_collection:
        structure_name = "{}_RNAfold_output".format(el.name)

        draft_folding_string, mfe = run_rnafold_for_one_seq(el.name, el.seq)

        r_t = RNAfold_result(structure_name, draft_folding_string, mfe)
        rnafold_collection.append(r_t)

    print(rnafold_collection)
    return rnafold_collection


if __name__ == '__main__':
    main()