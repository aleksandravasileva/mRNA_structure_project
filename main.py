import sys
from working_with_RNAfold import *
from working_with_mfold import *
import logging


def main():
    input_filename, output_folder = parse_args()

    # create output folder
    os.makedirs(output_folder, exist_ok=True)

    logger = logging.getLogger('mRNA_structure_project')
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler('mRNA_structure_project.log', mode='w')
    fh.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter_fh = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s'
                                     ' - %(message)s\n')
    formatter_ch = logging.Formatter('%(message)s\n')
    fh.setFormatter(formatter_fh)
    ch.setFormatter(formatter_ch)

    logger.addHandler(fh)
    logger.addHandler(ch)

    rna_collection = create_mRNA_structure_collection(input_filename)
    logger.debug('RNA structures collection:\n%s', rna_collection)
    logger.info('Created RNA structures collection')

    rnafold_collection = run_rnafold_for_collection(rna_collection)
    logger.debug('RNAfold predicted structures collection:\n%s',
                 rnafold_collection)

    mfold_collection = run_mfold_for_collection(rna_collection)
    logger.debug('RNAfold predicted structures collection:\n%s',
                 mfold_collection)


def parse_args():
    """Get paths to the input file and output directory
    and convert them to absolute paths
    """

    output_folder = "mRNA_structure_project_out"

    if len(sys.argv) == 3:
        output_folder = sys.argv[2]
    elif len(sys.argv) != 2:
        sys.exit("ERROR! Wrong number of arguments!\nUsage: ./main.py"
                 "input_filename (required) output_folder (optional)")

    input_filename = sys.argv[1]

    abs_input_filename = os.path.abspath(input_filename)
    abs_output_folder = os.path.abspath(output_folder)

    return abs_input_filename, abs_output_folder


def get_structure_from_input(input_lines, number):
    """Takes list of several lines from input file and counting number
    of the first line from this list.
    Selects first RNA structure from these lines.
    Returns named tuple of this RNA structure, its sequence and structure
    in dot-bracket form and number of line, at which function stopped.
    """
    RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                  "real_structure"])

    lines = iter(input_lines)

    for line in lines:
        name_of_seq = line[1:].strip()
        line = next(lines)
        number += 1
        first_sym = line[0]
        while first_sym.upper() < "A" or first_sym.upper() > "Z":
            line = next(lines)
            number += 1

        seq = line.strip()
        line = next(lines)
        number += 1

        first_sym = line[0]
        while first_sym != "(" and first_sym != ".":
            line = next(lines)
            number += 1

        structure = line.strip()
        break

    mRNA_tuple = RnaStructureInfo(name_of_seq, seq, structure)

    return mRNA_tuple, number


def create_mRNA_structure_collection(input_filename):
    """Takes input file, containing several mRNA_sructures with their names,
    nucleotide sequences and experimentally determined secondary structures
    in dot-bracket form.
    Example:
    >HIV-1_IRES
    AUGGGUGCGAGAGCGUCGGUAUU...
    ...(((((....))))).........
    >HIV-2_IRES
    AUGGGCGCGAGAAACUCCGUCUU...
    .....((.(((...)))))((((...

    Creates list, containing named tuples of experimentally proved RNA
    structures names, their sequences and structures in dot-bracket form.
    """

    RnaStructureInfo = namedtuple("RnaStructureInfo", ["name", "seq",
                                  "real_structure"])

    mrna_collection = []

    i = 0

    with open(input_filename, "r") as input_file:
        input_lines = input_file.readlines()
        while i < len(input_lines):
            try:
                line = input_lines[i]
                if not line.strip():
                    i += 1
                    continue
                if line[0] == ">":
                    one_struct_info = get_structure_from_input(input_lines[i:],
                                                               i)
                    mrna_collection.append(one_struct_info[0])
                    i = one_struct_info[1]
                i += 1
            except StopIteration:
                break

    return mrna_collection


if __name__ == '__main__':
    main()
