import sys
import rnafold
import mfold
import logging
import os
import datatypes
from rnaeval import get_mfe

def main():
    input_filename, output_folder = parse_args()

    #create main output folder
    os.makedirs(output_folder, exist_ok=True)

    # logger output will be written in the file in output_folder
    # named as output_folder
    name_of_logger = os.path.basename(output_folder)
    logger = create_logger(name_of_logger, output_folder)

    rna_collection = create_mRNA_structure_collection(input_filename)
    logger.info('Created RNA structures collection')
    logger.debug('RNA structures collection:\n%s', rna_collection)

    rnafold_collection = rnafold.get_structs_for_collection(rna_collection,
                                                            name_of_logger)
    logger.info('Created RNAfold predicted structures collection')
    logger.debug('RNAfold predicted structures collection:\n%s',
                 rnafold_collection)

    best_mfold_collection = mfold.get_best_structs_collection(rna_collection,
                                                              name_of_logger)
    logger.info('Created collection of the best mfold predicted structures')
    logger.debug('Best mfold predicted structures collection:\n%s',
                 best_mfold_collection)

    #create folder for forna files
    forna_folder = "{}/forna".format(output_folder)
    os.makedirs(forna_folder, exist_ok=True)

    create_forna_files_for_collection(rna_collection, forna_folder, "Real")

    logger.info('Created files of real structures for forna')

    create_forna_files_for_collection(rnafold_collection, forna_folder,
                                      "RNAfold")

    logger.info('Created files of RNAfold predictions for forna')

    create_forna_files_for_collection(best_mfold_collection, forna_folder,
                                      "mfold")

    logger.info('Created files of mfold predictions for forna')


def create_logger(name, output_folder):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    logging_file = output_folder + '/' + name + '.log'

    fh = logging.FileHandler(logging_file, mode='w')
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

    return logger


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

    mfe = get_mfe(seq, structure)
    length = len(seq)

    mRNA_tuple = datatypes.RnaStructureInfo(seq, structure, length, mfe)

    return mRNA_tuple, number, name_of_seq


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

    Returns dictionary, containing name of the structure as a key and
    named tuple of experimentally proved structure in dot-bracket form,
    its nucleotide sequence, length of this sequence and energy of the
    structure
    """

    mrna_collection = {}

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
                    mrna_collection[one_struct_info[2]] = one_struct_info[0]
                    i = one_struct_info[1]
                i += 1
            except StopIteration:
                break

    return mrna_collection


def create_forna_file(output_folder, origin, name, seq, structure):
    """Creates file that can be used as an input of RNA structures visualization
    tool - forna (http://nibiru.tbi.univie.ac.at/forna/).
    File contains name of the structure, nucleotide sequence and the structure
    itself in dot-bracket form.
    The name of the used predicting tool should be passed as an argument.
    """
    if origin == "Real":
        forna_file = '{}/{}_(Real).txt'.format(output_folder, name)
    else:
        forna_file = '{}/{}_({}_predicted).txt'.format(output_folder, name, origin)
    with open(forna_file, 'w') as output:
        if origin == "Real":
            output.write('>{}_Real'.format(name))
        else:
            output.write('>{}_{}_predicted'.format(name, origin))
        output.write('\n')
        output.write(seq)
        output.write('\n')
        output.write(structure)


def create_forna_files_for_collection(input_collection, output_folder, origin):
    """Creates forna files for the collection of predicted structures.
    """
    for key, value in input_collection.items():
        create_forna_file(output_folder, origin, key, value.seq,
                          value.structure)


if __name__ == '__main__':
    main()