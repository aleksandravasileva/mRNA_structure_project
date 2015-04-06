import sys
import os
import subprocess as sp

def main():
    abs_in_path, abs_out_path = parse_args()

    # create output folder
    os.makedirs("{}".format(abs_out_path), exist_ok=True)

    mrna_dict = create_dict(abs_in_path)

    rnafold_dict = run_proc_rnafold(mrna_dict)

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

def create_dict(in_put):
    """Takes input file and creates dictionary, containing name of the sequence
    as the key and the list with RNA sequence and its real structure in
    dot-bracket form as the value.
    """

    mrna_dict = dict()

    with open(in_put, "r") as input_file:
        line = input_file.readlines()
        i = 0

        while i < len(line):
            if line[i].strip() == "":
                i += 1
                continue

            if line[i][0] == ">":
                name_of_seq = line[i][1:].strip()
                i += 1

                first_sym = line[i][0]
                while ord(first_sym) < 65 or ord(first_sym) > 90:
                    i += 1
                seq = line[i].strip()
                i += 1

                first_sym = line[i][0]
                while first_sym != "(" and first_sym != ".":
                    i += 1
                structure = line[i].strip()

                mrna_dict[name_of_seq] = [seq, structure]
            i += 1
    return mrna_dict

def run_proc_rnafold(input_dict):
    """Run RNAfold
    Returns dictionary, containing (name of the sequence + "_RNAfold_output")
    as the key and list with structure in dot-bracket form and mfe
    as the value.
    """

    rna_fold_dict = dict()
    for key, value in input_dict.items():
        print("\nRunning RNAfold for {}...\n".format(key))

        rnafold_proc = sp.Popen(["RNAfold"], shell=True, stdin=value[0],
                                stdout=sp.PIPE)

        # И вот с этого момента ничего не хочет работать.
        # Я не понимаю, как перехватить правильно stdin
        
        # AttributeError: 'str' object has no attribute 'fileno'
        # Что бы вот это все значило...


        rnafold_result = rnafold_proc.communicate()[0].decode().splitlines()

        print("\nRNAfold worked successfully!\n")
        draft_folding_string = rnafold_result[2].split(" ")[0]
        mfe = rnafold_result[2].split(" ")[1]
        rnafold_key_name = "{}_RNAfold_output".format(key)
        rna_fold_dict[rnafold_key_name] = [draft_folding_string, mfe]

    print(rna_fold_dict)
    return rna_fold_dict


if __name__ == '__main__':
    main()