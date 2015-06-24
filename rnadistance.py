import subprocess as sp


def get_dist(first_structure, second_structure, dist_option):
    """Compare two RNA structures using RNAdistance tool
    """
    dist_input = first_structure + '\n' + second_structure

    rnadist_proc = sp.Popen("RNAdistance --distance={}".format(dist_option),
                            shell=True, stdin=sp.PIPE, stdout=sp.PIPE)

    rnadist_stdout = rnadist_proc.communicate(dist_input.encode())
    #RNAdistance returns distance as 'dist_option: distance_value \n'
    rnadist_result = int(rnadist_stdout[0].decode().strip().split()[1])

    return rnadist_result


def get_tree_dist(first_structure, second_structure):
    """Compare two RNA structures using RNAdistance tool
    (tree editing distance).
    """
    return get_dist(first_structure, second_structure, "f")


def get_n_tree_dist(first_structure, second_structure):
    """Compare two RNA structures using RNAdistance tool
    (tree editing distance).
    Returns normalized value.
    """
    # distance value normalizes by comparing it to the sum of the
    # distances between structures and roots of the trees

    dist_two_struc = get_tree_dist(first_structure, second_structure)
    dist_1_struc = get_tree_dist(first_structure, ".")
    dist_2_struc = get_tree_dist(".", second_structure)

    return dist_two_struc/(dist_1_struc+dist_2_struc)


def get_bp_dist(first_structure, second_structure):
    """Compare two RNA structures using RNAdistance tool
    (base pair distance).
    """
    return get_dist(first_structure, second_structure, "P")


def get_n_bp_dist(first_structure, second_structure):
    """Compare two RNA structures using RNAdistance tool
    (base pair distance).
    Returns normalized value
    """
    # distance value normalizes by the length of the sequence

    length = len(first_structure)

    return get_bp_dist(first_structure, second_structure)/length