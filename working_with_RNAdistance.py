import subprocess as sp


def run_rnadistance(pred_structure, real_structure):
    """Compare predicted structure with real structure using
     RNAdistance tool.
    """
    dist_input = pred_structure + '\n' + real_structure

    rnadist_proc = sp.Popen("RNAdistance", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)

    rnadist_stdout = rnadist_proc.communicate(dist_input.encode())
    #RNAdistance returns distance as 'f: distance_value'
    rnadist_result = rnadist_stdout[0].decode().strip().split()[1]

    return rnadist_result
