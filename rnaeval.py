import subprocess as sp

def get_mfe(sequence, structure):
    """Runs RNAeval for RNA sequence and its structure.
    Returns the mfe value.
    """
    rnaeval_input = sequence + '\n' + structure

    rnaeval_proc = sp.Popen("RNAeval", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)

    rnaeval_stdout = rnaeval_proc.communicate(rnaeval_input.encode())
    # RNAeval returns the following string:
    # 'sequence\nstructure (mfe_value)\n'
    mfe = float(rnaeval_stdout[0].decode().strip().split()[2][1:][:-1])

    return mfe