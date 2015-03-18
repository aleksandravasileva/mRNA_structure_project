import subprocess as sp

rnafold_proc = sp.Popen("RNAfold < /home/aleksandra/Документы/Bioinformatics_Institute/Spring/Python/mRNA_structure_project/Week_2/HIV-1_IRES.fasta", shell=True, stdin=sp.PIPE, stdout=sp.PIPE)
# Я чего-то так и не разобралась, как сделать так, чтобы fasta-файл считывался с клавиатуры и передавался в Popen

rnafold_proc.wait()
rnafold_result_0 = rnafold_proc.communicate()


rnafold_result = str(rnafold_result_0[0])

print(rnafold_result)
#b'>HIV-1_IRES\nAUGGGUGCGAGAGCGUCGGUAUUAAGCGGGGGAGAAUUAGAUAAAUGGGAAAAAAUUCGGUUAAGGCCAGGGGGAAAGAAACAAUAUAAACUAAAACAUAUAGUAUGGGCAAGCAGGGAGCUAGAACGAUUCGCAGUUAAUCCUGGCCUUUUAGAGACAUCAGAAGGCUGUAGACAAAUACUGGGACAGCUACAACCAUCCCUUCAGACAGGAUCAGAAGAACUUAGAUCAUUAUAUAAUACAAUAGCAGUCCUCUAUUGUGUGCAUCAAAGGAUAGAUGUAAAAGACACCAAGGAAGCCUUAGAUAAGAUAGAGGAAGGACAAAACAAAAGUAAGAAAAAGGCACAGCAAGCAGCAGCUGACACAGGAAACAACAGCCAGGUCAGCCAAAAUUACCCUAUAGUGCAGAACCUCCAGGGGCAAAUGGUACAUCAGGCCAUAUCACCUAGAACU\n.(((((((((.....)))((((((...((((((..........(((.(((.....))).)))((((((((((.....(((.(.......((((.......))))..(..(.(((.....))).)..)).)))........))))))))))...((....))....((((((...((.....))..)))))).......)))))).......((((...........))))......)))))).....(..((((((((((((((((((........))))).....)))).((((...))))......)))))))))..)...........................((.....)).((((((.(.((.........)).)))))))........((((..((........))..))))....(((((.......)))))..))))))..... (-88.70)\n'
#Нужна только часть с dot-brackets

rnafold_result = rnafold_result.split("'")
# Если не убрать кавычки в строке, ничего не работает

rnafold_result = rnafold_result[1].split('\\n')
print(rnafold_result)
#['>HIV-1_IRES', 'AUGGGUGCGAGAGCGUCGGUAUUAAGCGGGGGAGAAUUAGAUAAAUGGGAAAAAAUUCGGUUAAGGCCAGGGGGAAAGAAACAAUAUAAACUAAAACAUAUAGUAUGGGCAAGCAGGGAGCUAGAACGAUUCGCAGUUAAUCCUGGCCUUUUAGAGACAUCAGAAGGCUGUAGACAAAUACUGGGACAGCUACAACCAUCCCUUCAGACAGGAUCAGAAGAACUUAGAUCAUUAUAUAAUACAAUAGCAGUCCUCUAUUGUGUGCAUCAAAGGAUAGAUGUAAAAGACACCAAGGAAGCCUUAGAUAAGAUAGAGGAAGGACAAAACAAAAGUAAGAAAAAGGCACAGCAAGCAGCAGCUGACACAGGAAACAACAGCCAGGUCAGCCAAAAUUACCCUAUAGUGCAGAACCUCCAGGGGCAAAUGGUACAUCAGGCCAUAUCACCUAGAACU', '.(((((((((.....)))((((((...((((((..........(((.(((.....))).)))((((((((((.....(((.(.......((((.......))))..(..(.(((.....))).)..)).)))........))))))))))...((....))....((((((...((.....))..)))))).......)))))).......((((...........))))......)))))).....(..((((((((((((((((((........))))).....)))).((((...))))......)))))))))..)...........................((.....)).((((((.(.((.........)).)))))))........((((..((........))..))))....(((((.......)))))..))))))..... (-88.70)', '']

rnafold_result = rnafold_result[2].split(' ')
print(rnafold_result)
#['.(((((((((.....)))((((((...((((((..........(((.(((.....))).)))((((((((((.....(((.(.......((((.......))))..(..(.(((.....))).)..)).)))........))))))))))...((....))....((((((...((.....))..)))))).......)))))).......((((...........))))......)))))).....(..((((((((((((((((((........))))).....)))).((((...))))......)))))))))..)...........................((.....)).((((((.(.((.........)).)))))))........((((..((........))..))))....(((((.......)))))..)))))).....', '(-88.70)']
#Нужно еще убрать значения mfe в конце

rnafold_output = rnafold_result[0]

print(rnafold_output)

with open('/home/aleksandra/Документы/Bioinformatics_Institute/Spring/Python/mRNA_structure_project/Week_2/RNAfold_output.fasta', 'w') as RNAfold_file:
    RNAfold_file.write(rnafold_output)
