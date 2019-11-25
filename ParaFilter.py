from ete3 import Tree
import ParaModules
import argparse


usage = "A script to filter paralogs."
toolname = "ParaFilter"
footer = "Who \n Mattia Giacomelli (mattia.giacomelli@bristol.ac.uk); \n Paschalis Natsidis (p.natsidis@ucl.ac.uk); \n \nWhere \n Pisani Lab, Uni Bristol; \n Telford Lab, UCL;\n\
 ITN IGNITE; \n  \nWhen\n October 2019; \n\n"

parser = argparse.ArgumentParser(description = usage, prog = toolname, epilog = footer, formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.add_argument('-f', metavar = 'filename', dest = 'fasta_file', required = True,
                    help = 'full path to fastas file')
parser.add_argument('-t', metavar = 'filename', dest = 'trees_file', required = True,
                    help = 'full path to trees file')
parser.add_argument('-n', metavar = 'int', dest = 'characters', required = True,
                    help = 'determine how many characters of a fasta header consist a species ID')
parser.add_argument('-w', metavar = 'directory', dest = 'working_dir', required = True,
                    help = 'full path to working directory')

#parser.print_help()

args = parser.parse_args()

#READ USER INPUT
fasta_input = args.fasta_file
trees_input = args.trees_file
characters_number = int(args.characters)
wdir = args.working_dir

fastas_txt_file = open(fasta_input, "r")
trees_txt_file = open(trees_input, "r")

fastas_lines = fastas_txt_file.readlines()
trees_lines = trees_txt_file.readlines()

fastas = [x.strip() for x in fastas_lines]
trees = [x.strip() for x in trees_lines]


for i in range(len(fastas)):
    
    current_tree = Tree(wdir + trees[i])
    current_fasta = wdir + fastas[i]

    paralogs_tuples = ParaModules.GetParalogsTuples(current_tree, characters_number)    # GET TUPLES WITH PARALOGS OF THE TREE

    monophyletic = ParaModules.IsMonophyletic(current_tree, paralogs_tuples)            # GET NATURE OF EACH SET OF PARALOGS (MONOPHYLETIC ETC.

    to_remove = ParaModules.GetParalogsToRemove(current_tree, monophyletic)             # DETERMINE WHICH PARALOGS WILL BE REMOVED

    ParaModules.RemoveSequences(current_tree, wdir+current_fasta, to_remove)            # WRITE NEW FASTA WITH PARALOGS REMOVED