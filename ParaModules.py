
import collections
import ete3

def GetParalogsTuples(tree, number_of_characters):
    '''
    This function takes as input a gene tree and returns a list of tuples.
    Each tuple in this list contains a set of in-paralogs.
    Each set of in-paralogs will be tested for their monophyly in later stages.
    '''
    leaf_names = [leaf for leaf in tree.get_leaf_names()]
    leaf_names_sorted = sorted(leaf_names)
    ids_list = [tip[0:number_of_characters] for tip in leaf_names_sorted]
    paralogs_with_counts = dict(collections.Counter(ids_list))  # BUILD A DICTIONARY WITH GENE COUNTS
    list_of_paralogs = []
    for x in leaf_names_sorted:
        for key, values in paralogs_with_counts.items():
            if values > 1:
                if x.startswith(key):
                    list_of_paralogs.append(x)                  # IF AN ID HAS MORE THAN ONE GENE, KEEP IT
    final_list = []
    for prefix in set(ids_list):
        with_prefix = []
        for gene in list_of_paralogs:
            if prefix in gene:
                with_prefix.append(gene)
        if len(with_prefix) >= 1:
            final_list.append(with_prefix)
    final_sorted = sorted(final_list)
    paralogs_tuples = [tuple(i) for i in final_sorted]
    return paralogs_tuples                

def IsMonophyletic(tree, tuples_list):
    '''
    This function accepts as input a tree and a list of tuples.
    Each tuple contains a set of paralogs. Each set of paralogs will be tested for their monophyly.
    The output is a dictionary with the tuples as the keys, and their nature as the value
    '''
    nature_of_paralogs = []
    for entry in tuples_list:    
        is_monophyletic = tree.check_monophyly(entry, "name", unrooted=True)    # CHECK THE MONOPHYLY OF THE TUPLE CONTENTS
        nature_of_paralogs.append(list(entry) + [is_monophyletic[1]])           # CREATE A LIST THAT HOLDS THE INFO IN LAST ELEMENT
    monophyly_dict = {}         # THIS WILL HOLD THE TUPLES AND THEIR NATURE
    for entry in nature_of_paralogs:
        paralogs = entry[0:-1]     # ALL ELEMENTS EXCEPT THE LAST ARE THE PARALOG NAMES 
        monophyly_dict[tuple(paralogs)] = entry[-1] # THE LAST ELEMENT IS THE NATURE OF THE PARALOGS. PUT IT AS VALUE IN THE DICT
    return monophyly_dict

def GetParalogsToRemove(tree, monophyly_dict):
    '''
    This function takes as input the dictionary with the nature of different tuples of paralogs.
    For monophyletic tuples, it keeps all paralogs except the one with the shortest branch,
    For non-monophyletic tuples, it keeps everything.
    It returns a list with all paralogs that will be removed.
    '''
    # SEPARATE MONOPHYLETIC AND NON-MONOPHYLETIC PARALOGS  
    list_of_monophyletic_tuples = []
    list_of_non_monophyletic_tuples = []
    for k,v in monophyly_dict.items():
        if v=="monophyletic":
            list_of_monophyletic_tuples.append(k)
        else:
            list_of_non_monophyletic_tuples.append(k)
             
    paralogs_to_remove = []     # THIS WILL HOLD ALL THE PARALOGS TO BE REMOVED
    # FIRST PROCESS THE MONOPHYLETIC PARALOGS 
    for entry in list_of_monophyletic_tuples:
        paralog_distances = []
        for paralog in entry:
            paralog_distances.append(tree.get_distance(paralog))            # GET DISTANCE OF EVERY PARALOG FROM THE ROOT
            
        distances_dict = dict(zip(entry, paralog_distances))                # CREATE DICTIONARY WITH PARALOGS AND THEIR DISTANCES
        shortest_paralog = min(distances_dict, key=distances_dict.get)      # FIND THE PARALOG WITH THE SHORTEST BRANCH
        del distances_dict[shortest_paralog]                                # AND DELETE IT TO KEEP EVERYTHING ELSE
        paralogs_to_remove += list(distances_dict.keys())                   # THEN ADD WHAT'S LEFT TO THE PARALOGS TO BE REMOVED
        
    # THEN PROCESS THE NON-MONOPHYLETIC PARALOGS
    for entry in list_of_non_monophyletic_tuples:
        for paralog in list(entry):
            paralogs_to_remove.append(paralog)                              # ADD EVERYTHING TO THE PARALOGS TO BE REMOVED
    return paralogs_to_remove

def RemoveSequences(tree, fasta, paralogs_to_remove):
    # READ INPUT FASTA
    with open(fasta, "r") as f:  
        lines = f.readlines()
        stripped = [x.strip() for x in lines]
        new_stripped = [stripped[i:i+2] for i in range(0,len(stripped),2)]
        new_stripped_sorted = sorted(new_stripped, key = lambda x: x[0])
        final_stripped = [item for sublist in new_stripped_sorted for item in sublist]
    # WRITE NEW FASTA WITHOUT THE PARALOGS TO BE REMOVED
    output = open(fasta + ".new", "w")
    for i in range(0, len(final_stripped), 2):   
        if final_stripped[i][1:] not in paralogs_to_remove:
            output.write(final_stripped[i] + "\n")
            output.write(final_stripped[i+1] + "\n")     
