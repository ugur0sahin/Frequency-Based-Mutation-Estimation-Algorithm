import os
import pickle
import pandas as pd
from collections import Counter

# IMPORT DATA AS PICKLE BOTH CLASSES IN DICT AND PAIR-DICTIONARY
#CCLE DBS IMPORT

parent_dir = os.getcwd()
file_for_CCLE_main_dbs = open(parent_dir + "/dbs/CCLE_Highly_Frequent_dbs.csv")
CCLE_dbs = pd.read_csv(file_for_CCLE_main_dbs)

#LOAD PAIR DICT
f_1, f_2, f_3, f_4 = open("/Users/ugursahin/Downloads/RESULTS/RESULTS-1.pkl", "rb"), open(
    "/Users/ugursahin/Downloads/RESULTS/RESULTS-2.pkl", "rb"), open(
    "/Users/ugursahin/Downloads/RESULTS/RESULTS-3.pkl", "rb"), open(
    "/Users/ugursahin/Downloads/RESULTS/RESULTS-4.pkl", "rb")
dict_1, dict_2, dict_3, dict_4 = pickle.load(f_1), pickle.load(f_2), pickle.load(f_3), pickle.load(f_4)

prob_linkage_file = open(os.getcwd() + "/Prob_Linkage_filtered_rm_duplicates.tsv")
prob_linkage_dbs = pd.read_csv(prob_linkage_file,sep=";")


Grand_dictionary = dict();Grand_dictionary.update(dict_1);Grand_dictionary.update(dict_2);Grand_dictionary.update(dict_3);Grand_dictionary.update(dict_4)

def count_allele_to_gene(allele_dict):
    mutated_genes=list()
    for item in allele_dict.keys():
        mutated_genes.append(item.split(" - ")[0])
    return dict(Counter(mutated_genes))

def sort_dict_by_value(dict_int):
    dict_int = sorted(dict_int.items(), key=lambda x: x[1], reverse=True)

    dict_int_N = dict()
    for dict_item in dict_int:
        dict_int_N[dict_item[0]] = dict_item[1]
    return dict_int_N



def top_N_alleles(prob_linkage_dbs, threshold):
    # GETTING INTERACTED ALLELES TOP X
    object_ls = list()
    for index, row in prob_linkage_dbs.iterrows():
        if index < threshold:
            object_ls.append(row["object_1"].split(" - ")[0]);object_ls.append(row["object_2"].split(" - ")[0])
        else:
            print("\n")
            print(str(row["coefficient"]) + " - This is the threshold coefficient for the Top - " + str(index))
            print("\n")
            break

    return dict(Counter(object_ls))

def sort_frequency_by_gene(CCLE_dbs, Grand_ls_dict):
    gene_dbs_frequent_mut = dict()
    for to_gene_sym in Grand_ls_dict.keys():
        interest_gene = to_gene_sym.split(" - ")[0]
        CCLE_dbs_interested_gene = (CCLE_dbs[CCLE_dbs["Hugo_Symbol"] == interest_gene])
        gene_dbs_frequent_mut[interest_gene]=CCLE_dbs_interested_gene.shape[0]
    return gene_dbs_frequent_mut

if __name__ == '__main__':
    #IN TOP N PROB SCORE WHICH ALLELES THERE ARE BY COUNT
    threshold = 1000
    counted_top_N_dict = top_N_alleles(prob_linkage_dbs, threshold)

    print(sort_dict_by_value(counted_top_N_dict))

if __name__ != '__main__':
    #THIS IS FOR THE OVERALL FRQ WHICH GENES GET MORE MUTATION RATHER THAN OTHERS
    getting_mutation_gene_frequenctly = sort_frequency_by_gene(CCLE_dbs, Grand_dictionary)
    print(sort_dict_by_value(getting_mutation_gene_frequenctly))

