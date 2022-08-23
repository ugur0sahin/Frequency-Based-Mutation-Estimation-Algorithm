import os
import pickle
import pandas as pd
from collections import Counter

#DEFINE CLASS
class Mutant_Allele():

    def __init__(self, allele, patient_distribution, total_cases, belonged_gene, position, cases_ls, isDeleterious, isTCGAHs, isCOSMIC, defined_pathways=None):
        self.allele_name = allele
        self.patient_dist = patient_distribution
        self.total_cases = total_cases
        self.cases_ls=cases_ls
        self.defined_pathways=defined_pathways
        self.belonged_gene = belonged_gene
        self.position = position
        self.isDeleterious=isDeleterious
        self.isTCGAhotspot =isTCGAHs
        self.isCOSMIChotspot=isCOSMIC

    def proportion_distribution(self):
        dict_def = (dict(self.patient_dist))
        DF_sum = self.total_cases

        DF_dict = dict()
        for row in list(dict_def.items()):
            DF_dict[row[0]] = row[1] / DF_sum
        return DF_dict

    def get_cooccuring_prob(self,pair_allele,coocured_dictionary):
        prob_link = coocured_dictionary[pair_allele][self.allele_name]
        return prob_link

    def get_prob_vector_with_other_partners(self,coocured_dictionary):
        return coocured_dictionary[self.allele_name]

    def __str__(self):
        return str(self.allele_name)

parent_dir = os.getcwd()
main_CCLE_dbs = pd.read_csv(parent_dir +"/dbs/CCLE_Highly_Frequent_dbs.csv")

patient_diagnosis_database_file = open(parent_dir + "/dbs/DepMap-2018q3-celllines.csv")
main_DepMap_ID_dbs = pd.read_csv(patient_diagnosis_database_file)

objected_alleles = dict()
if __name__ != '__main__': #GENERATING OBJECT PKL FILE
    # IMPORT DATA AS PICKLE

    f_1, f_2, f_3, f_4 = open("/Users/ugursahin/Downloads/RESULTS/RESULTS-1.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-2.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-3.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-4.pkl", "rb")
    dict_1, dict_2, dict_3, dict_4 = pickle.load(f_1), pickle.load(f_2), pickle.load(f_3), pickle.load(f_4)
    Grand_dictionary = dict();Grand_dictionary.update(dict_1);Grand_dictionary.update(dict_2);Grand_dictionary.update(dict_3);Grand_dictionary.update(dict_4)

    for allele_name in Grand_dictionary.keys():
        position = allele_name.split(" - ")[1]


        main_DepMap_ID_dict = dict()
        for index, row in main_DepMap_ID_dbs.iterrows():
            main_DepMap_ID_dict[row["Broad_ID"]] = row["Primary Disease"]
        #main_DepMap_ID_dict it is the main case:diagnosis dict overall

        allele_belonged_database = main_CCLE_dbs[main_CCLE_dbs["Genome_Change"] == position]
        deleterious, cosmic_hs, tcga_hs = allele_belonged_database["isDeleterious"].to_list()[0], allele_belonged_database["isTCGAhotspot"].to_list()[0], allele_belonged_database["isCOSMIChotspot"].to_list()[0]

        cases_belonged_allele = allele_belonged_database["DepMap_ID"].to_list()

        list_of_diagnosis_for_interested_allele = list()
        for case_of_interest in cases_belonged_allele:
            try:
                diagnosis_of_int = main_DepMap_ID_dict[case_of_interest]
            except:
                diagnosis_of_int = "unkown_diagnosis"
            list_of_diagnosis_for_interested_allele.append(diagnosis_of_int)

        diagnosis_dist = dict(Counter(list_of_diagnosis_for_interested_allele))

        sum_of_total_patient_found = sum([x for x in diagnosis_dist.values()])

        objected_alleles[allele_name] = Mutant_Allele(allele_name,diagnosis_dist,sum_of_total_patient_found,allele_name.split(" - ")[0],allele_name.split(" - ")[1],cases_belonged_allele,deleterious,tcga_hs,cosmic_hs)
    file_objected_ls = open("Mutation_Info_as_Obj.pkl","wb")
    pickle.dump(objected_alleles,file_objected_ls)
    file_objected_ls.close()



if __name__ == '__main__': #TESTING OUTPUT PKL OBJECT FILE
    # IMPORT DATA AS PICKLE

    f_1, f_2, f_3, f_4 = open("/Users/ugursahin/Downloads/RESULTS/RESULTS-1.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-2.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-3.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-4.pkl", "rb")
    dict_1, dict_2, dict_3, dict_4 = pickle.load(f_1), pickle.load(f_2), pickle.load(f_3), pickle.load(f_4)

    Grand_dictionary = dict();Grand_dictionary.update(dict_1);Grand_dictionary.update(dict_2);Grand_dictionary.update(dict_3);Grand_dictionary.update(dict_4)


    file_obj_ls = open("Mutation_Info_as_Obj.pkl","rb")
    obj_dict = pickle.load(file_obj_ls)

    chosen_obj = obj_dict["MAGEA2 - g.chrX:151919147C>G"]
    print(chosen_obj.patient_dist)

if __name__ != '__main__':
    #IMPORT DATA

    f_1, f_2, f_3, f_4 = open("/Users/ugursahin/Downloads/RESULTS/RESULTS-1.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-2.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-3.pkl", "rb"), open(
        "/Users/ugursahin/Downloads/RESULTS/RESULTS-4.pkl", "rb")
    dict_1, dict_2, dict_3, dict_4 = pickle.load(f_1), pickle.load(f_2), pickle.load(f_3), pickle.load(f_4)

    Grand_dictionary = dict();Grand_dictionary.update(dict_1);Grand_dictionary.update(dict_2);Grand_dictionary.update(dict_3);Grand_dictionary.update(dict_4)

    parent_dir = os.getcwd()
    file_for_CCLE_main_dbs = open(parent_dir + "/dbs/CCLE_Highly_Frequent_dbs.csv")
    CCLE_dbs = pd.read_csv(file_for_CCLE_main_dbs)

    #CODE INIT
    CCLE_dbs_chosen_gene = CCLE_dbs[CCLE_dbs["Hugo_Symbol"] == "TP53"]
    total_allele_set_for_chosen_gene = list(set(CCLE_dbs_chosen_gene["Genome_Change"].to_list()))

    for allele in total_allele_set_for_chosen_gene:
        objected_dict = Grand_dictionary["TP53 - "+str(allele)]

        from hand_checking import sort_dict_by_value, count_allele_to_gene

        objected_dict_sorted = sort_dict_by_value(objected_dict)
        print('\n')
        print(sort_dict_by_value(count_allele_to_gene(objected_dict_sorted)))
        print("\n----------Cooccuring Dict------------\n")
        for pair_allele in objected_dict_sorted.items():
            print(pair_allele)

        print("_____________________________________________________________________________")
