import os

from mutation_estimation.case_vector_production import *

def create_mutation_spanset_Dataframe(sample_mutation_dict):
    span_set_overall_for_present_profile, span_set_overall_for_newly_added = sample_mutation_dict, \
                                                                             sample_mutation_dict

    span_set_overall_for_present_profile, span_set_overall_for_newly_added = pd.DataFrame(
        span_set_overall_for_present_profile.items(), columns=["index", "mutation_column"]), \
                                                                             pd.DataFrame(
                                                                                 span_set_overall_for_newly_added.items(),
                                                                                 columns=["index", "mutation_column"])

    span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null = span_set_overall_for_present_profile.set_index(
        "index"), \
                                                                                       span_set_overall_for_newly_added.set_index(
                                                                                           "index")
    return span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null

def plotting_heat_map(correlation_matrix, allele_lst):
    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    im = ax.imshow(correlation_matrix)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(allele_lst)), labels=allele_lst)
    ax.set_yticks(np.arange(len(allele_lst)), labels=allele_lst)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    ax.set_title("correlation_matrix of coocured Mutations of between the cases")
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    """
    ##!!!NOTE!!! - This block MUST run for a single time to sustain database consistency every generation create an unique database and overall_mutation_spanset
    import numpy as np
    import string, random
    import pickle

    sample_mutation_dict=dict()
    for iterator in range(1000):
        letters = string.ascii_letters
        letter_code = ''.join(random.choice(letters.upper()) for i in range(4))

        letters = string.digits
        number_code = ''.join(random.choice(letters).replace("0","1") for i in range(2))

        for number in number_code:
            position = np.random.randint(len(letter_code))
            letter_code = letter_code[:position] + str(number) + letter_code[position:]
        sample_mutation_dict[iterator]=letter_code

    # generate patient database RANDOM

    patient_db = np.random.randint(2, size=(115, 1000))

    file_of_case_dbs, file_of_overall_span_set_loc = open("dbs/file_of_case_dbs.pkl","wb"), open("dbs/file_of_overall_span_set_loc.pkl","wb")
    pickle.dump(sample_mutation_dict, file_of_overall_span_set_loc), pickle.dump(patient_db, file_of_case_dbs)
    file_of_case_dbs.close(), file_of_overall_span_set_loc.close()
    """

    ##TRIAL-1 Committed Date
    """
    import pickle

    file_of_overall_spanset, file_of_case_dbs = open(os.getcwd() + "/dbs/file_of_overall_span_set_loc.pkl", "rb"), open(os.getcwd() + "/dbs/file_of_case_dbs.pkl", "rb")
    sample_mutation_dict, cases_dbs = pickle.load(file_of_overall_spanset), pickle.load(file_of_case_dbs)

    mutually_inclusive_mutation_set = ["RASF1", "PERK"]
    sample_mutation_dict[len(sample_mutation_dict.keys())], sample_mutation_dict[len(sample_mutation_dict.keys())] = mutually_inclusive_mutation_set[0], mutually_inclusive_mutation_set[1]



    ### Import Mutually Inclusive Mutations


    high_freq_vector_of_both = np.random.randint(10, size=cases_dbs.shape[0])
    high_freq_vector_of_both[ high_freq_vector_of_both < 7] = 0
    high_freq_vector_of_both[high_freq_vector_of_both >= 7] = 1



    cases_dbs = np.column_stack((cases_dbs,high_freq_vector_of_both))
    #Let's little bit differentiate vector for anouther mutation

    new_hf_bot = np.array(high_freq_vector_of_both)
    for i in range(15):
        if high_freq_vector_of_both[np.random.randint(cases_dbs.shape[0])] == 1:
            high_freq_vector_of_both[np.random.randint(cases_dbs.shape[0])] = 0
        else:
            high_freq_vector_of_both[np.random.randint(cases_dbs.shape[0])] = 1


    frequency_of_AND_situation_for_two_case = np.logical_and(high_freq_vector_of_both,
                                                             new_hf_bot)
    frequency_of_OR_situation_for_two_case = np.logical_or(high_freq_vector_of_both, new_hf_bot)


    cases_dbs = np.column_stack((cases_dbs,high_freq_vector_of_both))

    from mutation_estimation.case_vector_production import coocuring_calculation

    chosen_mutation_as_ls = ["32NGWH","2TL5CX","PERK","FCT94Z"]
    step_prob_dist = dict()
    for mutation_will_get in sample_mutation_dict.values():
        if mutation_will_get not in chosen_mutation_as_ls:
            span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null = create_mutation_spanset_Dataframe(sample_mutation_dict)
            results_coocured_freq = coocuring_calculation(mutation_will_get, chosen_mutation_as_ls,span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null, cases_dbs)
            step_prob_dist[mutation_will_get] = results_coocured_freq
    print(step_prob_dist)
    """
    ##TRIAL-2 08.19.22

    import pickle

    file_of_overall_spanset, file_of_case_dbs = open(os.getcwd() + "/dbs/CCLE_Frequent_Mutation_ls.pkl", "rb"), open(os.getcwd() + "/dbs/CCLE_Case_Matrix_Frequent.pkl", "rb")
    sample_mutation_set, cases_dbs_non_binary = pickle.load(file_of_overall_spanset), pickle.load(file_of_case_dbs)


    sample_mutation_dict, counter = dict(), 0
    for cluster in sample_mutation_set:
        sample_mutation_dict[counter] = str(cluster[0]) + " - " +str(cluster[1])
        counter += 1


    case_ls = np.asarray(cases_dbs_non_binary[::,:1],dtype=object).tolist()
    cases_dbs = cases_dbs_non_binary[::,1:]
    from mutation_estimation.case_vector_production import coocuring_calculation


    sample_mutation_set_str = list()
    for cluster in sample_mutation_set:
        sample_mutation_set_str.append(cluster[0] + " - " + cluster[1])
    #print(sample_mutation_set_str)

    #chosen_mutation_as_ls = ["MT-CYB - g.chrM:14798T>C"]

    """
    chosen_mutation_as_ls, grand_prop_dict = list(), dict()
    for addin_mutation_to_calculate in list(sample_mutation_dict.values()):
        chosen_mutation_as_ls.append(addin_mutation_to_calculate)
        step_prob_dist = dict()
        for mutation_will_get in sample_mutation_dict.values():
            if mutation_will_get not in chosen_mutation_as_ls:
                span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null = create_mutation_spanset_Dataframe(
                    sample_mutation_dict)
                results_coocured_freq = coocuring_calculation(mutation_will_get, chosen_mutation_as_ls,
                                                              span_set_overall_for_present_profile_Null,
                                                              span_set_overall_for_newly_added_Null, cases_dbs)
                step_prob_dist[mutation_will_get] = results_coocured_freq

        chosen_mut_vector = np.asarray(list(step_prob_dist.values()))
        np.savetxt(addin_mutation_to_calculate, chosen_mut_vector, delimiter=",")
        grand_prop_dict[addin_mutation_to_calculate] = step_prob_dist
        chosen_mutation_as_ls =list()
        print(grand_prop_dict)

    case_dictionary_pickle_file = open("RESULTS.pkl","wb")
    pickle.dump(grand_prop_dict, case_dictionary_pickle_file)
    case_dictionary_pickle_file.close()
    """

    f_1, f_2, f_3,f_4   = open("/Users/ugursahin/Downloads/RESULTS/RESULTS-1.pkl","rb"), open("/Users/ugursahin/Downloads/RESULTS/RESULTS-2.pkl", "rb"), open("/Users/ugursahin/Downloads/RESULTS/RESULTS-3.pkl", "rb"), open("/Users/ugursahin/Downloads/RESULTS/RESULTS-4.pkl", "rb")

    dict_1, dict_2, dict_3, dict_4 = pickle.load(f_1), pickle.load(f_2), pickle.load(f_3), pickle.load(f_4)

    Grand_dictionary = dict();Grand_dictionary.update(dict_1);Grand_dictionary.update(dict_2);Grand_dictionary.update(dict_3);Grand_dictionary.update(dict_4)

    #print(list(Grand_dictionary.keys()))

    mutation_matrix = np.zeros((len(list(Grand_dictionary.keys())) , len(list(Grand_dictionary.keys()))))
    for i in range(len(list(Grand_dictionary.keys()))):
        for j in range(len(list(Grand_dictionary.keys()))):
            if i != j :
                mutation_matrix[i,j] = Grand_dictionary[list(Grand_dictionary.keys())[i]][list(Grand_dictionary.keys())[j]]
            else:
                mutation_matrix[i, j] = 1
    #print(mutation_matrix)
    np.savetxt("CCLE_Mutual_Cases_on_Frequent_Mutations.csv", mutation_matrix, delimiter=",")


    #ADDING LABELS AS ROW AND COLUMN
    mutation_matrix = np.asarray(mutation_matrix,dtype=object)

    print(Grand_dictionary["BRCA2 - g.chr13:32913079_32913080insA"]["MT-CYB - g.chrM:14798T>C"])
    print(mutation_matrix[4][0])

    """
    mutation_ls = list(Grand_dictionary.keys())
    row_vector_mutation_ls = np.asarray(list(Grand_dictionary.keys()), dtype=object)
    mutation_matrix = np.row_stack((row_vector_mutation_ls, mutation_matrix))

    mutation_ls.insert(0,"Mutation/Mutation")
    #print(mutation_ls)
    colon_vector_mutation_ls = np.asarray(mutation_ls, dtype=object)
    #print(colon_vector_mutation_ls)
    mutation_matrix = np.column_stack((colon_vector_mutation_ls, mutation_matrix))

    #print(mutation_matrix)
    """
    file = open("CCLE_Mutual_Cases_on_Frequent_Mutations.pkl","wb")
    pickle.dump(mutation_matrix, file)


    #print(len(list(Grand_dictionary.values())))
    #print(len(list(sample_mutation_dict.values())))


    import seaborn as sns
    import matplotlib.pyplot as plt


    #print(mutation_matrix.shape)


    mutation_matrix = np.asarray(mutation_matrix, dtype=float)
    df = pd.DataFrame(mutation_matrix, columns=list(Grand_dictionary.keys()), index=list(Grand_dictionary.keys()))

    print(df)

    object_1_ls, object_2_ls, coefficient = list(), list(), list()
    for cluster in Grand_dictionary.items():
        ech = cluster[0]
        list = cluster[1].items()
        for i in list:
            if ech is not i[0]:
                object_1_ls.append(ech), object_2_ls.append(i[0]), coefficient.append(i[1])

    interaction_dict_dbs = { "object_1":object_1_ls , "object_2":object_2_ls , "coefficient":coefficient}
    interaction_dbs = pd.DataFrame(interaction_dict_dbs)
    interaction_dbs.drop_duplicates().to_csv("Prob_Linkage_filtered_rm_duplicates.csv",index=True)

    print(interaction_dbs)


    """
    sns.clustermap(df)
    plt.show()
    """











