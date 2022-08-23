import pandas as pd
import numpy as np

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

def transform_chosen_profile(span_set_overall, chosen_mutations_as_ls):
    location_set, vector_chosen_mutations = list(), np.zeros(len(span_set_overall.index))
    for mutation in chosen_mutations_as_ls:
        location = span_set_overall.loc[span_set_overall["mutation_column"] == mutation]
        vector_chosen_mutations[location.index-1] = 1
        location_set.append(location)

    span_set_overall.insert(int(1),"Presence",vector_chosen_mutations, True)

    return span_set_overall

def get_frequency_from_case_db(case_db,defined_span_set_overall):
    present_mutation_profile_vector = np.asarray(defined_span_set_overall["Presence"].tolist())
    #matrix vector product part
    result_colomun_as_vector = np.dot(case_db,present_mutation_profile_vector)
    return result_colomun_as_vector, present_mutation_profile_vector

## CORE-1 END

## SHELL-1 INIT

def LinkageD_calculation(mutation_newly_added, chosen_mutations_as_ls, span_set_overall_for_present_profile,span_set_overall_for_newly_added, case_database):

    # Calculation for the obtain presence status of present_profile and newly added.
    defined_span_set_overall_present_profile, defined_span_set_overall_newly_added = transform_chosen_profile(span_set_overall_for_present_profile, chosen_mutations_as_ls), transform_chosen_profile(span_set_overall_for_newly_added, [mutation_newly_added])

    # In this part results_output_state_matrix from case db and present_vector_of_interest is obtained from the get_frequency_from_case function.
    vector_for_case_present_profile, binary_state_of_interest_vector_of_present_profile = get_frequency_from_case_db(case_database,defined_span_set_overall_present_profile)
    vector_for_case_newly_added, binary_state_of_interest_vector_of_newly_added_mutation= get_frequency_from_case_db(case_database, defined_span_set_overall_newly_added)

    #print(vector_for_case_newly_added)
    #print(vector_for_case_present_profile)

    # This part all non-appropriate (under maximum number as a threshold) is replaced with 0 to do calculation logaical based
    vector_for_case_present_profile[vector_for_case_present_profile < np.sum(binary_state_of_interest_vector_of_present_profile)], vector_for_case_newly_added[vector_for_case_newly_added < np.sum(binary_state_of_interest_vector_of_newly_added_mutation)] = 0, 0

    """
    #Calculation of the and/or states to find coocurence 
    Main formula is p(A and B) / p(A)p(B)
    """

    frequency_of_AND_situation_for_two_case = np.logical_and(vector_for_case_newly_added, vector_for_case_present_profile)
    probAandB = np.sum(frequency_of_AND_situation_for_two_case) / frequency_of_AND_situation_for_two_case.size

    probA = np.sum(vector_for_case_newly_added)/vector_for_case_newly_added.size
    probB = np.sum(vector_for_case_present_profile)/vector_for_case_present_profile.size

    probAandB_Indep = probA*probB

    return probAandB - probAandB_Indep


if __name__ == '__main__':
    import pickle

    file_of_mutation_ls = open("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Frequent_Mutation_ls.pkl","rb")
    mutation_ls = pickle.load(file_of_mutation_ls)
    mutation_ls_str, counter, sample_mutation_dict = list(), 0 , dict()
    for cluster in mutation_ls:
        counter += 1
        mutation_ls_str.append(str(cluster[0] + " - " + cluster[1]))
        sample_mutation_dict[counter] = str(cluster[0] + " - " + cluster[1])


    file_of_case_mutation_matrix = open("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Case_Matrix_Frequent.pkl","rb")
    cases_dbs  = pickle.load(file_of_case_mutation_matrix)

    cases_dbs = cases_dbs[::, 1:]

    """
    chosen_mutation_as_ls, mutation_will_get = ["MT-CYB - g.chrM:14798T>C"], "SLC25A40 - g.chr7:87465597_87465598insA"

    span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null = create_mutation_spanset_Dataframe(sample_mutation_dict)

    results_coocured_freq = coocuring_calculation(mutation_will_get, chosen_mutation_as_ls,
                                                  span_set_overall_for_present_profile_Null,
                                                  span_set_overall_for_newly_added_Null, cases_dbs)
    """

    chosen_mutation_as_ls, grand_prop_dict = list(), dict()
    for addin_mutation_to_calculate in list(sample_mutation_dict.values()):
        chosen_mutation_as_ls.append(addin_mutation_to_calculate)
        step_prob_dist = dict()
        for mutation_will_get in sample_mutation_dict.values():
            if mutation_will_get not in chosen_mutation_as_ls:
                span_set_overall_for_present_profile_Null, span_set_overall_for_newly_added_Null = create_mutation_spanset_Dataframe(
                    sample_mutation_dict)
                results_coocured_freq = LinkageD_calculation(mutation_will_get, chosen_mutation_as_ls,
                                                              span_set_overall_for_present_profile_Null,
                                                              span_set_overall_for_newly_added_Null, cases_dbs)
                step_prob_dist[mutation_will_get] = results_coocured_freq

        chosen_mut_vector = np.asarray(list(step_prob_dist.values()))
        np.savetxt(addin_mutation_to_calculate, chosen_mut_vector, delimiter=",")
        grand_prop_dict[addin_mutation_to_calculate] = step_prob_dist
        chosen_mutation_as_ls = list()
        print(grand_prop_dict)

    case_dictionary_pickle_file = open("RESULTS.pkl", "wb")
    pickle.dump(grand_prop_dict, case_dictionary_pickle_file)
    case_dictionary_pickle_file.close()