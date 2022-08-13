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
    import pickle

    file_of_overall_spanset, file_of_case_dbs = open(os.getcwd() + "\\dbs\\file_of_overall_span_set_loc.pkl", "rb"), open(os.getcwd() + "\\dbs\\file_of_case_dbs.pkl", "rb")
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








