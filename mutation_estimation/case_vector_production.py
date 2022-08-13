import numpy as np
import pandas as pd

## CORE-1 INIT

def transform_chosen_profile(span_set_overall, chosen_mutations_as_ls):
    location_set, vector_chosen_mutations = list(), np.zeros(len(span_set_overall.index))
    for mutation in chosen_mutations_as_ls:
        location = span_set_overall.loc[span_set_overall["mutation_column"] == mutation]
        vector_chosen_mutations[location.index] = 1
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

def coocuring_calculation(mutation_newly_added, chosen_mutations_as_ls, span_set_overall_for_present_profile,span_set_overall_for_newly_added, case_database):

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
    Main formula is p(A and B) / p(A or B)
    """



    frequency_of_AND_situation_for_two_case = np.logical_and(vector_for_case_newly_added, vector_for_case_present_profile)
    frequency_of_OR_situation_for_two_case = np.logical_or(vector_for_case_newly_added, vector_for_case_present_profile)

    return np.sum(frequency_of_AND_situation_for_two_case)/np.sum(frequency_of_OR_situation_for_two_case)

if __name__ == '__main__':
    """
    #Dataframes should defined to the different templates to prevent overlapping

    span_set_overall_for_present_profile, span_set_overall_for_newly_added = { 0:"A", 1:"B", 2:"C", 3:"D", 4:"E", 5:"F", 6:"G",7:"H",8:"J",9:"K"},\
                                                                             {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H", 8: "J", 9: "K"}
    span_set_overall_for_present_profile, span_set_overall_for_newly_added = pd.DataFrame(span_set_overall_for_present_profile.items(), columns=["index","mutation_column"]), \
                                                                             pd.DataFrame(span_set_overall_for_newly_added.items(), columns=["index", "mutation_column"])
    span_set_overall_for_present_profile, span_set_overall_for_newly_added =span_set_overall_for_present_profile.set_index("index"), \
                                                                            span_set_overall_for_newly_added.set_index("index")

    chosen_mutations_as_ls = ["A","B","F"]

    #defined_span_set_overall = transform_chosen_profile(span_set_overall, chosen_mutations_as_ls)

    #generate patient database RANDOM

    patient_db= np.random.randint(2,size=(100000,10))

    #frequency_of_interest = get_frequency_from_case_db(patient_db,defined_span_set_overall)
    #print(frequency_of_interest)

    print(coocuring_calculation("C", chosen_mutations_as_ls, span_set_overall_for_present_profile,span_set_overall_for_newly_added, patient_db))

    """





