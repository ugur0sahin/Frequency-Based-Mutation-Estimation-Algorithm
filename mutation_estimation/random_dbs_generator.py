##!!!NOTE!!! - This block MUST run for a single time to sustain database consistency every generation create an unique database and overall_mutation_spanset

import numpy as np
import string, random
import pickle


def random_mutation_spanset_generator(number_of_mutation, path_to_export = "C:\Users\Uğur\PycharmProject\Frequency_Based_Mutation_Estimation\dbs"):
    sample_mutation_dict=dict()
    for iterator in range(number_of_mutation):
        letters = string.ascii_letters
        letter_code = ''.join(random.choice(letters.upper()) for i in range(4))

        letters = string.digits.replace("0","1")
        number_code = ''.join(random.choice(letters) for i in range(2))

        for number in number_code:
            position = np.random.randint(len(letter_code))
            letter_code = letter_code[:position] + str(number) + letter_code[position:]
        sample_mutation_dict[iterator]=letter_code



    file_of_overall_span_set_loc =  open(path_to_export + "\\file_of_overall_span_set_loc.pkl","wb")
    pickle.dump(sample_mutation_dict, file_of_overall_span_set_loc)
    file_of_overall_span_set_loc.close()

    print("Sample mutation dbs is created with " + str(number_of_mutation) + " number of different mutations !")

    return sample_mutation_dict

    # generate patient database RANDOM

def random_patient_dbs_generator(number_of_cases, number_of_different_mutations, path_to_export = "C:\Users\Uğur\PycharmProject\Frequency_Based_Mutation_Estimation\dbs"):
    file_of_case_dbs = open(path_to_export + "\\file_of_case_dbs.pkl","wb", )
    patient_dbs = np.random.randint(2, size=(number_of_cases, number_of_different_mutations))
    pickle.dump(patient_dbs, file_of_case_dbs)
    file_of_case_dbs.close()
    return patient_dbs