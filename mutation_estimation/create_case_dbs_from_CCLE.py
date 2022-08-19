import pandas as pd
import numpy as np
import pickle as pkl

def get_filter_CCLE_dbs(file):
    CCLE_file = open(file,"r")
    CCLE_db = pd.read_csv(CCLE_file)

    #Filtering Useless Columns
    del CCLE_db["RNAseq_AC"]; del CCLE_db["SangerWES_AC"]; del CCLE_db["WGS_AC"]
    del CCLE_db["NCBI_Build"]; del CCLE_db["CGA_WES_AC"]; del CCLE_db["HC_AC"]
    del CCLE_db["RD_AC"]; del CCLE_db["ExAC_AF"]; del CCLE_db["COSMIChsCnt"]

    print("Filtered database is constructed !")
    return CCLE_db

def implement_patient_matrix(pd_dbs, dir_extract_overall_mutation_ls ="/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Frequent_Mutation_ls.pkl", dir_extract_case_mutation_matrix="/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Case_Matrix_Frequent.pkl"):

    overall_mutations = list()
    for index, row in pd_dbs.iterrows():
        overall_mutations.append((str(row.Hugo_Symbol) , str(row.Genome_Change)))

    overall_mutations = list(set(overall_mutations))
    print("All mutations uniquely indicated and transerred to list -> " + str(overall_mutations))

    #Fulfill the patient Matrix.
    case_ls = list(set(pd_dbs["DepMap_ID"].tolist()))
    pd_dbs = pd_dbs.set_index("Hugo_Symbol")

    mutation_matrix = np.zeros((len(case_ls), len(overall_mutations)))
    mutation_matrix =np.asarray(mutation_matrix, dtype=object)
    case_vector = np.asarray(case_ls,dtype=object)

    case_mutation_matrix = np.column_stack((case_vector,mutation_matrix))


    for case_index in range(len(case_ls)):
        for gene_index, index_info in pd_dbs.iterrows():
            if case_ls[case_index] == index_info["DepMap_ID"]:
                index_of_mutation = overall_mutations.index((gene_index,index_info["Genome_Change"]))
                case_mutation_matrix[case_index,index_of_mutation+1] = 1
        print("The Case " + str(case_ls[case_index]) + " is implemented to the Matrix ! ")


    case_mutation_matrix_dump_gen_file = open(dir_extract_case_mutation_matrix, "wb")
    pkl.dump(case_mutation_matrix,case_mutation_matrix_dump_gen_file)
    case_mutation_matrix_dump_gen_file.close()
    print("Case Mutation Matrix is written to path " + str(dir_extract_case_mutation_matrix) +" ! ")

    mutation_ls_dump_gen_file = open(dir_extract_overall_mutation_ls, "wb")
    pkl.dump(overall_mutations,mutation_ls_dump_gen_file)
    mutation_ls_dump_gen_file.close()
    print("Overall Mutation Vector is written to path " + str(dir_extract_overall_mutation_ls) + " ! ")

    """
    file_read = open("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Case_Matrix_Frequent.pkl", "rb")
    case_mutation_matrix_read = pkl.load(file_read)
    print(case_mutation_matrix_read)
    """

"""
if __name__ == '__main__':
    #CCLE_dbs_ps = get_filter_CCLE_dbs("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Highly_Frequent_dbs.csv")
    #implement_patient_matrix(CCLE_dbs_ps)

    f = open("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Frequent_Mutation_ls.pkl","rb")
    print(pkl.load(f))

    f_1 = open("/Users/ugursahin/Downloads/Frequency-Based-Mutation-Estimation-Algorithm-main/dbs/CCLE_Case_Matrix_Frequent.pkl","rb")
    print(pkl.load(f_1))
"""