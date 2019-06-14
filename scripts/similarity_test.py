import argparse
import pandas as pd

def jaccard_coefficient(search_compound, true_compound, search_compounds, true_compounds):
    """ Generate jaccard similarity coefficient between two compounds """
    tempA = search_compounds.loc[search_compounds[search_compound] == 1]
    tempA_index = tempA.index
    tempB = true_compounds.loc[true_compounds[true_compound] == 1]
    tempB_index = tempB.index
    similarity = (len(set(tempA_index).intersection(set(tempB_index))))/(len(set(tempA_index).union(set(tempB_index))))
    return float(similarity)

def similarity_search(search_compounds, true_compounds):
    """ Parallelize generation of jaccard similarity matrix """
    import multiprocessing as mp
    inputs = [(i, j, search_compounds, true_compounds) for i in search_compounds.columns for j in true_compounds.columns]
    pool = mp.Pool(mp.cpu_count())
    outputs = pool.starmap(jaccard_coefficient, inputs)
    pool.terminate()
    pool.join()
    d = {k:{} for k in set([i[0] for i in inputs])}
    for i in range(len(inputs)):
        d[inputs[i][0]].update({inputs[i][1]:outputs[i]})
    jaccard_similarity = pd.DataFrame(data=d)
    
    matched_compounds = pd.DataFrame(columns=['Search compound', 'Matched compounds'])
    for col in jaccard_similarity.columns:
        true_matches = jaccard_similarity[jaccard_similarity[col] == 1.0].index.tolist()
        if len(true_matches) > 0:
            matched_compounds.loc[len(matched_compounds)] = [col, ','.join(true_matches)] 
    jaccard_similarity.drop(matched_compounds['Search compound'], axis=1, inplace=True)
    
    sorted_candidate_compounds = pd.DataFrame(columns=['Search compound', 'Closest compounds', 'Cumulative similarity'])
    for col in jaccard_similarity.columns:
        closest_compounds = list(jaccard_similarity[jaccard_similarity[col]==jaccard_similarity[col].max()].index)
        sum_similarity = jaccard_similarity[col].sum()
        sorted_candidate_compounds.loc[len(sorted_candidate_compounds)] = [col, ','.join(closest_compounds), sum_similarity] 
    sorted_candidate_compounds = sorted_candidate_compounds.sort_values(by='Cumulative similarity', ascending=False)
    
    return matched_compounds, sorted_candidate_compounds

if __name__ == '__main__':
    """ Load two SMILES matrices, generate jaccard similarity coefficients between entries, and format output files  """
    parser = argparse.ArgumentParser(description='Similarity search between two fingerprint matrices.')
    parser.add_argument('-s', type=str, help='path to search set')
    parser.add_argument('-t', type=str, help='path to true set')
    parser.add_argument('-o', type=str, help='output prefix')
    args = parser.parse_args()
    search_compounds = pd.read_csv(args.s, index_col=0)
    true_compounds = pd.read_csv(args.t, index_col=0)
    matched_compounds, sorted_candidate_compounds = similarity_search(search_compounds, true_compounds)
    output_prefix = args.o
    matched_compounds.to_csv('../data/{}_matched_compounds.csv'.format(output_prefix), index=False)
    sorted_candidate_compounds.to_csv('../data/{}_sorted_candidate_compounds.csv'.format(output_prefix), index=False)
