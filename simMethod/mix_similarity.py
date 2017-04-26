#!/usr/bin/env python
# -*- coding:utf-8 -*-
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import csv, os, math
import multiprocessing as mp
from scipy.stats import spearmanr
scaff_path = '/home/ll/work/simMethod/single_scaff.csv'
frag_path = '/home/ll/work/simMethod/single_smile_sidechain2.csv'

def generate_scaff_dict(scaff_path):
    reader = csv.reader(open(scaff_path))
    scaff_dict = {}
    reader.next()
    for line in reader:
        scaff_dict[line[0]] = [line[1], line[2], line[3]]
        #id, scaff, number of mol heavy atoms, number of scaffold heavy atoms
    return scaff_dict#id, scaff
def generate_chain_dict(frag_path):
    reader = csv.reader(open(frag_path))
    chain_dict = {}
    reader.next()
    for line in reader:
        chain_dict[line[1]] = line[2].split(',')
        #id, frag_list
    return chain_dict
scaff_dict = generate_scaff_dict(scaff_path)
chain_dict = generate_chain_dict(frag_path)
longbits = 16384
fpdict = {}
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
def CalculateFP(fp_name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)
    return fpdict[fp_name](m)
def getBulkSimilarity(fp, fp_list, sort=False):
    if sort is False:
        return DataStructs.BulkTanimotoSimilarity(fp, fp_list)
    else:
        return sorted(DataStructs.BulkTanimotoSimilarity(fp, fp_list), reverse=True)

def calculate_score(ref_id, query_list_id, fpname):
    ref_id = 'CHEMBL' + ref_id
    query_list_id = ['CHEMBL' + i for i in query_list_id]

    ref_scaff = scaff_dict[ref_id][0]#scaff smile
    query_scaff_ls = [scaff_dict[x][0] for x in query_list_id]

    ref_chain = chain_dict[ref_id]#ref mol chain list
    query_chain_ls = [(x, chain_dict[x]) for x in query_list_id]#list of list(query_id, query_chain)
    scaff_sim = np.zeros(4)
    chain_sim = np.zeros(4)
    if ref_scaff is not '':#if ref mol has scaff
        rscaff_fp = CalculateFP(fpname, ref_scaff)#ref scaff fingerprint
        qscaff_fp_ls = [CalculateFP(fpname, q) for q in query_scaff_ls]#list of query mol's scaff fingerprint
        scaff_score_list = getBulkSimilarity(rscaff_fp, qscaff_fp_ls, sort=False)#list of similarity score with ref
        # print scaff_score_list
        for i in range(4):
            qid = query_list_id[i]
            if scaff_dict[qid][-1] is not '':#have no scaffold
                scaff_sim[i] = (float(scaff_dict[qid][-1]) / float(scaff_dict[qid][-2])) * scaff_score_list[i]
                #weight * scaff similarity
    if ref_chain is not '':#if mol has side chain
        rchain_fp_ls = [CalculateFP(fpname, r) for r in ref_chain]#ref's list of chain's fingerprint

        for chain_ls in query_chain_ls:#chain_ls is a query mol's all side chain,which is a list
            score = 0.0
            for q in chain_ls[1]:#q is a chain of a query mol's side chain, is a smile
                print q
                qchain_fp = CalculateFP(fpname, q)#one side chain's fp
                print getBulkSimilarity(qchain_fp, rchain_fp_ls)
                a_chain_score = getBulkSimilarity(qchain_fp, rchain_fp_ls, sort=True)[0]#rank and select maxnumber
                weight = float(Chem.MolFromSmiles(q).GetNumHeavyAtoms()) / float(scaff_dict[chain_ls[0]][-2])
                # print weight, a_chain_score,'--'
                score += (weight * a_chain_score)
            chain_sim[query_chain_ls.index(chain_ls)] = score
    print 'scaff:', scaff_sim
    print 'chain:', chain_sim
    similarity = np.array(scaff_sim) + np.array(chain_sim)
    print 'score:', similarity
    return similarity.tolist()
calculate_score('8614', ['274589','8336','12694','23975'],'lecfp4')
# def getCorrelations(sim_array):
#     ref_corr = range(4, 0, -1)
#     if len(set(sim_array)) == 1:
#         corr = (0, None)
#     else:
#         corr = spearmanr(ref_corr, sim_array)
#     return corr[0]
#
# def worker(file_path, file_ls, result_path):
#     for file in file_ls:
#         inputfile = os.path.join(file_path, file)
#         outputfile = os.path.join(result_path, file)
#         data = np.loadtxt(inputfile, dtype=str)
#         f = open(outputfile, 'a')
#         for i, series in enumerate(data):
#             ref = series[0]
#             query_ls = [series[x] for x in range(1, 5)]
#             similarity = calculate_score(ref, query_ls, 'lecfp4')
#             #similarity is a list of similarity score between query and ref mol
#             corr = getCorrelations(similarity)
#             content = '%s,%s,%s,%s,%s,%s' % (i, corr, similarity[0], similarity[1], similarity[2], similarity[3])
#             print content
#             f.write(content + '\n')
#         f.close()
# def mp_run(nprocs):
#     file_path = '/home/ll/work/work-before/literature-based similarity/similaritybenchmark-master/SingleAssay/dataset/'
#     result_path = '/home/ll/work/simMethod/newsim_result2'
#     file_list = os.listdir(file_path)
#     #print len(file_list)
#     if not os.path.exists(result_path):
#         os.mkdir(result_path)
#
#     chunk_size = int(math.ceil(len(file_list) / float(nprocs)))
#     procs = list()
#     for i in range(nprocs):
#         print i
#         p = mp.Process(
#             target=worker,
#             args=(file_path, file_list[chunk_size * i:chunk_size * (i + 1)], result_path)
#         )
#         procs.append(p)
#         p.start()
#     for p in procs:
#         p.join()
# mp_run(8)
#
#
