#!/usr/bin/env python
# -*- coding:utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np
import csv, os, math
import multiprocessing as mp
from rdkit import DataStructs
from scipy.stats import spearmanr
path = '/home/ll/work/simMethod/single_smile.csv'
def get_mol(path):
    reader = csv.reader(open(path))
    reader.next()
    mol_dict = {}
    for line in reader:
        # print line[0], line[1]
        mol_dict[line[0]] = Chem.MolFromSmiles(line[1])
    return mol_dict
mol_dict = get_mol(path)
print 'mol_dict is ok!'
#generate core and sidechains
def scaffold(mol):
    core = MurckoScaffold.GetScaffoldForMol(mol)
    # framework = MurckoScaffold.MakeScaffoldGeneric(core)
    side_chains = AllChem.DeleteSubstructs(mol, core)
    return core, side_chains

#difine lecfp4 fingerprint
longbits = 16384
fpdict = {}
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
FPcalculator = fpdict['lecfp4']
def getBulkSimilarity(fp, fp_list, sort=False):
    if sort is False:
        return DataStructs.BulkTanimotoSimilarity(fp, fp_list)
    else:
        return sorted(DataStructs.BulkTanimotoSimilarity(fp, fp_list), reverse=True)
def similarity(ref_id, query_id_ls):
    ref = mol_dict['CHEMBL' + ref_id]#ref molecule
    query_id_ls = ['CHEMBL' + query_id for query_id in query_id_ls]
    query_ls = [mol_dict[i] for i in query_id_ls]#query mol list

    # print query
    ref_core, ref_chains = scaffold(ref)
    core_chain_ls = [scaffold(q) for q in query_ls]
    # print core_chain_ls
    scaff_sim = np.zeros(4)
    chain_sim = np.zeros(4)
    if ref_core is not '':
        ref_core_fp = FPcalculator(ref_core)
        query_core_fp = [FPcalculator(query_core[0]) for query_core in core_chain_ls]
        tanimoto = getBulkSimilarity(ref_core_fp, query_core_fp)

        for i in range(4):
            scaff_sim[i] = (float(core_chain_ls[i][0].GetNumHeavyAtoms()) / query_ls[i].GetNumHeavyAtoms()) * tanimoto[i]

    if ref_chains is not '':
        rchain_fp_ls = [FPcalculator(Chem.MolFromSmiles(r)) for r in Chem.MolToSmiles(ref_chains).split('.')]#ref's list of chain's fingerprint

        for q in range(4):
            qsim = 0.0
            q_chains = Chem.MolToSmiles(core_chain_ls[q][1]).split('.')
            q_chains = [Chem.MolFromSmiles(i) for i in q_chains]#a query mol's fragmol list

            chainfp_ls = [FPcalculator(chain) for chain in q_chains]#q query mol's fragfp list

            for chainfp in chainfp_ls:
                sim = getBulkSimilarity(chainfp, rchain_fp_ls, sort=True)[0]

                weight = float(q_chains[chainfp_ls.index(chainfp)].GetNumHeavyAtoms()) / query_ls[q].GetNumHeavyAtoms()
                qsim += (weight * sim)
            chain_sim[q] = qsim

    similarity = np.array(scaff_sim) + np.array(chain_sim)
    # print 'score:', similarity,type(similarity)
    return similarity.tolist()
# similarity('8614', ['274589','8336','12694','23975'])
def getCorrelations(sim_array):
    ref_corr = range(4, 0, -1)
    if len(set(sim_array)) == 1:
        corr = (0, None)
    else:
        corr = spearmanr(ref_corr, sim_array)
    return corr[0]

def worker(file_path, file_ls, result_path):
    for file in file_ls:
        inputfile = os.path.join(file_path, file)
        outputfile = os.path.join(result_path, file)
        data = np.loadtxt(inputfile, dtype=str)
        f = open(outputfile, 'a')
        for i, series in enumerate(data):
            ref_id = series[0]
            query_id_ls = [series[x] for x in range(1, 5)]

            tanimoto = similarity(ref_id, query_id_ls)

            corr = getCorrelations(tanimoto)
            content = '%s,%s,%s,%s,%s,%s' % (i, corr, tanimoto[0], tanimoto[1], tanimoto[2], tanimoto[3])
            print 'corr:', corr
            f.write(content + '\n')
        f.close()
def mp_run(nprocs):
    file_path = '/home/ll/work/work-before/literature-based similarity/similaritybenchmark-master/SingleAssay/dataset/'
    result_path = '/home/ll/work/simMethod/newsim_result2'
    file_list = os.listdir(file_path)
    #print len(file_list)
    if not os.path.exists(result_path):
        os.mkdir(result_path)

    chunk_size = int(math.ceil(len(file_list) / float(nprocs)))
    procs = list()
    for i in range(nprocs):
        print i
        p = mp.Process(
            target=worker,
            args=(file_path, file_list[chunk_size * i:chunk_size * (i + 1)], result_path)
        )
        procs.append(p)
        p.start()
    for p in procs:
        p.join()
mp_run(8)
