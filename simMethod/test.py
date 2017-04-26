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
        mol_dict[line[0]] = Chem.MolFromSmiles(line[1])# mol_dict{chembl_id: mol}
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
def ref_fp(ref_id):
    ref = mol_dict['CHEMBL' + ref_id]
    ref_core, ref_chains = scaffold(ref)
    ref_core_fp = FPcalculator(ref_core)
    rchain_fp_ls = [FPcalculator(Chem.MolFromSmiles(r)) for r in Chem.MolToSmiles(ref_chains).split('.')]
    return ref_core, ref_chains, ref_core_fp, rchain_fp_ls
def similarity(ref, query_id):
    #ref molecule
    ref_core, ref_core_fp,ref_chains,rchain_fp_ls = ref
    query = mol_dict['CHEMBL' + query_id]# query molecule
    query_core, query_chains = scaffold(query)
    scaff_sim = 0.0
    chain_sim = 0.0
    if Chem.MolToSmiles(ref_core) is not '': #if reference mol has core

        query_core_fp = FPcalculator(query_core)
        tanimoto = DataStructs.FingerprintSimilarity(ref_core_fp, query_core_fp)
        scaff_sim = (float(query_core.GetNumHeavyAtoms()) / query.GetNumHeavyAtoms()) * tanimoto
    if Chem.MolToSmiles(ref_chains) is not '': #if reference mol has chains

        q_chains = Chem.MolToSmiles(query_chains).split('.')
        q_chains = map(Chem.MolFromSmiles, q_chains)
        for chain in q_chains:# a chain from a query mol
            chainfp = FPcalculator(chain)
            temp = []
            for r in rchain_fp_ls: # a chain from a reference mol
                sim = DataStructs.FingerprintSimilarity(r, chainfp)
                temp.append(sim)
            weight = float(chain.GetNumHeavyAtoms()) / query.GetNumHeavyAtoms()
            chain_sim += (weight * max(temp))
    # print 'scaff:', scaff_sim
    # print 'chain:', chain_sim
    similarity = scaff_sim + chain_sim
    # print 'score:', similarity
    return similarity

# similarity('8614', '23975')
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

            ref = ref_fp(series[0])
            tanimoto1 = similarity(ref, series[1])
            tanimoto2 = similarity(ref, series[2])
            tanimoto3 = similarity(ref, series[3])
            tanimoto4 = similarity(ref, series[4])
            sim_ls = [tanimoto1, tanimoto2, tanimoto3, tanimoto4]
            corr = getCorrelations(sim_ls)
            print 'corr:', corr
            content = '%s,%s,%s,%s,%s,%s' % (i, corr, tanimoto1, tanimoto2, tanimoto3, tanimoto4)
            # print content
            f.write(content + '\n')
        f.close()
def mp_run(nprocs):
    file_path = '/home/ll/work/work-before/literature-based similarity/similaritybenchmark-master/SingleAssay/dataset/'
    result_path = '/home/ll/work/simMethod/newsim_result3'
    file_list = os.listdir(file_path)[:8]
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
# mp_run(8)
