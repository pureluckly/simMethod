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

#difine lecfp4 fingerprint
longbits = 16384
fpdict = {}
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
FPcalculator = fpdict['lecfp4']
#generate core and sidechains
def generate_scaffold(mol, framework=False):
    """
    return core_scaffold, side_chains
    """
    core = MurckoScaffold.GetScaffoldForMol(mol)
    side_chains = AllChem.DeleteSubstructs(mol, core)
    # core fp

    if framework:
        core = MurckoScaffold.MakeScaffoldGeneric(core)
    # side chain
    side_mols = Chem.GetMolFrags(side_chains, asMols=True)
    return core, side_mols

def generate_fp(mol, framework=False):
    """
    return ((scaffold_weight, scaffold_fp),[(side_chain_weight, fp)...])
    """
    total_atom_num = float(mol.GetNumHeavyAtoms())
    core_fp = None
    core_weight = 0
    side_fp_list = []
    core, side_mols = generate_scaffold(mol, framework)

    if core:
        core_weight = core.GetNumHeavyAtoms()/total_atom_num
        core_fp = FPcalculator(core)
    for side_mol in side_mols:
        side_fp_list.append((
            side_mol.GetNumHeavyAtoms()/total_atom_num,
            FPcalculator(side_mol)
            ))
    return (core_weight, core_fp), side_fp_list

def similarity(fp1, fp2):
    """
    fp1 and fp2 should in form of the result of generate_fp
    """
    (scaffold_weight1, scaffold_fp1), side_fps1 = fp1
    (scaffold_weight2, scaffold_fp2), side_fps2 = fp2

    if scaffold_fp1 and scaffold_fp2:
        scaffold_sim = DataStructs.FingerprintSimilarity(scaffold_fp1, scaffold_fp2)
        scaffold_sim *= (scaffold_weight1 + scaffold_weight2)/2
    else:
        scaffold_sim = 0
    side_chain_sim = 0
    side_chain_weight = 0
    for s_weight1, s_fp1 in side_fps1:
        weight = 0
        sim = 0
        for s_weight2, s_fp2 in side_fps2:
            tmp_sim = DataStructs.FingerprintSimilarity(s_fp1, s_fp2)
            if tmp_sim > sim:
                sim = tmp_sim
                weight = (s_weight1 + s_weight2)/2
        side_chain_sim += sim * weight
        side_chain_weight += weight
    return scaffold_sim + side_chain_sim

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
            fp1 = generate_fp(mol_dict['CHEMBL' + series[0]])
            fp2 = generate_fp(mol_dict['CHEMBL' + series[1]])
            fp3 = generate_fp(mol_dict['CHEMBL' + series[2]])
            fp4 = generate_fp(mol_dict['CHEMBL' + series[3]])
            fp5 = generate_fp(mol_dict['CHEMBL' + series[4]])
            tanimoto1 = similarity(fp1, fp2)
            tanimoto2 = similarity(fp1, fp3)
            tanimoto3 = similarity(fp1, fp4)
            tanimoto4 = similarity(fp1, fp5)
            sim_ls = [tanimoto1, tanimoto2, tanimoto3, tanimoto4]
            corr = getCorrelations(sim_ls)
            print 'corr:', corr
            content = '%s,%s,%s,%s,%s,%s' % (i, corr, tanimoto1, tanimoto2, tanimoto3, tanimoto4)
            # print content
            f.write(content + '\n')
        f.close()
def mp_run(nprocs):
    file_path = '/home/ll/work/work-before/literature-based similarity/similaritybenchmark-master/SingleAssay/dataset/'
    result_path = '/home/ll/work/simMethod/newsim_result5'
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
mp_run(8)
