#!/usr/bin/env python
# -*- coding:utf-8 -*-
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import Chem
import csv
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import pandas as pd
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

longbits = 16384
fpdict = {}
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
FPcalculator = fpdict['lecfp4']

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
        scaffold_sim *= (scaffold_weight1 + scaffold_weight2)
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
                weight = s_weight1 + s_weight2
        side_chain_sim += sim * weight
        side_chain_weight += weight
    return (scaffold_sim + side_chain_sim) / (scaffold_weight1 + scaffold_weight2 + side_chain_weight)
