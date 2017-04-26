#!/usr/bin/env python
# -*- coding:utf-8 -*-
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
import csv
from rdkit.Chem import DataStructs
import pandas as pd
import numpy as np
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
def delete_bonds(mol,bonds):
    em = Chem.EditableMol(mol)
    for b in bonds:
        em.RemoveBond(b[0], b[1])
        newAtomA = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[0],newAtomA,Chem.BondType.SINGLE)
        newAtomB = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[1],newAtomB,Chem.BondType.SINGLE)
    modifiedMol = em.GetMol()
    Chem.SanitizeMol(modifiedMol,Chem.SanitizeFlags.SANITIZE_PROPERTIES|Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
    fragmented_smi = Chem.MolToSmiles(modifiedMol,True)
    return fragmented_smi
smart_b0_b1 = Chem.MolFromSmarts('[R]!@!=!#[AR0]')
def cut_core(core):
    bond_ls = []
    all_bonds = core.GetSubstructMatches(smart_b0_b1)
    if all_bonds:
        fragment = delete_bonds(core, all_bonds)
        fragments = fragment.split('.')
        molfrag = [Chem.MolFromSmiles(f) for f in fragments]
        return molfrag
    else:
        return core
inputfile = '/home/ll/work/work-before/literature-based similarity/similaritybenchmark-master/SingleAssay/dataset/999.txt'
data = np.loadtxt(inputfile, dtype=str)
