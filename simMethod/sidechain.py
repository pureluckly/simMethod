
# -*- coding:utf-8 -*-
from rdkit import Chem
import pandas as pd
import csv
from itertools import combinations
def delete_bonds(mol,bonds):
    #use the same parent mol object and create editable mol
    em = Chem.EditableMol(mol)
    #loop through the bonds to delete
    #print "Breaking bonds between atoms: ",bonds
    for b in bonds:
        #remove the bond
        em.RemoveBond(b[0], b[1])
        #now add attachement points
        newAtomA = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[0],newAtomA,Chem.BondType.SINGLE)
        newAtomB = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[1],newAtomB,Chem.BondType.SINGLE)

    #should be able to get away without sanitising mol
    #as the valencies should be okay
    modifiedMol = em.GetMol()

    #do not do a full sanitization, but do find rings and calculate valences:
    Chem.SanitizeMol(modifiedMol,Chem.SanitizeFlags.SANITIZE_PROPERTIES|Chem.SanitizeFlags.SANITIZE_SYMMRINGS)

    fragmented_smi = Chem.MolToSmiles(modifiedMol,True)
    #print fragmented_smi
    return fragmented_smi

smart = Chem.MolFromSmarts('[R]!=!#!@[AR0]')
cycle_smart = Chem.MolFromSmarts('[R]@[R]')
def generate_bondls(mol):
    bond_ls = []
    all_bonds = mol.GetSubstructMatches(smart)
    for bond in all_bonds:
        frag = delete_bonds(mol, [bond])#cut one bond if generate one chain is suit bond
        fragments = frag.split('.')
        for f in fragments:
            fm = Chem.MolFromSmiles(f)
            if fm.HasSubstructMatch(cycle_smart) is False:
                bond_ls.append(bond)
    return bond_ls
def filter_frag(frag_smi):
    out_fragments = []
    fragments = frag_smi.split('.')
    for f in fragments:
        mol_f = Chem.MolFromSmiles(f)
        if mol_f is not None:
            if mol_f.HasSubstructMatch(cycle_smart) is False:
                out_fragments.append(f)
    return out_fragments
def generate_frag(mol_smi):
    mol = Chem.MolFromSmiles(mol_smi)
    if "." in mol_smi:
        frags = list(Chem.GetMolFrags(mol, asMols=True))
        frags.sort(key=lambda x: x.GetNumHeavyAtoms(), reverse=True)
        mol = frags[0]
    if mol.HasSubstructMatch(cycle_smart) is False:#all is chain,return one fragment
        return mol_smi
    else:
        bond_ls = generate_bondls(mol)
        if len(bond_ls) == 0:#all is scaffold,return None
            return None
        else:#chain and scaffold mixed,return sets of fragments
            fragmented_smi = delete_bonds(mol, bond_ls)
            out_fragments = filter_frag(fragmented_smi)
            return ','.join(out_fragments)
# inpath = '/home/ll/Desktop/testmol.csv'
# outpath = '/home/ll/Desktop/tempfrag.csv'
inpath = '/home/ll/work/simMethod/single_smile.csv'
outpath = '/home/ll/work/simMethod/single_smile_sidechain2.csv'
reader = csv.reader(open(inpath))
reader.next()
id_ls = []
frag_ls = []
for line in reader:
    id_ls.append(line[0])
    frag_ls.append(generate_frag(line[1]))
out_df = pd.DataFrame({'chembl_id':id_ls, 'frag':frag_ls})
out_df.to_csv(outpath)

