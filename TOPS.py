#==============================================================================
# author          : Alfonso Pérez-Garrido
# date            : 01-10-2023
# version         : 0.1
# python_version  : 3.7
#==============================================================================
import os
import argparse
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from rdkit.Chem import AllChem
import networkx as nx
import pandas as pd
from tkinter import messagebox
import Calculadora as calc

numero = 15
path = ""
descriptores = ['Std', 'Dip', 'Dip2', 'Hyd', 'Pols', 'Mol', 'Pol', 'Van', 'Gas', 'Ato', 'Ab-R2', 'Ab-pi2H', 'Ab-sumA2H', 'Ab-sumB2H', 'Ab-sumB20', 'Ab-logL16']
descriptores_meta = ' '.join(descriptores)
smarts = {}

def load_smarts():
    smarts['Abraham_3'] = calc._ReadPatts_Ab(os.path.join(get_script_path(),'SMARTS/Abraham_3.txt'))
    smarts['Abraham_4'] = calc._ReadPatts_Abalpha(os.path.join(get_script_path(),'SMARTS/Abraham_4.txt'))
    smarts['Crippen'] = calc._ReadPatts(os.path.join(get_script_path(),'SMARTS/Crippen.txt'))
    smarts['Dip'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Dip.txt'))
    smarts['Pols'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pols.txt'))
    smarts['Pol'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pol.txt'))
    smarts['Std'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Std.txt'))

def get_script_path():
    return os.path.abspath(os.path.dirname(__file__))

def get_suffix(des_type):
    if des_type == 'ato':
        return 'uato0'
    return 'u0'

def get_molname(mol, id_field_set, id):
    tmp = ""
    if id_field_set == 'title':
        tmp = mol.GetProp('_Name')
        if tmp == "":
            tmp = 'MolID_' + str(id)
            id += 1
    elif id_field_set == 'gen' or id_field_set is None:
        tmp = 'MolID_' + str(id)
        id += 1
    return tmp, id

def calcular_ato(mol, nombres, numero, verbose):
    if verbose == 1:
        print('calculating descriptors...')

    overall = np.zeros(numero*len(nombres)+1)

    mol_H = Chem.AddHs(mol, explicitOnly=True)
    mol_sinH = Chem.RemoveHs(mol_H, False, True, False)
    mol_ch_sinH = Chem.RemoveHs(mol_H, False, True, True)
    Chem.SanitizeMol(mol_H)
    Chem.SanitizeMol(mol_sinH)
    at = mol_sinH.GetAtoms()
    bonds = mol_sinH.GetBonds()
    at_ch = mol_ch_sinH.GetAtoms()

    fila = np.array([mol_sinH.GetNumAtoms()])

    if 'Std' in nombres:
        order_Std, patts_Std = smarts['Std']
        valores_std = calc._pyGetContribs(mol_sinH, patts_Std, order_Std)  # Standard distances
        for enlace in bonds:
            enlace.SetProp('Std1', str(round(valores_std[enlace.GetIdx()], 6)))
            at_inicial = at[enlace.GetBeginAtomIdx()]
            at_final = at[enlace.GetEndAtomIdx()]
            if (at_final.HasProp('Std')==0):
                at_final.SetDoubleProp('Std',(valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            else:
                at_final.SetDoubleProp('Std', at_final.GetDoubleProp('Std') + (valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            if (at_inicial.HasProp('Std')==0):
                at_inicial.SetDoubleProp('Std',(valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))
            else:
                at_inicial.SetDoubleProp('Std', at_inicial.GetDoubleProp('Std') + (valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))

        des_std = calc.calcular(mol_sinH, 'Std', numero,'ato')
        fila = np.append(fila, des_std)

    if 'Dip' in nombres:
        order_Dip, patts_Dip = smarts['Dip']
        valores_dip = calc._pyGetContribs(mol_sinH, patts_Dip, order_Dip)  # Dipole moments
        for enlace in bonds:
            enlace.SetProp('Dip1', str(round(valores_dip[enlace.GetIdx()], 6)))
            at_inicial = at[enlace.GetBeginAtomIdx()]
            at_final = at[enlace.GetEndAtomIdx()]

            if (at_final.HasProp('Dip')==0):
                at_final.SetDoubleProp('Dip',(valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            else:
                at_final.SetDoubleProp('Dip', at_final.GetDoubleProp('Dip') + (valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            if (at_inicial.HasProp('Dip')==0):
                at_inicial.SetDoubleProp('Dip',(valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))
            else:
                at_inicial.SetDoubleProp('Dip', at_inicial.GetDoubleProp('Dip') + (valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))

        des_dip = calc.calcular(mol_sinH, 'Dip', numero,'ato')
        fila = np.append(fila, des_dip)

    if 'Dip2' in nombres:
        valores_dip2 = calc.calcular_dipolos2(mol_sinH)  # Dipole moments2
        for enlace in bonds:
            enlace.SetProp('Dip21', str(round(valores_dip2[enlace.GetIdx()], 6)))
            at_inicial = at[enlace.GetBeginAtomIdx()]
            at_final = at[enlace.GetEndAtomIdx()]
            if (at_final.HasProp('Dip2') == 0):
                at_final.SetDoubleProp('Dip2', (valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            else:
                at_final.SetDoubleProp('Dip2', at_final.GetDoubleProp('Dip2') + (
                            valores_std[enlace.GetIdx()] / (2*at_final.GetDegree())))
            if (at_inicial.HasProp('Dip2') == 0):
                at_inicial.SetDoubleProp('Dip2', (valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))
            else:
                at_inicial.SetDoubleProp('Dip2', at_inicial.GetDoubleProp('Dip2') + (
                            valores_std[enlace.GetIdx()] / (2*at_inicial.GetDegree())))

        des_dip2 = calc.calcular(mol_sinH, 'Dip2', numero,'ato')
        fila = np.append(fila, des_dip2)

    if 'Hyd' in nombres:
        order, patts = smarts['Crippen']
        valores_hyd_MR = calc._pyGetAtomContribs(mol_sinH, patts, order, verbose)  # Hydrof and MR Crippen
        if 'Hyd' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Hyd', valores_hyd_MR[i][0])
            des_hyd = calc.calcular(mol_sinH, 'Hyd', numero,'ato')
            fila = np.append(fila, des_hyd)

    if 'Pols' in nombres:
        order_Pols, patts_Pols = smarts['Pols']
        valores_Pols = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pols, order_Pols, verbose)  # Polar surface area
        for i in range(len(at)):
            at[i].SetDoubleProp('Pols', valores_Pols[i])
        des_pols = calc.calcular(mol_sinH, 'Pols', numero,'ato')
        fila = np.append(fila, des_pols)

    if 'Mol' in nombres:
        for i in range(len(at)):
            at[i].SetDoubleProp('Mol', valores_hyd_MR[i][1])

        des_mol = calc.calcular(mol_sinH, 'Mol', numero,'ato')
        fila = np.append(fila, des_mol)

    if 'Pol' in nombres:
        order_Pol, patts_Pol = smarts['Pol']
        valores_Pol = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pol, order_Pol, verbose)  # Polarizabilities
        for i in range(len(at)):
            at[i].SetDoubleProp('Pol', valores_Pol[i])

        des_pol = calc.calcular(mol_sinH, 'Pol', numero,'ato')
        fila = np.append(fila, des_pol)

    if 'Van' in nombres:
        for i in range(len(at)):
            at[i].SetDoubleProp('Van', Chem.GetPeriodicTable().GetRvdw(at[i].GetAtomicNum()))

        des_van = calc.calcular(mol_sinH, 'Van', numero,'ato')
        fila = np.append(fila, des_van)

    if 'Gas' in nombres:
        # AllChem.ComputeGasteigerCharges(mol_sinH)
        # stored on each atom are stored a computed property ( under the name _GasteigerCharge)
        AllChem.ComputeGasteigerCharges(mol_ch_sinH)
        if (np.isnan(at_ch[0].GetDoubleProp('_GasteigerCharge'))):
            # debe haber pasado anter por la curacion con KNIME para tener el SMILES original a los P pero no a los organometálicos
            if (mol.HasProp('Original_SMILES')):
                mol_smi = Chem.MolFromSmiles(mol.GetProp('Original_SMILES'))
                mol_smi_H = Chem.AddHs(mol_smi, explicitOnly=True)
                mol_ch_sinH = Chem.RemoveHs(mol_smi_H, False, True, True)
                AllChem.ComputeGasteigerCharges(mol_ch_sinH)
                at_ch = mol_ch_sinH.GetAtoms()
        for i in range(len(at)):
            at_ch[i].SetDoubleProp('Gas', at_ch[i].GetDoubleProp('_GasteigerCharge'))

        des_gas = calc.calcular(mol_ch_sinH, 'Gas', numero,'ato')
        fila = np.append(fila, des_gas)

    if 'Ato' in nombres:
        for i in range(len(at)):
            at[i].SetDoubleProp('Ato', at[i].GetMass())

        des_ato = calc.calcular(mol_sinH, 'Ato', numero,'ato')
        fila = np.append(fila, des_ato)

    if 'Ab-R2' or 'api2' in nombres:
        order_Ab, patts_Ab = smarts['Abraham_3']
        valores_Ab = calc._pyGetAtomContribs_Ab2(mol_sinH, patts_Ab, order_Ab, verbose)  # Abraham properties
        if 'Ab-R2' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-R2', valores_Ab[i][0])

            des = calc.calcular(mol_sinH, 'Ab-R2', numero,'ato')
            fila = np.append(fila, des)

        if 'Ab-pi2H' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-pi2H', valores_Ab[i][1])

            des = calc.calcular(mol_sinH, 'Ab-pi2H', numero,'ato')
            fila = np.append(fila, des)

    if 'Ab-sumA2H' in nombres:
        order_Aba, patts_Aba = smarts['Abraham_4']
        valores_Aba = calc._pyGetAtomContribs_Abalpha(mol_sinH, patts_Aba, order_Aba,
                                                 verbose)  # Abraham properties alpha
        for i in range(len(at)):
            at[i].SetDoubleProp('Ab-sumA2H', valores_Aba[i])

        des = calc.calcular(mol_sinH, 'Ab-sumA2H', numero,'ato')
        fila = np.append(fila, des)

    if 'b2h' or 'b2o' or 'l16' in nombres:
        if 'Ab-sumB2H' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-sumB2H', valores_Ab[i][2])

            des = calc.calcular(mol_sinH, 'Ab-sumB2H', numero,'ato')
            fila = np.append(fila, des)

        if 'Ab-sumB20' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-sumB20', valores_Ab[i][3])

            des = calc.calcular(mol_sinH, 'Ab-sumB20', numero,'ato')
            fila = np.append(fila, des)

        if 'Ab-logL16' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-logL16', valores_Ab[i][4])

            des = calc.calcular(mol_sinH, 'Ab-logL16', numero,'ato')
            fila = np.append(fila, des)

    overall = np.vstack([overall, fila])
    overall = np.delete(overall, 0, 0)

    return overall

def calcular_bond(moleculas, nombres, numero, verbose):
    if verbose == 1:
        print('calculating descriptors...')
    overall = np.array([])
    overall = np.append(overall, [0])
    for n in nombres:
        overall = np.append(overall, [[0 for i in range(1, numero + 1)]])
    m = 0
    for mol in moleculas:
        m = m + 1
        if verbose == 1:
            print('\tMolecula %d' % (m))
        mol_H = Chem.AddHs(mol, explicitOnly=True)
        mol_sinH = Chem.RemoveHs(mol_H, False, True, False)
        mol_ch_sinH = Chem.RemoveHs(mol_H, False, True, True)
        Chem.SanitizeMol(mol_H)
        Chem.SanitizeMol(mol_sinH)
        at = mol_sinH.GetAtoms()
        bonds = mol_sinH.GetBonds()
        at_ch = mol_ch_sinH.GetAtoms()
        bonds_ch = mol_ch_sinH.GetBonds()
        fila = np.array([])
        fila = np.append(fila, [mol_sinH.GetNumBonds()])
        fila_std = np.array([])
        fila_dip = np.array([])
        fila_dip2 = np.array([])
        fila_hyd = np.array([])
        fila_mr = np.array([])
        fila_pls = np.array([])
        fila_vdw = np.array([])
        fila_chg = np.array([])
        fila_atw = np.array([])
        fila_pol = np.array([])
        fila_ar2 = np.array([])
        fila_api2 = np.array([])
        fila_b2h = np.array([])
        fila_b2o = np.array([])
        fila_l16 = np.array([])
        fila_a2h = np.array([])
        if 'Std' in nombres:
            order_Std, patts_Std = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Std.txt'))
            valores_std = calc._pyGetContribs(mol_sinH, patts_Std, order_Std)  # Standard distances
            c = 0
            for enlace in bonds:
                c = c + 1
                enlace.SetProp('Std1', str(round(valores_std[enlace.GetIdx()], 6)))
            des_std = calc.calcular(mol_sinH, 'Std', numero, 'bond')
            fila_std = np.append(fila_std, des_std)

        if 'Dip' in nombres:
            order_Dip, patts_Dip = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Dip.txt'))
            valores_dip = calc._pyGetContribs(mol_sinH, patts_Dip, order_Dip)  # Dipole moments
            c = 0
            for enlace in bonds:
                c = c + 1
                enlace.SetProp('Dip1', str(round(valores_dip[enlace.GetIdx()], 6)))
            des_dip = calc.calcular(mol_sinH, 'Dip', numero,'bond')
            fila_dip = np.append(fila_dip, des_dip)

        if 'Dip2' in nombres:
            valores_dip2 = calc.calcular_dipolos2(mol_sinH)  # Dipole moments2
            c = 0
            for enlace in bonds:
                c = c + 1
                enlace.SetProp('Dip21', str(round(valores_dip2[enlace.GetIdx()], 6)))
            des_dip2 = calc.calcular(mol_sinH, 'Dip2', numero,'bond')
            fila_dip2 = np.append(fila_dip2, des_dip2)

        if 'Hyd' or 'Mol' in nombres:
            order, patts = calc._ReadPatts(os.path.join(get_script_path(),'SMARTS/Crippen.txt'))
            valores_hyd_MR = calc._pyGetAtomContribs(mol_sinH, patts, order, verbose)  # Hydrof and MR Crippen
            if 'Hyd' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Hyd', valores_hyd_MR[i][0])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor_Hyd = (at_inicial.GetDoubleProp('Hyd') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Hyd') / at_final.GetDegree())
                    enlace.SetProp('Hyd1', str(round(valor_Hyd, 6)))
                des_hyd = calc.calcular(mol_sinH, 'Hyd', numero,'bond')
                fila_hyd = np.append(fila_hyd, des_hyd)
            if 'Mol' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Mol', valores_hyd_MR[i][1])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor_mol = (at_inicial.GetDoubleProp('Mol') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Mol') / at_final.GetDegree())
                    enlace.SetProp('Mol1', str(round(valor_mol, 6)))
                des_mol = calc.calcular(mol_sinH, 'Mol', numero,'bond')
                fila_mr = np.append(fila_mr, des_mol)

        if 'Pols' in nombres:
            order_Pols, patts_Pols = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pols.txt'))
            valores_Pols = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pols, order_Pols, verbose)  # Polar surface area
            for i in range(len(at)):
                at[i].SetDoubleProp('Pols', valores_Pols[i])
            c = 0
            for enlace in bonds:
                c = c + 1
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_Pols = (at_inicial.GetDoubleProp('Pols') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Pols') / at_final.GetDegree())
                enlace.SetProp('Pols1', str(round(valor_Pols, 6)))
            des_pols = calc.calcular(mol_sinH, 'Pols', numero,'bond')
            fila_pls = np.append(fila_pls, des_pols)

        if 'Van' in nombres:
            c = 0
            for enlace in bonds:
                c = c + 1
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_van = (Chem.GetPeriodicTable().GetRvdw(at_inicial.GetAtomicNum()) / at_inicial.GetDegree()) + (
                Chem.GetPeriodicTable().GetRvdw(at_final.GetAtomicNum()) / at_final.GetDegree())
                enlace.SetProp('Van1', str(round(valor_van, 6)))
            des_van = calc.calcular(mol_sinH, 'Van', numero,'bond')
            fila_vdw = np.append(fila_vdw, des_van)

        if 'Gas' in nombres:
            #AllChem.ComputeGasteigerCharges(mol_sinH)
            # stored on each atom are stored a computed property ( under the name _GasteigerCharge)
            AllChem.ComputeGasteigerCharges(mol_ch_sinH)
            if np.isnan(at_ch[0].GetDoubleProp('_GasteigerCharge')):
                # debe haber pasado anter por la curacion con KNIME para tener el SMILES original
                mol_smi = Chem.MolFromSmiles(mol.GetProp('Original_SMILES'))
                mol_smi_H = Chem.AddHs(mol_smi, explicitOnly=True)
                mol_ch_sinH = Chem.RemoveHs(mol_smi_H, False, True, True)
                AllChem.ComputeGasteigerCharges(mol_ch_sinH)
                at_ch = mol_ch_sinH.GetAtoms()
                bonds_ch = mol_ch_sinH.GetBonds()
            c = 0
            for enlace in bonds_ch:
                c = c + 1
                valor_Ch = (at_ch[enlace.GetBeginAtomIdx()].GetDoubleProp('_GasteigerCharge') / at_ch[
                    enlace.GetBeginAtomIdx()].GetDegree()) + (
                                   at_ch[enlace.GetEndAtomIdx()].GetDoubleProp('_GasteigerCharge') / at_ch[
                               enlace.GetEndAtomIdx()].GetDegree())
                enlace.SetProp('Gas1', str(round(valor_Ch, 6)))
            des_gas = calc.calcular(mol_ch_sinH, 'Gas', numero, 'bond')
            fila_chg = np.append(fila_chg, des_gas)

        if 'Ato' in nombres:
            c = 0
            for enlace in bonds:
                c = c + 1
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_mass = (at_inicial.GetMass() / at_inicial.GetDegree()) + (
                        at_final.GetMass() / at_final.GetDegree())
                enlace.SetProp('Ato1', str(round(valor_mass, 6)))
            des_ato = calc.calcular(mol_sinH, 'Ato', numero,'bond')
            fila_atw = np.append(fila_atw, des_ato)

        if 'Pol' in nombres:
            order_Pol, patts_Pol = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pol.txt'))
            valores_Pol = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pol, order_Pol, verbose)  # Polarizabilities
            for i in range(len(at)):
                at[i].SetDoubleProp('Pol', valores_Pol[i])
            c = 0
            for enlace in bonds:
                c = c + 1
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_Pol = (at_inicial.GetDoubleProp('Pol') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Pol') / at_final.GetDegree())
                enlace.SetProp('Pol1', str(round(valor_Pol, 6)))
            des_pol = calc.calcular(mol_sinH, 'Pol', numero,'bond')
            fila_pol = np.append(fila_pol, des_pol)

        if 'Ab-R2' or 'api2' or 'b2h' or 'b2o' or 'l16' in nombres:
            order_Ab, patts_Ab = calc._ReadPatts_Ab(os.path.join(get_script_path(),'SMARTS/Abraham_3.txt'))
            valores_Ab = calc._pyGetAtomContribs_Ab2(mol_sinH, patts_Ab, order_Ab, verbose)  # Abraham properties
            if 'Ab-R2' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-R2', valores_Ab[i][0])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-R2') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-R2') / at_final.GetDegree())
                    enlace.SetProp('Ab-R21', str(round(valor, 6)))
                des = calc.calcular(mol_sinH, 'Ab-R2', numero,'bond')
                fila_ar2 = np.append(fila_ar2, des)

            if 'Ab-pi2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-pi2H', valores_Ab[i][1])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor= (at_inicial.GetDoubleProp('Ab-pi2H') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-pi2H') / at_final.GetDegree())
                    enlace.SetProp('Ab-pi2H1', str(round(valor, 6)))
                des = calc.calcular(mol_sinH, 'Ab-pi2H', numero, 'bond')
                fila_api2 = np.append(fila_api2, des)

            if 'Ab-sumB2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB2H', valores_Ab[i][2])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-sumB2H') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-sumB2H') / at_final.GetDegree())
                    enlace.SetProp('Ab-sumB2H1', str(round(valor, 6)))
                des = calc.calcular(mol_sinH, 'Ab-sumB2H', numero, 'bond')
                fila_b2h = np.append(fila_b2h, des)

            if 'Ab-sumB20' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB20', valores_Ab[i][3])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-sumB20') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-sumB20') / at_final.GetDegree())
                    enlace.SetProp('Ab-sumB201', str(round(valor, 6)))
                des = calc.calcular(mol_sinH, 'Ab-sumB20', numero, 'bond')
                fila_b2o = np.append(fila_b2o, des)

            if 'Ab-logL16' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-logL16', valores_Ab[i][4])
                c = 0
                for enlace in bonds:
                    c = c + 1
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-logL16') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-logL16') / at_final.GetDegree())
                    enlace.SetProp('Ab-logL161', str(round(valor, 6)))
                des = calc.calcular(mol_sinH, 'Ab-logL16', numero, 'bond')
                fila_l16 = np.append(fila_l16, des)

        if 'Ab-sumA2H' in nombres:
            order_Aba, patts_Aba = calc._ReadPatts_Abalpha(os.path.join(get_script_path(),'SMARTS/Abraham_4.txt'))
            valores_Aba = calc._pyGetAtomContribs_Abalpha(mol_sinH, patts_Aba, order_Aba, verbose)  # Abraham properties alpha
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-sumA2H', valores_Aba[i])
            c = 0
            for enlace in bonds:
                c = c + 1
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor = (at_inicial.GetDoubleProp('Ab-sumA2H') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Ab-sumA2H') / at_final.GetDegree())
                enlace.SetProp('Ab-sumA2H1', str(round(valor, 6)))
            des = calc.calcular(mol_sinH, 'Ab-sumA2H', numero, 'bond')
            fila_a2h = np.append(fila_a2h, des)

        fila = np.append(fila, fila_std if len(fila_std) > 0 else [])
        fila = np.append(fila, fila_dip if len(fila_dip) > 0 else [])
        fila = np.append(fila, fila_dip2 if len(fila_dip2) > 0 else [])
        fila = np.append(fila, fila_hyd if len(fila_hyd) > 0 else [])
        fila = np.append(fila, fila_pls if len(fila_pls) > 0 else [])
        fila = np.append(fila, fila_mr if len(fila_mr) > 0 else [])
        fila = np.append(fila, fila_pol if len(fila_pol) > 0 else [])
        fila = np.append(fila, fila_vdw if len(fila_vdw) > 0 else [])
        fila = np.append(fila, fila_chg if len(fila_chg) > 0 else [])
        fila = np.append(fila, fila_atw if len(fila_atw) > 0 else [])
        fila = np.append(fila, fila_ar2 if len(fila_ar2) > 0 else [])
        fila = np.append(fila, fila_api2 if len(fila_api2) > 0 else [])
        fila = np.append(fila, fila_a2h if len(fila_a2h) > 0 else [])
        fila = np.append(fila, fila_b2h if len(fila_b2h) > 0 else [])
        fila = np.append(fila, fila_b2o if len(fila_b2o) > 0 else [])
        fila = np.append(fila, fila_l16 if len(fila_l16) > 0 else [])
        overall = np.vstack([overall, fila])

    overall = np.delete(overall, 0, 0)
    return overall

def main_params(in_fname, out_fname,des_names,des_type, numero, id_field_set, activity_field_set, verbose):

    moleculas = Chem.SDMolSupplier(in_fname, False, False, False)

    # setup path and create working directory
    archivo = os.path.split(os.path.abspath(in_fname))
    path = os.path.join(archivo[0], os.path.join("TOPSMODE"+archivo[1].split(".")[0]), "descriptores_"+des_type)

    if not os.path.exists(path):
        os.makedirs(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")

    # build array with column names
    titulo = ['mol']
    if activity_field_set is not None and len(activity_field_set) > 0:
        titulo.append(activity_field_set)
    titulo.append(get_suffix(des_type))

    prefix = 'uato' if des_type == 'ato' else 'u'
    for n in des_names:
        for i in range(1, numero + 1):
            a = prefix + '(' + n + ')' + str(i)
            titulo.append(a)

    # load the content of the SMARTS files into memory
    load_smarts()

    # build data matrix
    id = 1
    nombres = np.array([])
    act = np.array([])
    descriptores = np.array([])

    for mol in moleculas:
        # get the molecule's name and save it for later
        tmp, id = get_molname(mol, id_field_set, id)
        nombres = np.append(nombres, tmp)

        # save activity
        if activity_field_set is not None and activity_field_set != '':
            act = np.vstack([act, mol.GetProp(activity_field_set)])

        # calculate descriptors and store the results in an array
        if des_type == 'ato':
            aux = calcular_ato(mol, des_names, numero, verbose)
        else:
            aux = calcular_bond(mol, des_names, numero, verbose) # TODO

        if len(descriptores) == 0:
            descriptores = aux
        else:
            descriptores = np.vstack((descriptores, aux))

    # FIXME: para que es -> Nombre de las moleculas
    if activity_field_set is not None:
        nombres = np.append(nombres, act, axis=1)

    nombres = nombres.reshape(len(nombres), 1)
    descriptores = np.hstack((nombres, descriptores))

    # transform into dataframe and export
    df = pd.DataFrame(data=descriptores, columns=titulo)
    df.to_csv(os.path.join(path, out_fname), index=False, sep=';')

def main():

    parser = argparse.ArgumentParser(description='Calculate TOPSMODE descriptors for single molecules.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input file (allowed formats: sdf) with standardized structures, '
                             'molecules or reactions should have titles.')
    parser.add_argument('-o', '--out', metavar='TOPSMODE.csv', required=True,
                        help='output file with calculated descriptors.')
    parser.add_argument('-t', '--des_type', metavar='ato|bond', default='bond',
                        help='descriptors type: ato for atomic calculation and bond for bond calculation.')
    parser.add_argument('-d', '--descriptors',
                        metavar=descriptores_meta,
                        default=descriptores, nargs='*',
                        help='descriptors to build. Default all.')
    parser.add_argument('-n', '--num_des', default=15,
                        help='Integer value. Number of descriptors. Default: 15.')
    parser.add_argument('-w', '--id_field_set', metavar='field_set', default='gen',
                       help='name of unique ID for compounds (sdf). gen - auto-generated names and titles - sdf titles will be used'
                            'If omitted for sdf molecule titles will be used or auto-generated names; ')
    parser.add_argument('-a', '--field_activity', metavar='activity_field_set', default=None,
                        help='name of field with activity values for compounds (sdf). none - activity values is ommited')
    parser.add_argument('-v', '--verbose', default=0,
                        help='Integer value. 0 - print no details. 1 and more - verbose output. Default: 0.')

    args = vars(parser.parse_args())
    in_fname = args['in']
    out_fname = args['out']
    des_names = args['descriptors']
    des_type = args['des_type']
    numero = args['num_des']
    id_field_set = args['id_field_set']
    activity_field_set = args['field_activity']
    verbose = args['verbose']

    main_params(in_fname, out_fname, des_names, des_type, numero, id_field_set, activity_field_set, verbose)


if __name__ == '__main__':
    main()
