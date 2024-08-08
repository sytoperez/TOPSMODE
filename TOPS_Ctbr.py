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
import csv
from numpy import genfromtxt
import pandas as pd
import Calculadora as calc

numero = 15
path = ""
smarts = {}

def load_smarts():
    smarts['Abraham_3'] = calc._ReadPatts_Ab(os.path.join(get_script_path(),'SMARTS/Abraham_3.txt'))
    smarts['Abraham_4'] = calc._ReadPatts_Abalpha(os.path.join(get_script_path(),'SMARTS/Abraham_4.txt'))
    smarts['Crippen'] = calc._ReadPatts(os.path.join(get_script_path(),'SMARTS/Crippen.txt'))
    smarts['Dip'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Dip.txt'))
    smarts['Pols'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pols.txt'))
    smarts['Pol'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Pol.txt'))
    smarts['Std'] = calc._ReadPatts_Pol(os.path.join(get_script_path(),'SMARTS/Std.txt'))

def get_molname(mol, id_field_set, id):
    tmp = ""
    if id_field_set == 'title':
        tmp = mol.GetProp('_Name')
    elif id_field_set == 'gen' or id_field_set is None:
        tmp = 'MolID_' + str(id)
        id += 1
    return tmp, id

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")

def get_script_path():
    return os.path.abspath(os.path.dirname(__file__))

def plot_contributions(moleculas, in_file, type_set, total_only):
    # load source data
    df = pd.read_csv(in_file + ".csv", sep=';')
    df['molecule'] = np.nan 

    data = df.to_numpy()
    lista = np.copy(data)
    headers = list(df.columns.values[2:])

    data = data[:, 2:]
    lista = lista[:, 1]
    c = 0

    # pre-compute drawing parameters
    ps = Draw.ContourParams()
    ps.fillGrid = True
    ps.gridResolution = 0.1
    ps.extraGridPadding = 0.5
    # ps.setScale = False

    # draw the molecules
    for m, mol in enumerate(moleculas, start=1):
        mol_H = Chem.AddHs(mol, explicitOnly=True)
        mol_sinH = Chem.RemoveHs(mol_H, False, True, False)

        locs = []
        if type_set == 'bond':
            cantidad = mol_sinH.GetNumBonds()
            for e, bond in enumerate(mol_sinH.GetBonds(), start=1):
                if e != lista[c + e - 1]:
                    print('Error datos no coinciden. Bond ' + str(e) + '!=' + lista[c + e - 1])
                    return
                idx1 = bond.GetBeginAtomIdx()
                idx2 = bond.GetEndAtomIdx()
                pos1 = mol_sinH.GetConformer().GetAtomPosition(idx1)
                pos2 = mol_sinH.GetConformer().GetAtomPosition(idx2)
                p = (pos1 - pos2)
                centrox = max(pos1.x, pos2.x) - abs(p.x / 2)
                centroy = max(pos1.y, pos2.y) - abs(p.y / 2)
                locs.append(Geometry.Point2D(centrox, centroy))
        else:
            cantidad = mol_sinH.GetNumAtoms()
            for e, atom in enumerate(mol_sinH.GetAtoms(), start=1):
                if e != lista[c + e - 1]:
                    print('Error datos no coinciden. Atom ' + str(e) + '!=' + lista[c + e - 1])
                    return
                pos = mol_sinH.GetConformer().GetAtomPosition(atom.GetIdx())
                locs.append(Geometry.Point2D(pos.x, pos.y))

        sigma = 0.3 * (locs[0] - locs[1]).Length()
        sigmas = [sigma] * cantidad

        if total_only:
            i = len(headers) - 1
            out_file = os.path.join(os.path.dirname(in_file + ".csv"), str(headers[i]) + '_' + str(m) + '.png')
            
            d = Draw.MolDraw2DCairo(600, 600)
            d.ClearDrawing()
            Draw.ContourAndDrawGaussians(d, locs, list(data[c:(c + cantidad), i]), sigmas, nContours=10, params=ps)
            d.drawOptions().clearBackground = False
            d.DrawMolecule(mol_sinH)
            d.WriteDrawingText(out_file)
        else:
            for i in range(0, len(headers)):
                out_file = os.path.join(os.path.dirname(in_file + ".csv"), str(headers[i]) + '_' + str(m) + '.png')

                d = Draw.MolDraw2DCairo(600, 600)
                d.ClearDrawing()
                Draw.ContourAndDrawGaussians(d, locs, list(data[c:(c + cantidad), i]), sigmas, nContours=10, params=ps)
                d.drawOptions().clearBackground = False
                d.DrawMolecule(mol_sinH)
                d.WriteDrawingText(out_file)

        c = c + cantidad

def extraer_ctrb(mol, prop, n, tipo):
    c = 0
    if (tipo=='bond'):
        ea = calc.adj_enlace(mol)
        for enlace in mol.GetBonds():
            ea[c, c] = round(float(enlace.GetProp(prop + str(1))), 6)
            c = c + 1
        multipl = ea
        for j in range(1, n):
            multipl = np.dot(multipl, ea)

        diagonal = multipl.diagonal()
    else:
        ea = np.array(Chem.GetAdjacencyMatrix(mol), "d")
        des = np.empty([n])
        for atomo in mol.GetAtoms():
            ea[c, c] = atomo.GetDoubleProp(prop)
            c = c + 1
        des[0] = np.sum(ea.diagonal())
        multipl = ea
        for j in range(1, n):
            multipl = np.dot(multipl, ea)
        diagonal = multipl.diagonal()
    return diagonal


def calcular_ato(moleculas, variables, coeff, names, linear_set, verbose):
    print('calculating contributions...')

    titulo = ["molecule", "atom", "nombre"] + [p for p in variables] + ["total"]
    tabla_final = np.array([])

    for m, mol in enumerate(moleculas, start = 1):
        if verbose:
            print('\tMolecula %d' % m)

        mol_H = Chem.AddHs(mol, explicitOnly=True)
        mol_sinH = Chem.RemoveHs(mol_H, False, True, False)
        mol_ch_sinH = Chem.RemoveHs(mol_H, False, True, True)
        Chem.SanitizeMol(mol_H)
        Chem.SanitizeMol(mol_sinH)
        at = mol_sinH.GetAtoms()
        bonds = mol_sinH.GetBonds()
        at_ch = mol_ch_sinH.GetAtoms()
        nombres = variables.keys()

        if 'solo' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('solo', 1)

        if 'Std' in nombres:
            order_Std, patts_Std = smarts['Std']
            valores_std = calc._pyGetContribs(mol_sinH, patts_Std, order_Std)  # Standard distances
            for enlace in bonds:
                enlace.SetProp('Std1', str(round(valores_std[enlace.GetIdx()], 6)))
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                if (at_final.HasProp('Std') == 0):
                    at_final.SetDoubleProp('Std', (valores_std[enlace.GetIdx()] / at_final.GetDegree()))
                else:
                    at_final.SetDoubleProp('Std', at_final.GetDoubleProp('Std') + (
                                valores_std[enlace.GetIdx()] / at_final.GetDegree()))
                if (at_inicial.HasProp('Std') == 0):
                    at_inicial.SetDoubleProp('Std', (valores_std[enlace.GetIdx()] / at_inicial.GetDegree()))
                else:
                    at_inicial.SetDoubleProp('Std', at_inicial.GetDoubleProp('Std') + (
                                valores_std[enlace.GetIdx()] / at_inicial.GetDegree()))
                
        if 'Dip' in nombres:
            order_Dip, patts_Dip = smarts['Dip']
            valores_dip = calc._pyGetContribs(mol_sinH, patts_Dip, order_Dip)  # Dipole moments
            for enlace in bonds:
                enlace.SetProp('Dip1', str(round(valores_dip[enlace.GetIdx()], 6)))
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                if (at_final.HasProp('Dip') == 0):
                    at_final.SetDoubleProp('Dip', (valores_dip[enlace.GetIdx()] / at_final.GetDegree()))
                else:
                    at_final.SetDoubleProp('Dip', at_final.GetDoubleProp('Dip') + (
                            valores_dip[enlace.GetIdx()] / at_final.GetDegree()))
                if (at_inicial.HasProp('Dip') == 0):
                    at_inicial.SetDoubleProp('Dip', (valores_dip[enlace.GetIdx()] / at_inicial.GetDegree()))
                else:
                    at_inicial.SetDoubleProp('Dip', at_inicial.GetDoubleProp('Dip') + (
                            valores_dip[enlace.GetIdx()] / at_inicial.GetDegree()))

        if 'Dip2' in nombres:
            valores_dip2 = calc.calcular_dipolos2(mol_sinH)  # Dipole moments2
            for enlace in bonds:
                enlace.SetProp('Dip21', str(round(valores_dip2[enlace.GetIdx()], 6)))
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                if (at_final.HasProp('Dip2') == 0):
                    at_final.SetDoubleProp('Dip2', (valores_dip2[enlace.GetIdx()] / at_final.GetDegree()))
                else:
                    at_final.SetDoubleProp('Dip2', at_final.GetDoubleProp('Dip2') + (
                            valores_dip2[enlace.GetIdx()] / at_final.GetDegree()))
                if (at_inicial.HasProp('Dip2') == 0):
                    at_inicial.SetDoubleProp('Dip2', (valores_dip2[enlace.GetIdx()] / at_inicial.GetDegree()))
                else:
                    at_inicial.SetDoubleProp('Dip2', at_inicial.GetDoubleProp('Dip2') + (
                            valores_dip2[enlace.GetIdx()] / at_inicial.GetDegree()))

        if 'Hyd' or 'Mol' in nombres:
            order, patts = smarts['Crippen']
            valores_hyd_MR = calc._pyGetAtomContribs(mol_sinH, patts, order, verbose)  # Hydrof and MR Crippen
            if 'Hyd' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Hyd', valores_hyd_MR[i][0])

            if 'Mol' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Mol', valores_hyd_MR[i][1])
                    
        if 'Pols' in nombres:
            order_Pols, patts_Pols = smarts['Pols']
            valores_Pols = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pols, order_Pols, verbose)  # Polar surface area
            for i in range(len(at)):
                at[i].SetDoubleProp('Pols', valores_Pols[i])
        
        if 'Van' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Van', Chem.GetPeriodicTable().GetRvdw(at[i].GetAtomicNum()))

        if 'Gas' in nombres:
            # AllChem.ComputeGasteigerCharges(mol_sinH)
            # stored on each atom under the name _GasteigerCharge
            AllChem.ComputeGasteigerCharges(mol_ch_sinH)
            if (np.isnan(at_ch[0].GetDoubleProp('_GasteigerCharge'))):
                # to avoid NaN due Gasteiger rdkit calculation of P. We must have a sdf canonical_smile
                if (mol.HasProp('Original_SMILES')):
                    mol_smi = Chem.MolFromSmiles(mol.GetProp('Original_SMILES'))
                    mol_smi_H = Chem.AddHs(mol_smi, explicitOnly=True)
                    mol_ch_sinH = Chem.RemoveHs(mol_smi_H, False, True, True)
                    AllChem.ComputeGasteigerCharges(mol_ch_sinH)
                    at_ch = mol_ch_sinH.GetAtoms()
            for i in range(len(at)):
                at_ch[i].SetDoubleProp('Gas', at_ch[i].GetDoubleProp('_GasteigerCharge'))

        if 'Ato' in nombres:
            for i in range(len(at)):
                at[i].SetDoubleProp('Ato', at[i].GetMass())

        if 'Pol' in nombres:
            order_Pol, patts_Pol = smarts['Pol']
            valores_Pol = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pol, order_Pol, verbose)  # Polarizabilities
            for i in range(len(at)):
                at[i].SetDoubleProp('Pol', valores_Pol[i])

        if 'Ab-R2' or 'api2' or 'b2h' or 'b2o' or 'l16' in nombres:
            order_Ab, patts_Ab = smarts['Abraham_3']
            valores_Ab = calc._pyGetAtomContribs_Ab2(mol_sinH, patts_Ab, order_Ab, verbose)  # Abraham properties
            if 'Ab-R2' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-R2', valores_Ab[i][0])

            if 'Ab-pi2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-pi2H', valores_Ab[i][1])

            if 'Ab-sumB2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB2H', valores_Ab[i][2])

            if 'Ab-sumB20' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB20', valores_Ab[i][3])

            if 'Ab-logL16' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-logL16', valores_Ab[i][4])

        if 'Ab-sumA2H' in nombres:
            order_Aba, patts_Aba = smarts['Abraham_4']
            valores_Aba = calc._pyGetAtomContribs_Abalpha(mol_sinH, patts_Aba, order_Aba, verbose)  # Abraham properties alpha
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-sumA2H', valores_Aba[i])

        c = 0
        n = names[m - 1]
        columna1 = [str(n)]*mol_sinH.GetNumAtoms()
        columna2 = list(range(1, mol_sinH.GetNumAtoms()+1))
        columna3 = [atomo.GetSymbol() for atomo in mol_sinH.GetAtoms()]
        tabla = np.column_stack((columna1, columna2, columna3))

        if linear_set == 'lin':
            total = np.array([float(coeff[c]) / mol_sinH.GetNumAtoms()]*mol_sinH.GetNumAtoms())
        else:
            total = np.zeros(mol_sinH.GetNumAtoms())

        c += 1
        for p in variables:
            inter = np.zeros(mol_sinH.GetNumAtoms(), dtype=np.int8)

            for idx in variables[p]:
                if p == 'Gas':
                    valores = extraer_ctrb(mol_ch_sinH, p, idx, 'ato')
                elif p == 'solo':
                    valores = np.ones(mol_sinH.GetNumAtoms(), dtype=np.int8)
                else:
                    valores = extraer_ctrb(mol_sinH, p, idx, 'ato')

                # if sum(valores)!=0:
                #    valores = valores / sum(valores) # ver si incluimos standarizacion o os coeff sin std
                v2 = valores * float(coeff[c])
                v2.astype(float)
                inter = inter + v2
                total = total + v2
                c += 1

            tabla = np.hstack((tabla, np.atleast_2d(inter).T))

        tabla = np.hstack((tabla, np.atleast_2d(total).T))
 
        if len(tabla_final) == 0:
            tabla_final = tabla
        else:
            tabla_final = np.vstack((tabla_final, tabla))

    tabla_final = np.vstack((titulo, tabla_final))
    return tabla_final


def calcular_bond(moleculas, variables, coeff, names, linear_set, verbose):
    print('calculating contributions...')

    titulo = ["molecule", "bond"] + [p for p in variables] + ["total"]    
    tabla_final = np.array([])

    for m, mol in enumerate(moleculas, start=1):
        if verbose:
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
        nombres = variables.keys()

        if 'solo' in nombres:
            for enlace in bonds:
                enlace.SetProp('solo', str(1))

        if 'Std' in nombres:
            order_Std, patts_Std = smarts['Std']
            valores_std = calc._pyGetContribs(mol_sinH, patts_Std, order_Std)  # Standard distances
            for enlace in bonds:
                enlace.SetProp('Std1', str(round(valores_std[enlace.GetIdx()], 6)))

        if 'Dip' in nombres:
            order_Dip, patts_Dip = smarts['Dip']
            valores_dip = calc._pyGetContribs(mol_sinH, patts_Dip, order_Dip)  # Dipole moments
            for enlace in bonds:
                enlace.SetProp('Dip1', str(round(valores_dip[enlace.GetIdx()], 6)))

        if 'Dip2' in nombres:
            valores_dip2 = calc.calcular_dipolos2(mol_sinH)  # Dipole moments2
            for enlace in bonds:
                enlace.SetProp('Dip21', str(round(valores_dip2[enlace.GetIdx()], 6)))

        if 'Hyd' or 'Mol' in nombres:
            order, patts = smarts['Crippen']
            valores_hyd_MR = calc._pyGetAtomContribs(mol_sinH, patts, order, verbose)  # Hydrof and MR Crippen
            if 'Hyd' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Hyd', valores_hyd_MR[i][0])

                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor_Hyd = (at_inicial.GetDoubleProp('Hyd') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Hyd') / at_final.GetDegree())
                    enlace.SetProp('Hyd1', str(round(valor_Hyd, 6)))

            if 'Mol' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Mol', valores_hyd_MR[i][1])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor_mol = (at_inicial.GetDoubleProp('Mol') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Mol') / at_final.GetDegree())
                    enlace.SetProp('Mol1', str(round(valor_mol, 6)))


        if 'Pols' in nombres:
            order_Pols, patts_Pols = smarts['Pols']
            valores_Pols = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pols, order_Pols, verbose)  # Polar surface area
            for i in range(len(at)):
                at[i].SetDoubleProp('Pols', valores_Pols[i])
            for enlace in bonds:
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_Pols = (at_inicial.GetDoubleProp('Pols') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Pols') / at_final.GetDegree())
                enlace.SetProp('Pols1', str(round(valor_Pols, 6)))

        if 'Van' in nombres:
            for enlace in bonds:
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_van = (Chem.GetPeriodicTable().GetRvdw(at_inicial.GetAtomicNum()) / at_inicial.GetDegree()) + (
                Chem.GetPeriodicTable().GetRvdw(at_final.GetAtomicNum()) / at_final.GetDegree())
                enlace.SetProp('Van1', str(round(valor_van, 6)))

        if 'Gas' in nombres:
            AllChem.ComputeGasteigerCharges(
                mol_ch_sinH)  # stored on each atom are stored a computed property ( under the name _GasteigerCharge)
            if (np.isnan(at_ch[0].GetDoubleProp('_GasteigerCharge'))):
                mol_smi = Chem.MolFromSmiles(mol.GetProp('Original_SMILES'))
                mol_smi_H = Chem.AddHs(mol_smi, explicitOnly=True)
                mol_ch_sinH = Chem.RemoveHs(mol_smi_H, False, True, True)
                AllChem.ComputeGasteigerCharges(mol_ch_sinH)
                at_ch = mol_ch_sinH.GetAtoms()
                bonds_ch = mol_ch_sinH.GetBonds()
            for enlace in bonds_ch:
                valor_Ch = (at_ch[enlace.GetBeginAtomIdx()].GetDoubleProp('_GasteigerCharge') / at_ch[
                    enlace.GetBeginAtomIdx()].GetDegree()) + (
                                   at_ch[enlace.GetEndAtomIdx()].GetDoubleProp('_GasteigerCharge') / at_ch[
                               enlace.GetEndAtomIdx()].GetDegree())
                enlace.SetProp('Gas1', str(round(valor_Ch, 6)))

        if 'Ato' in nombres:
            for enlace in bonds:
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_mass = (at_inicial.GetMass() / at_inicial.GetDegree()) + (
                        at_final.GetMass() / at_final.GetDegree())
                enlace.SetProp('Ato1', str(round(valor_mass, 6)))

        if 'Pol' in nombres:
            order_Pol, patts_Pol = smarts['Pol']
            valores_Pol = calc._pyGetAtomContribs_Pol(mol_sinH, patts_Pol, order_Pol, verbose)  # Polarizabilities
            for i in range(len(at)):
                at[i].SetDoubleProp('Pol', valores_Pol[i])
            for enlace in bonds:
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor_Pol = (at_inicial.GetDoubleProp('Pol') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Pol') / at_final.GetDegree())
                enlace.SetProp('Pol1', str(round(valor_Pol, 6)))

        if 'Ab-R2' or 'Ab-pi2H' or 'Ab-sumB2H' or 'Ab-sumB20' or 'Ab-logL16' in nombres:
            order_Ab, patts_Ab = smarts['Abraham_3']
            valores_Ab = calc._pyGetAtomContribs_Ab2(mol_sinH, patts_Ab, order_Ab, verbose)  # Abraham properties
            if 'Ab-R2' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-R2', valores_Ab[i][0])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-R2') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-R2') / at_final.GetDegree())
                    enlace.SetProp('Ab-R21', str(round(valor, 6)))

            if 'Ab-pi2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-pi2H', valores_Ab[i][1])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor= (at_inicial.GetDoubleProp('Ab-pi2H') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-pi2H') / at_final.GetDegree())
                    enlace.SetProp('Ab-pi2H1', str(round(valor, 6)))

            if 'Ab-sumB2H' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB2H', valores_Ab[i][2])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor= (at_inicial.GetDoubleProp('Ab-sumB2H') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-sumB2H') / at_final.GetDegree())
                    enlace.SetProp('Ab-sumB2H1', str(round(valor, 6)))

            if 'Ab-sumB20' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-sumB20', valores_Ab[i][3])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-sumB20') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-sumB20') / at_final.GetDegree())
                    enlace.SetProp('Ab-sumB201', str(round(valor, 6)))

            if 'Ab-logL16' in nombres:
                for i in range(len(at)):
                    at[i].SetDoubleProp('Ab-logL16', valores_Ab[i][4])
                for enlace in bonds:
                    at_inicial = at[enlace.GetBeginAtomIdx()]
                    at_final = at[enlace.GetEndAtomIdx()]
                    valor = (at_inicial.GetDoubleProp('Ab-logL16') / at_inicial.GetDegree()) + (
                            at_final.GetDoubleProp('Ab-logL16') / at_final.GetDegree())
                    enlace.SetProp('Ab-logL161', str(round(valor, 6)))

        if 'Ab-sumA2H' in nombres:
            order_Aba, patts_Aba = smarts['Abraham_4']
            valores_Aba = calc._pyGetAtomContribs_Abalpha(mol_sinH, patts_Aba, order_Aba, verbose)  # Abraham properties alpha
            for i in range(len(at)):
                at[i].SetDoubleProp('Ab-sumA2H', valores_Aba[i])
            for enlace in bonds:
                at_inicial = at[enlace.GetBeginAtomIdx()]
                at_final = at[enlace.GetEndAtomIdx()]
                valor = (at_inicial.GetDoubleProp('Ab-sumA2H') / at_inicial.GetDegree()) + (
                        at_final.GetDoubleProp('Ab-sumA2H') / at_final.GetDegree())
                enlace.SetProp('Ab-sumA2H1', str(round(valor, 6)))

        c = 0
        n = names[m - 1]
        columna1 = [str(n)]*mol_sinH.GetNumBonds()
        columna2 = list(range(1, mol_sinH.GetNumBonds()+1))
        tabla = np.column_stack((columna1, columna2))

        if linear_set=='lin':
            total = np.array([float(coeff[c]) / mol_sinH.GetNumBonds()]*mol_sinH.GetNumBonds())
        else:
            total = np.zeros(mol_sinH.GetNumBonds())

        c += 1
        for p in variables:
            inter = np.zeros(mol_sinH.GetNumBonds(), dtype=np.int8)
            for idx in variables[p]:
                if p == 'Gas':
                    valores = extraer_ctrb(mol_ch_sinH, p, idx, 'bond')
                elif p == 'solo':
                    valores = np.ones(mol_sinH.GetNumBonds(), dtype=np.int8)
                else:
                    valores = extraer_ctrb(mol_sinH, p, idx, 'bond')
                #if sum(valores)!=0:
                #    valores = valores / sum(valores) # ver si incluimos standarizacion o os coeff sin std

                v2 = valores * float(coeff[c])
                v2.astype(float)
                inter = inter + v2
                total = total + v2
                c += 1

            tabla = np.hstack((tabla, np.atleast_2d(inter).T))
        tabla = np.hstack((tabla, np.atleast_2d(total).T))

        if len(tabla_final) == 0:
            tabla_final = tabla
        else:
            tabla_final = np.vstack((tabla_final, tabla))

    tabla_final = np.vstack((titulo, tabla_final))
    return tabla_final


def main_params(in_fname, in_model, out_fname, id_field_set, type_set, linear_set, data_only, verbose):
    propiedades = ['Std', 'Dip', 'Dip2', 'Hyd', 'Pols', 'Mol', 'Pol', 'Van', 'Gas', 'Ato', 'Ab-R2', 'Ab-pi2H',
                   'Ab-sumA2H', 'Ab-sumB2H', 'Ab-sumB20', 'Ab-logL16', 'solo']

    # load the input file and ensure the format is correct
    df = pd.read_csv(in_model, sep=';')
    df = df.reindex(['n', 'variables', 'coeff'], axis=1)
    df['n'] = df['n'].astype('int')

    modelos = df.values.tolist()
    for idx, model in enumerate(modelos):
        model.insert(0, f'modelo_{idx+1}')

    # load molecules
    moleculas = Chem.SDMolSupplier(in_fname, False, False, False)
     
    # load the content of the SMARTS files into memory
    load_smarts()

    # calculate molecules' names
    nombres = np.array([])
    id = 1
    for mol in moleculas:
        tmp, id = get_molname(mol, id_field_set, id)
        nombres = np.append(nombres, tmp)

    # build output files' path
    archivo = os.path.split(os.path.abspath(in_fname))
    path = os.path.dirname(os.path.abspath(in_fname)) #+'/'

    """
    if (type_set == 'bond'):
        path = os.path.join(
            path, 
            "TOPSMODE" + archivo[1].split(".")[0], 
            "Contr_", 
            in_fname.rsplit('/', 1)[1].rsplit(".", 1)[0], 
            in_model.rsplit('/', 1)[1].rsplit(".", 1)[0]
        )
    else:
        path = os.path.join(
            path, 
            "TOPSMODE" + archivo[1].split(".")[0], 
            "Contr_ato", 
            in_fname.rsplit('/', 1)[1].rsplit(".", 1)[0], 
            in_model.rsplit('/', 1)[1].rsplit(".", 1)[0]
        )
    """
    path = os.path.join(
        path, 
        "TOPSMODE" + archivo[1].split(".")[0], 
        "Contr_" + type_set, 
        in_fname.rsplit('/', 1)[1].rsplit(".", 1)[0], 
        in_model.rsplit('/', 1)[1].rsplit(".", 1)[0]
    )

    # iterate over the models to process it
    variables = []
    model_counter = 0
    totales = []

    for modelo in modelos:
        """
        modelo[0]: nombre
        modelo[1]: n
        modelo[2]: variables
        modelo[3]: coeff
        """
        model_counter += 1
        model_name = modelo[0]
        n = modelo[1]
        variables = modelo[2].split('|')
        coeff = modelo[3].split('|')
 
        if verbose:
            print('Cargando datos '+ model_name)

        var = []
        des_names = []
        prop_dict = {}
        # FIXME -> para que sirve?
        for i in range(0, n):
            variables[i] = variables[i].replace('(', '|').replace(')', '|').split('|')

            if len(variables[i]) == 1: # u0
                variables[i].append('solo')
                variables[i].append(str(0))
                variables[i][0] = variables[i][0][:-1] # FIXME: xq? eliminar last character

            if not variables[i][1] in prop_dict.keys():
                prop_dict[variables[i][1]] = list()

            prop_dict[variables[i][1]].append(int(variables[i][2]))

        # FIXME -> para que sirve?
        create_dir(path)
        path_model = os.path.join(path, model_name)

        #cumple = os.path.exists(os.path.join(path_model, "TOPSMODE_contr.csv"))        

        if not os.path.exists(os.path.join(path_model, "TOPSMODE_contr.csv")):
            #if variables[0][0]=='u' or variables[0][0]=='u0':
            if type_set == 'bond':
                contribuciones = calcular_bond(moleculas, prop_dict, coeff, nombres, linear_set, verbose=0)  # primera es el intercepto
            else:
                contribuciones = calcular_ato(moleculas, prop_dict, coeff, nombres, linear_set, verbose=0)  # primera es el intercepto

            # export contributions to csv and txt files
            create_dir(path_model)

            df = pd.DataFrame(data=contribuciones[1:], columns=contribuciones[0])
            df.to_csv(os.path.join(path_model, out_fname + ".csv"), index=False, sep=';')
            df.to_csv(os.path.join(path_model, out_fname + ".txt"), index=False, sep=',')

            # aggregate the totals for the late summary
            if len(totales) == 0:
                totales = df
                totales = totales.rename({'total': model_name}, axis=1)
            else:
                totales[model_name] = df['total']

            totales[model_name] = totales[model_name].astype('float64')

            # plot contributions
            if data_only == 'total':
                if type_set == 'ato':
                    plot_contributions(moleculas, os.path.join(path_model, out_fname), 'ato', True)
                else:
                    plot_contributions(moleculas, os.path.join(path_model, out_fname), 'bond', True)

    # calculate the summary of the processed models as the model N+1
    path_summary = os.path.join(path, f"modelo_{model_counter+1}")

    if not os.path.exists(path_summary) and len(modelos) > 1:
        os.makedirs(path_summary)

        # FIXME: esto no debería estar fuera del if? si la carpeta existe no se recalcula el sumatorio

        # remove useless columns for the summary
        to_delete = []
        if 'nombre' in totales.columns:
            to_delete.append('nombre')
        if 'solo' in totales.columns:
            to_delete.append('solo')
        totales = totales.drop(to_delete, axis=1)

        # summarize figures
        totales['total'] = totales.iloc[:, 2:len(list(totales))].sum(axis=1)
        totales['total'] = totales['total']/model_counter
        totales.to_csv(os.path.join(path_summary, out_fname + ".csv"), index=False, header=True, sep=';')
        totales.to_csv(os.path.join(path_summary, out_fname + ".txt"), index=False, header=True, sep=',')

        if data_only == 'total':
            # FIXME: xq no comprueba type_set?
            #if variables[0][0] == 'uato' or variables[0][0] == 'uato0':
            if type_set == 'ato':
                plot_contributions(moleculas, os.path.join(path_summary, out_fname), 'ato', True)
            else:
                plot_contributions(moleculas, os.path.join(path_summary, out_fname), 'bond', True)

        modelos.append(['modelo_'+str(model_counter+1), 'consenso', '-', '-'])

    return nombres


def main():

    parser = argparse.ArgumentParser(description='Calculate TOPSMODE contribution for models.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input file (allowed formats: sdf) with standardized structures')
    parser.add_argument('-m', '--im', metavar='input.csv', required=True,
                        help='File name (with full path) for contributions. Should contain at least these columns (named): "n", "variables" (separated with |), "coeff" (separated with |).'
                             'If the coefficient are of linear model the first one must be the intercept')
    parser.add_argument('-o', '--out', metavar='TOPSMODE.csv', required=True,
                        help='output file with calculated contributions.')
    parser.add_argument('-w', '--id_field_set', metavar='field_set', default='gen',
                       help='name of unique ID for compounds (sdf). gen - auto-generated names and titles - sdf titles will be used')
    parser.add_argument('-t', '--type_set', metavar='bond|ato', default='bond',
                        help='descriptors type: ato for atomic calculation and bond for bond calculation.')
    parser.add_argument('-l', '--linear', metavar='lin|non', default='lin',
                        help='Model technique used (lin-linear, non-nonlinear). Default linear')
    parser.add_argument('-d', '--data_only', metavar='data|total', default='total',
                        help='only output .csv with data without structures in .png')
    parser.add_argument('-v', '--verbose', default=0,
                        help='Integer value. 0 - print no details. 1 and more - verbose output. Default: 0.')

    args = vars(parser.parse_args())
    opt_mix_ordered = None
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "im": in_model = v
        if o == "out": out_fname = v
        if o == "id_field_set": id_field_set = v
        if o == "type_set": type_set = v
        if o == "linear": linear_set = v
        if o == "data_only": data_only = v
        if o == "verbose": verbose = int(v)
    main_params(in_fname, in_model, out_fname, id_field_set, type_set, linear_set, data_only, verbose)


if __name__ == '__main__':
    main()