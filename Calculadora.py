#==============================================================================
# author          : Alfonso Pérez-Garrido
# date            : 01-10-2023
# version         : 0.1
# python_version  : 3.7
#==============================================================================

from rdkit import Chem
import numpy as np

def _ReadPatts_Ab(fileName):
  """ *Internal Use Only*
    parses the pattern list from the data file
  """
  patts = {}
  order = []
  for line in open(fileName,'r').readlines():
    if line[0] != '#':
      splitLine = line.split('\t')
      if len(splitLine) >= 8 and splitLine[0] != '':
        sma = splitLine[1]
        if sma!='SMARTS':
          sma.replace('"', '')
          try:
            p = Chem.MolFromSmarts(sma)
          except:
            pass
          else:
            if p:
              if len(splitLine[0])>1 and splitLine[0][1] not in 'S0123456789':
                cha = splitLine[0][:2]
              else:
                cha = splitLine[0][0]
              if splitLine[2] != '':
                r2  = float(splitLine[2])
              else:
                r2  = 0.0
              if splitLine[3] != '':
                pi2h = float(splitLine[3])
              else:
                pi2h = 0.0
              if splitLine[4] != '':
                b2h = float(splitLine[4])
              else:
                b2h = 0.0
              if splitLine[5] != '':
                b2o = float(splitLine[5])
              else:
                b2o = 0.0
              if splitLine[6] != '':
                l16 = float(splitLine[6])
              else:
                l16 = 0.0

              if splitLine[7] != '':
                tipo = float(splitLine[7])
              else:
                tipo = 0.0

              if cha not in order:
                order.append(cha)

              l = patts.get(cha,[])
              l.append((sma,p,r2, pi2h, b2h, b2o, l16,tipo))
              patts[cha] = l
            else:
              print ('Problems parsing smarts: %s' %(sma))
  return order,patts


def _ReadPatts_Abalpha(fileName):
  """ *Internal Use Only*
    parses the pattern list from the data file
  """
  patts = {}
  order = []
  for line in open(fileName,'r').readlines():
    if line[0] != '#':
      splitLine = line.split('\t')
      if len(splitLine)>=4 and splitLine[0] != '':
        sma = splitLine[1]
        if sma!='SMARTS':
          sma.replace('"','')
          try:
            p = Chem.MolFromSmarts(sma)
          except:
            pass
          else:
            if p:
              if len(splitLine[0])>1 and splitLine[0][1] not in 'S0123456789':
                cha = splitLine[0][:2]
              else:
                cha = splitLine[0][0]
              if splitLine[2] != '':
                a2 = float(splitLine[2])
              else:
                a2 = 0.0
              if splitLine[3] != '':
                tipo = float(splitLine[3])
              else:
                tipo = 0.0

              if cha not in order:
                order.append(cha)

              l = patts.get(cha,[])
              l.append((sma,p,a2,tipo))
              patts[cha] = l
            else:
              print ('Problems parsing smarts: %s' %(sma))
  return order,patts


def _ReadPatts_Pol(fileName):
  """ *Internal Use Only*
    parses the pattern list from the data file
  """
  patts = {}
  order = []
  for line in open(fileName,'r').readlines():
    if line[0] != '#':
      splitLine = line.split('\t')
      if len(splitLine)>=3 and splitLine[0] != '':
        sma = splitLine[1]
        if sma!='SMARTS':
          sma.replace('"','')
          try:
            p = Chem.MolFromSmarts(sma)
          except:
            pass
          else:
            if p:
              if len(splitLine[0])>1 and splitLine[0][1] not in 'S0123456789':
                cha = splitLine[0][:2]
              else:
                cha = splitLine[0][0]
              Pol = float(splitLine[2])
              if cha not in order:
                order.append(cha)
              l = patts.get(cha, [])
              l.append((sma, p, Pol))
              patts[cha] = l
            else:
              print ('Problems parsing smarts: %s' %(sma))
  return order,patts


def _ReadPatts(fileName):
  """ *Internal Use Only*
    parses the pattern list from the data file
  """
  patts = {}
  order = []
  for line in open(fileName,'r').readlines():
    if line[0] != '#':
      splitLine = line.split('\t')
      if len(splitLine)>=4 and splitLine[0] != '':
        sma = splitLine[1]
        if sma!='SMARTS':
          sma.replace('"','')
          try:
            p = Chem.MolFromSmarts(sma)
          except:
            pass
          else:
            if p:
              if len(splitLine[0])>1 and splitLine[0][1] not in 'S0123456789':
                cha = splitLine[0][:2]
              else:
                cha = splitLine[0][0]
              logP = float(splitLine[2])
              if splitLine[3] != '':
                mr = float(splitLine[3])
              else:
                mr = 0.0
              if cha not in order:
                order.append(cha)
              l = patts.get(cha,[])
              l.append((sma,p,logP,mr))
              patts[cha] = l
            else:
              print ('Problems parsing smarts: %s' %(sma))
  return order,patts


def _pyGetAtomContribs(mol,patts=None,order=None, verbose=1):
  """ *Internal Use Only*
    calculates atomic contributions to the LogP and MR values
    if the argument *force* is not set, we'll use the molecules stored
    _crippenContribs value when possible instead of re-calculating.
  **Note:** Changes here affect the version numbers of MolLogP and MolMR
    as well as the VSA descriptors in Chem.MolSurf
  """

  nAtoms = mol.GetNumAtoms()
  atomContribs = [(0., 0.)]*nAtoms
  doneAtoms=[0]*nAtoms
  nAtomsFound=0
  done = False
  for cha in order:
    pattVect = patts[cha]
    for sma,patt,logp,mr in pattVect:
      #print 'try:',entry[0]
      for match in mol.GetSubstructMatches(patt,False,False):
        firstIdx = match[0]
        if not doneAtoms[firstIdx]:
          doneAtoms[firstIdx] = 1
          atomContribs[firstIdx] = (logp, mr)
          if verbose:
            print ('\tAtom %d: %s %4.4f %4.4f'%(match[0], sma, logp, mr))
          nAtomsFound+=1
          if nAtomsFound >= nAtoms:
            done = True
            break
    if done: break
  return atomContribs


def _pyGetAtomContribs_Ab(mol, patts=None, order=None, verbose=0):

      nAtoms = mol.GetNumAtoms()
      atomContribs = [(0., 0., 0., 0., 0.)] * nAtoms
      doneAtoms = [0] * nAtoms
      nAtomsFound = 0
      done = False
      for cha in order:
          pattVect = patts[cha]

          for sma, patt, r2, pi2h, b2h, b2o, l16 in pattVect:
              # print 'try:',entry[0]
              for match in mol.GetSubstructMatches(patt, False, False):
                  firstIdx = match[0]
                  if not doneAtoms[firstIdx]:
                      doneAtoms[firstIdx] = 1
                      atomContribs[firstIdx] = (r2, pi2h, b2h, b2o, l16)
                      if verbose:
                          print('\tAtom %d: %s %4.4f %4.4f %4.4f %4.4f %4.4f' % (match[0], sma, r2, pi2h, b2h, b2o, l16))
                      nAtomsFound += 1
                      if nAtomsFound >= nAtoms:
                          done = True
                          break
          if done: break
      return atomContribs


def eliminar_duplicados(b):
    my_list = list(b)
    para_borrar=[0]*len(b)
    for i in range(len(b) - 1):
        for j in range(i+1, len(b)):
            if sorted(b[i]) == sorted(b[j]):
                para_borrar[j]=1
    for i in range(len(para_borrar) - 1, 0, -1):
        if para_borrar[i] == 1:
            my_list.pop(i)
    b = tuple(my_list)
    return b


def eliminar_1_2_iguales(b):
    my_list = list(b)
    para_borrar=[0]*len(b)
    for i in range(len(b)-1):
        for j in range(i+1, len(b)):
            if (b[i][0] == b[j][0] or b[i][0] == b[j][len(b[j])-1]) and (b[i][len(b[i])-1] == b[j][0] or b[i][len(b[i])-1] == b[j][len(b[j])-1]):
                para_borrar[j]=1
    for i in range(len(para_borrar) - 1, 0, -1):
        if para_borrar[i] == 1:
            my_list.pop(i)
    b = tuple(my_list)
    return b


def _pyGetAtomContribs_Ab2(mol, patts = None, order = None, verbose=0):

    natoms = mol.GetNumAtoms()
    atomcontribs = [(0., 0., 0., 0., 0.)] * natoms
    doneatoms = [0] * natoms
    doneatoms2 = [0] * natoms
    for cha in order:
        pattvect = patts[cha]
        n=0
        for sma, patt, r2, pi2h, b2h, b2o, l16, tipo in pattvect:
            #print(n)
            n = n+1
            coincidencias = mol.GetSubstructMatches(patt, False, False)
            #print(coincidencias)
            if tipo == 3: # simetric moiety asignar la mitad al primero y al último y eliminar fragmentos duplicados
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    res = tuple(map(sum, zip(atomcontribs[match[0]], (r2 / 2, pi2h / 2, b2h / 2, b2o / 2, l16 / 2))))
                    res2 = tuple(map(sum, zip(atomcontribs[match[len(match)-1]], (r2 / 2, pi2h / 2, b2h / 2,
                                                                                      b2o / 2, l16 / 2))))
                    atomcontribs[match[0]] = res
                    atomcontribs[match[len(match)-1]] = res2
                    if verbose:
                        print('\tAtomos %d - %d: E%s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f' %
                                  (match[0], match[len(match)-1], n, sma, r2, pi2h, b2h, b2o, l16, tipo))
            elif tipo == 5: # simetric moiety asignar la mitad al primero y al último y eliminar fragmentos duplicados o con primero y ultimo coinciden
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                if len(coincidencias) > 2:
                    coincidencias = eliminar_1_2_iguales(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if (coincidencias[0][0] == coincidencias[1][0] or coincidencias[0][0] == coincidencias[1][len(coincidencias[1])-1]) and (coincidencias[0][len(coincidencias[0])-1] == coincidencias[1][0] or coincidencias[0][len(coincidencias[0])-1] == coincidencias[1][len(coincidencias[1])-1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    res = tuple(map(sum, zip(atomcontribs[match[0]], (r2 / 2, pi2h / 2, b2h / 2, b2o / 2, l16 / 2))))
                    res2 = tuple(map(sum, zip(atomcontribs[match[len(match)-1]], (r2 / 2, pi2h / 2, b2h / 2,
                                                                                      b2o / 2, l16 / 2))))
                    atomcontribs[match[0]] = res
                    atomcontribs[match[len(match)-1]] = res2
                    if verbose:
                        print('\tAtomos %d - %d: E%s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f' %
                                  (match[0], match[len(match)-1], n, sma, r2, pi2h, b2h, b2o, l16, tipo))
            elif tipo==2:
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    for idx in range(len(match)):
                        res = tuple(map(sum, zip(atomcontribs[match[idx]], (
                        r2 / len(match), pi2h / len(match), b2h / len(match), b2o / len(match), l16 / len(match)))))
                        atomcontribs[match[idx]] = res
                        if verbose:
                            print('\tAtom %d: E%s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f' % (match[idx], n, sma, r2, pi2h, b2h, b2o, l16, tipo))
            else:
                for match in coincidencias:
                    if tipo == 4:#contr solo primero y una sola vez
                        firstIdx = match[0]
                        #print("\tAtom %d: %4.4f" % (match[0], tipo))
                        if not doneatoms2[firstIdx]:
                            for idx in range(len(match)):
                                doneatoms2[match[idx]] = 1
                            #print(doneatoms2)
                            res = tuple(map(sum, zip(atomcontribs[firstIdx], (r2 , pi2h , b2h, b2o, l16))))
                            atomcontribs[firstIdx] = res
                            if verbose:
                                print('\tAtom %d: E%s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f' % (match[0], n, sma, r2, pi2h, b2h, b2o, l16,tipo))
                    elif tipo == 1: #contr solo primero sin repetición
                        firstIdx = match[0]
                        #print("\tAtom %d: %4.4f" % (match[0], tipo))
                        if not doneatoms[firstIdx]:
                            doneatoms[firstIdx] = 1
                            res = tuple(map(sum, zip(atomcontribs[firstIdx], (r2, pi2h, b2h, b2o, l16))))
                            atomcontribs[firstIdx] = res
                            if verbose:
                                print('\tAtom %d: E%s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f' % (match[0], n, sma, r2, pi2h, b2h, b2o, l16,tipo))

    return atomcontribs


def _pyGetAtomContribs_Abalpha(mol, patts=None, order=None, verbose=0):

    natoms = mol.GetNumAtoms()
    atomcontribs = [0] * natoms
    doneatoms = [0] * natoms
    doneatoms2= [0] * natoms
    for cha in order:
        pattVect = patts[cha]
        n=0
        for sma, patt, a2, tipo in pattVect:
            # print 'try:',entry[0]
            #print(sma)
            n = n + 1
            coincidencias = mol.GetSubstructMatches(patt, False, False)
            #print(coincidencias)
            if tipo == 3:  # simetric moiety asignar la mitad al primero y al último y eliminar fragmentos duplicados
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    atomcontribs[match[0]] += a2/2
                    atomcontribs[match[len(match) - 1]] += a2/2
                    if verbose:
                        print('\tAtomos %d - %d: E%s %s %4.4f %4.4f' %(match[0], match[len(match) - 1], n, sma, a2, tipo))
            elif tipo == 6:  # simetric moiety asignar la mitad al primero y al último y eliminar fragmentos duplicados
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                if len(coincidencias) > 2:
                    coincidencias = eliminar_1_2_iguales(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if (coincidencias[0][0] == coincidencias[1][0] or coincidencias[0][0] == coincidencias[1][
                        len(coincidencias[1]) - 1]) and (
                            coincidencias[0][len(coincidencias[0]) - 1] == coincidencias[1][0] or coincidencias[0][
                        len(coincidencias[0]) - 1] == coincidencias[1][len(coincidencias[1]) - 1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    atomcontribs[match[0]] += a2/2
                    atomcontribs[match[len(match) - 1]] += a2/2
                    if verbose:
                        print('\tAtomos %d - %d: E%s %s %4.4f %4.4f' %(match[0], match[len(match) - 1], n, sma, a2, tipo))
            elif tipo==2:# simetrico todos
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    for idx in range(len(match)):
                        atomcontribs[match[idx]] += a2/len(match)
                        if verbose:
                            print('\tAtom %d: E%s %s %4.4f %4.4f' % (match[idx], n, sma, a2, tipo))
            elif tipo == 5: #solo primero
                if len(coincidencias) > 2:
                    coincidencias = eliminar_duplicados(coincidencias)
                elif len(coincidencias) == 2:
                    my_list = list(coincidencias)
                    if sorted(coincidencias[0]) == sorted(coincidencias[1]):
                        my_list.pop(1)
                    coincidencias = tuple(my_list)
                #print(coincidencias)
                for match in coincidencias:
                    atomcontribs[match[0]] += a2
                    if verbose:
                        print('\tAtom %d: E%s %s %4.4f %4.4f' % (match[0], n, sma, a2, tipo))
            else:
                for match in coincidencias:
                    if tipo == 4:  # contr solo primero y una sola vez
                        firstIdx = match[0]
                        #print("\tAtom %d: %4.4f" % (match[0], tipo))
                        if not doneatoms2[firstIdx]:
                            for idx in range(len(match)):
                                doneatoms2[match[idx]] = 1
                            #print(doneatoms2)
                            atomcontribs[firstIdx] += a2
                            if verbose:
                                print('\tAtom %d: E%s %s %4.4f %4.4f' % (match[0], n, sma, a2, tipo))
                    elif tipo == 1:  # contr solo primero sin repetición
                        firstIdx = match[0]
                        #print("\tAtom %d: %4.4f" % (match[0], tipo))
                        if not doneatoms[firstIdx]:
                            doneatoms[firstIdx] = 1
                            atomcontribs[firstIdx] += a2
                            if verbose:
                                print('\tAtom %d: E%s %s %4.4f %4.4f' % (match[0], n, sma, a2, tipo))
    #print(atomcontribs)
    return atomcontribs


def _pyGetContribs(mol, patts, order):
    nbonds = mol.GetNumBonds()
    bondContribs = [0.] * nbonds
    doneBonds = [0] * nbonds
    nBondsFound = 0
    done = False
    for cha in order:
        pattVect = patts[cha]
        for sma, patt, valor in pattVect:
            # print 'try:',entry[0]
            if len(mol.GetSubstructMatches(patt, False, False))>0:
                for match in mol.GetSubstructMatches(patt, False, False):
                    aid1 = match[0]
                    aid2 = match[1]
                    firstIdx = mol.GetBondBetweenAtoms(aid1, aid2).GetIdx()
                    #for bond in patt.GetBonds():
                    #   aid1 = match[bond.GetBeginAtomIdx()]
                    #   aid2 = match[bond.GetEndAtomIdx()]
                    #    firstIdx = mol.GetBondBetweenAtoms(aid1, aid2).GetIdx()
                    if not doneBonds[firstIdx]:
                        doneBonds[firstIdx] = 1
                        bondContribs[firstIdx] = valor
                        #print('\tEnlace %d: %s %4.4f' % (firstIdx, sma, valor))
                        nBondsFound += 1
                        if nBondsFound >= nbonds:
                            done = True
                            break
        if done: break
    return bondContribs


def _pyGetAtomContribs_Pol(mol, patts=None, order=None, verbose=0):

    nAtoms = mol.GetNumAtoms()
    atomContribs = [0.] * nAtoms
    doneAtoms = [0] * nAtoms
    nAtomsFound = 0
    done = False
    for cha in order:
        pattVect = patts[cha]
        for sma, patt, pol in pattVect:
            # print 'try:',entry[0]
            for match in mol.GetSubstructMatches(patt, False, False):
                firstIdx = match[0]
                #print(match)
                if not doneAtoms[firstIdx]:
                    doneAtoms[firstIdx] = 1
                    atomContribs[firstIdx] = (pol)
                    if verbose:
                        print('\tAtom %d: %s %4.4f ' % (match[0], sma, pol))
                    nAtomsFound += 1
                    if nAtomsFound >= nAtoms:
                        done = True
                        break
        if done: break
    return atomContribs


def suma(matriz):
    salida = 0
    for i in range(len(matriz)):
        salida = salida+matriz[i, i]
    return salida


def calcular_dipolos2(mol):
    enlaces = mol.GetBonds()
    contrib = [0.]*mol.GetNumBonds()
    for enlace in enlaces:
        atomo1 = enlace.GetBeginAtom()
        atomo2 = enlace.GetEndAtom()
        if atomo1.GetSymbol() == 'F' or atomo2.GetSymbol() == 'F':
            contrib[enlace.GetIdx()] = 1.39
        if atomo1.GetSymbol() == 'Cl' or atomo2.GetSymbol() == 'Cl':
            contrib[enlace.GetIdx()] = 1.47
        if atomo1.GetSymbol() == 'Br' or atomo2.GetSymbol() == 'Br':
            contrib[enlace.GetIdx()] = 1.42
        if atomo1.GetSymbol() == 'C' or atomo2.GetSymbol() == 'C':
            if atomo1.GetSymbol() == 'O' or atomo2.GetSymbol() == 'O':
                if enlace.GetBondTypeAsDouble() == 2 or enlace.GetBondTypeAsDouble() == 1.5:
                    contrib[enlace.GetIdx()] = 2.40
                else:
                    contrib[enlace.GetIdx()] = 0.70
        if atomo1.GetSymbol() == 'C' or atomo2.GetSymbol() == 'C':
            if atomo1.GetSymbol() == 'N' or atomo2.GetSymbol() == 'N':
                if enlace.GetBondTypeAsDouble() == 2 or enlace.GetBondTypeAsDouble() == 1.5:
                    contrib[enlace.GetIdx()] = 1.40
                elif enlace.GetBondTypeAsDouble() == 3:
                    contrib[enlace.GetIdx()] = 3.10
                else:
                    contrib[enlace.GetIdx()] = 0.45
        if atomo1.GetSymbol() == 'C' or atomo2.GetSymbol() == 'C':
            if atomo1.GetSymbol() == 'S' or atomo2.GetSymbol() == 'S':
                if enlace.GetBondTypeAsDouble() == 2:
                    contrib[enlace.GetIdx()] = 2.00
                else:
                    contrib[enlace.GetIdx()] = 0.80
    return contrib

def extraer_valores(mol, property):
    num = mol.GetNumBonds()
    c = 0
    p = ""
    for enlace in mol.GetBonds():
        if c == num:
            p = p + enlace.GetProp(property)
        else:
            p = p + enlace.GetProp(property) + ", "
        c = c + 1
    return p

def calcular(mol, prop, n, tipo):
    c = 0
    if tipo == 'bond':
        ea = adj_enlace(mol)
        des = np.empty([n])
        for enlace in mol.GetBonds():
            ea[c, c] = round(float(enlace.GetProp(prop + str(1))), 6)
            c = c + 1
        des[0] = np.sum(ea.diagonal())
        multipl = ea
        for j in range(1, n):
            multipl = np.dot(multipl, ea)
            des[j] = np.sum(multipl.diagonal())
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
            des[j] = np.sum(multipl.diagonal())

    return des

def adj_enlace(mol):
    c = mol.GetNumBonds()
    matriz = np.zeros((c, c))
    enlaces = mol.GetBonds()
    for i in range(c):
        for j in range(i+1, c):
            if enlaces[i].GetBeginAtomIdx() == enlaces[j].GetBeginAtomIdx() or enlaces[i].GetBeginAtomIdx() == enlaces[
                j].GetEndAtomIdx() or enlaces[i].GetEndAtomIdx() == enlaces[j].GetBeginAtomIdx() or enlaces[
                i].GetEndAtomIdx() == enlaces[j].GetEndAtomIdx():
                matriz[i, j] = 1
                matriz[j, i] = 1
    return matriz
