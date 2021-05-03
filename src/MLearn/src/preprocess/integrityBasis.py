import numpy as np
import pdb

# This module compute integrity bases for the following set
# {S, \Omega, grad(P), grad(K)}
# S = rate of strain tensor,
# \Omega=rotation tensor
# grad(P) = pressur gradient
# grad(K) = gradient of turbulent kinetic energy


def field_to_bases(S_all, W_all, gradP_all, gradK_all):
    """ 
    Compute all basis given raw tensors and vectors

    input:
        S_all: Strain rate tensor array (nCell by 3 by 3)
        W_all: Rotation rate tensor array (nCell by 3 by 3)
        gradP_all: pressure gradient array (nCell by 3 by 1)
        gradK_all: TKE gradient array (nCell by 3 by 1)
    output:
        featureMatrix: nCell feature vectors (nCell by 47)
    """

    nCell = S_all.shape[0]
    nFeature = 47
    featureMatrix = np.zeros([nCell, nFeature])
    # assert shapes are correct
    assert S_all.shape == (nCell, 3, 3), "Please give array of (nCell by 3 by 3)"
    assert W_all.shape == (nCell, 3, 3), "Please give array of (nCell by 3 by 3)"
    assert gradP_all.shape == (nCell, 3, 1), "Please give array of (nCell by 3 by 1)"
    assert gradK_all.shape == (nCell, 3, 1), "Please give array of (nCell by 3 by 1)" 

    blist = []
    
    iCell = 0
    featureDescription = ['strainRate', 'rotationRate', 'gradP', 'gradK']
    for S, W, gradP, gradK  in zip(S_all, W_all, gradP_all, gradK_all):               
        Ap = vector_to_antisymm(gradP)
        Ak = vector_to_antisymm(gradK)
        bb, desSL = computer_tensor_bases(S, W, Ap, Ak, featureDescription)
        #for i in range(47):
            #print i+1, desSL[i]
        featureMatrix[iCell, :] = bb
        iCell = iCell + 1
    return featureMatrix, desSL


def vector_to_antisymm(v):
    """
    map a give vector $v$ to an antisymmetric tensor $A$
    A = - I $\times$ v, where $\times$ indicate cross product
    input:
        vector v
    output:
        antisymmetric tensor A
    """

    A = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    A = -A;
    
    return A

def computer_tensor_bases(S, A1, A2, A3, des):
    """
    compute bases for set {S, A_1, A_2, A_3}
    following the Table A1 in Johnson, Handbook of Fluid Dynamics
    input: S, A1, A2, A3
        S = symmetric tensor
        A_1, A_2, A_3 = antisymmetric tensor
        deslist = a list of words description the characteristics of the S, A1, A2, A3
    output:
        an array of 45 bases, with the order specified below
        an arry of set with description of all the features
    """
    
    # It is easier to append to a list by + operation
    # will convert to array later
    blist = []
    desSetList = []
    
    # ns, na = (1, 0);   b[0:1];  Total terms: 2
    # S^2, S^3
    # Note: trace(S) = 0 for incompressible flows
    Sp2 = np.dot(S, S)
    Sp3 = np.dot(Sp2, S)
    blist += [Sp2.trace(),  Sp3.trace()]
    desSetList += [set([des[0]])] * 2
    
    #pdb.set_trace()
    # ns, na = (0, 1);  b[2:4];  Total terms: 3
    # A1^2, A2^2, A3^2
    for Ai, desI in zip([A1, A2, A3], des[1:]):
        Aip2 = np.dot(Ai, Ai)
        blist += [Aip2.trace()]
        desSetList += [ set([desI]) ]

    # ns, na = (1, 1); b[5:13];  Total terms: 9
    # Ai^2 * S, Ai^2 * S^2, Ai^2 * S * Ai * S^2 
    for Ai, desI in zip([A1, A2, A3], des[1:]):
        #print Ai, desI
        Aip2_S = Ai.dot(Ai).dot(S)
        Aip2_Sp2 = Aip2_S.dot(S)
        Aip2_S_Ai_Sp2 = Aip2_S.dot(Ai).dot(Sp2)
        blist += [Aip2_S.trace(), Aip2_Sp2.trace(), Aip2_S_Ai_Sp2.trace()]
        desSetList += [ set( [desI, des[0]]) ] * 3

    # ns, na = (0, 2); b[14:16];  Total terms: 3
    # Ai * Aj
    for Ai, Aj, desI in [(A1, A2, des[1:3]), (A2, A3, des[2:]), (A1, A3, [des[ii] for ii in [1, 3]])]:
        blist += [np.dot(Ai, Aj).trace()]
        desSetList += [set(desI)]
    
    # ns, na = (1, 2); b[17:40]; Total terms: 24
    # & indicate cyclic permulation on A index,
    # e.g., Ai^2 * Aj (&) also include Aj^2 * Aj
    # Compute the following
    # Ai * Aj * S
    # Ai * Aj * S^2
    # Ai^2 * Aj * S (&)
    # Ai^2 * Aj * S^2 (&)
    # Ai^2 * S * Aj * S^2 (&)
    for Am, An, desI in [(A1, A2, des[1:3]), (A1, A3, [des[ii] for ii in [1, 3]]), (A2, A3, des[2:])]:
        # three different subsets
        Ai_Aj_S =  Am.dot(An).dot(S)
        Ai_Aj_Sp2 = Ai_Aj_S.dot(S)
        blist += [Ai_Aj_S.trace(), Ai_Aj_Sp2.trace()]
        desSetList += [ set(desI + [des[0]]) ] * 2
        for Ai, Aj in [(Am, An), (An, Am)]:
        # two different "cyclic permutations"
            Aip2_Aj_S = Ai.dot(Ai).dot(Aj).dot(S) 
            Aip2_Aj_Sp2 = Ai.dot(Ai).dot(Aj).dot(S).dot(S)
            Aip2_S_Aj_Sp2 = Ai.dot(Ai).dot(S).dot(Aj).dot(S).dot(S)
            blist += [Aip2_Aj_S.trace(), Aip2_Aj_Sp2.trace(),  Aip2_S_Aj_Sp2.trace()]

            desSetList += [ set(desI + [des[0]]) ] * 3

    # ns, na = (0, 3); b[41]
    # A1* A2 * A3
    blist += [A1.dot(A2).dot(A3).trace()]
    desSetList += [ set(des[1:]) ]

    # ns, na = (1, 3); b[42:46]; Total terms: 5
    # A1 * A2 * A3 * S
    # A1 * A3 * A2 * S
    # A1 * A2 * A3 * Sp2
    # A1 * A3 * A2 * Sp2
    # A1 * A2 * S * A3 * Sp2
    A1_A2_A3_S = A1.dot(A2).dot(A3).dot(S)
    A1_A3_A2_S = A1.dot(A3).dot(A2).dot(S)
    A1_A2_A3_Sp2 = A1_A2_A3_S.dot(S)
    A1_A3_A2_Sp2 = A1_A3_A2_S.dot(S)
    A1_A2_S_A3_Sp2 = A1.dot(A2).dot(S).dot(A3).dot(S).dot(S)

    blist += [
        A1_A2_A3_S.trace(), 
        A1_A3_A2_S.trace(),
        A1_A2_A3_Sp2.trace(),
        A1_A3_A2_Sp2.trace(),
        A1_A2_S_A3_Sp2.trace()
        ]
    # all the 5 entries above involves all features
    desSetList += [ set(des) ] * 5

    bb = np.array(blist)

    return bb, desSetList
