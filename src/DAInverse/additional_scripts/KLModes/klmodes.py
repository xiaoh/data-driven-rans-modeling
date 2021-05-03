
import numpy as np
import os

outDir = 'pehill_base-observation/0'
varList = ['K', 'Eta', 'Xi'] #varListAll = ['K', 'Eta', 'Xi', 'VA', 'VB', 'VC']

def get_header(name,N):
    header = \
        '/*--------------------------------*- C++ -*----------------------------------*\\\n'  + \
        '| =========                 |                                                 |  \n' + \
        '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n' + \
        '|  \\\\    /   O peration     | Version:  2.4.x                                 |\n' + \
        '|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n' + \
        '|    \\\\/     M anipulation  |                                                 |\n' + \
        '\\*---------------------------------------------------------------------------*/ \n' + \
        'FoamFile\n' + \
        '{\n' + \
        '    version     2.0;\n' + \
        '    format      ascii;\n' + \
        '    class       volScalarField;\n' + \
        '    location    "0";\n' + \
        '    object      ' + name + ';\n' + \
        '}\n' + \
        '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n' + \
        '\n' + \
        'dimensions      [0 2 -2 0 0 0 0];\n' + \
        '\n' + \
        'internalField   nonuniform List<scalar>\n' + \
        str(N) + '\n' + \
        '(\n'
    return header

footer = ')' +'\n' + \
    ';' + '\n' + \
    '\n' + \
    'boundaryField' + '\n' + \
    '    {' + '\n' + \
    '    topWall' + '\n' + \
    '    {' + '\n' + \
    '        type            zeroGradient;' + '\n' + \
    '    }' + '\n' + \
    '    bottomWall' + '\n' + \
    '    {' + '\n' + \
    '        type            zeroGradient;' + '\n' + \
    '    }' + '\n' + \
    '    inlet' + '\n' + \
    '    {' + '\n' + \
    '        type            zeroGradient;' + '\n' + \
    '    }' + '\n' + \
    '    outlet' + '\n' + \
    '    {' + '\n' + \
    '        type            fixedValue;' + '\n' + \
    '        value           uniform 0;' + '\n' + \
    '    }' + '\n' + \
    '    defaultFaces' + '\n' + \
    '    {' + '\n' + \
    '        type            empty;' + '\n' + \
    '    }' + '\n' + \
    '}' + '\n' + \
    '\n' + \
    '// ************************************************************************* //'

for var in varList:
    file = './randomData_' + var + '/KLModes.dat'
    KLmodes = np.loadtxt(file)
    Ncells,Nmodes = KLmodes.shape
    for i in range(Nmodes):
        name = 'KLMode_' + var + '_' + str(i+1)
        header = get_header(name,Ncells)
        f = open(outDir+os.sep+name, 'w')
        f.write(header)
        for j in range(Ncells):
            f.write(str(KLmodes[j,i]) + '\n')
        f.write(footer)
        f.close()
