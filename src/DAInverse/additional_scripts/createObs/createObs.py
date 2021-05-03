
import numpy as np

obs = np.loadtxt('obsdata.txt')
Ubulk = 0.028

obs[:,3:] *= Ubulk
print('Number of observations: ' + str(obs.shape[0]))


header = \
'/*--------------------------------*- C++ -*----------------------------------*\\\n'  + \
'| =========                 |                                                 |  \n' + \
'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n' + \
'|  \\\\    /   O peration     | Version:  2.1.x                                 |\n' + \
'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n' + \
'|    \\\\/     M anipulation  |                                                 |\n' + \
'\\*---------------------------------------------------------------------------*/ \n' + \
'FoamFile\n' + \
'{\n' + \
'    version     2.0;\n' + \
'    format      ascii;\n' + \
'    class       dictionary;\n' + \
'    location    "constant";\n' + \
'    object      obsLocations;\n' + \
'}\n' + \
'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n' + \
'\n' + \
'(\n' 

footer = ');'


#f = open('./constant/obsLocations', 'w')
f = open('./obsLocations', 'w')
f.write(header)
for i in range(len(obs)):
    f.write('(' + str(obs[i,0]) + ' ' + str(obs[i,1]) + ' ' + str(obs[i,2]) + ')\n')
f.write(footer)
f.close()

#f = open('observationData/obsVelocity','w')
f = open('obsVelocity','w')
for i in range(len(obs)):
    for j in [3,4,5]:
        f.write(str(obs[i,j])+'\n')
f.close()

#f = open('observationData/obsVPlot','w')
f = open('obsVPlot','w')
for i in range(len(obs)):
    f.write(str(obs[i,0]) + ' ' + str(obs[i,1]) + ' ' + str(obs[i,2]) + ' ' + str(obs[i,3]) + ' ' + str(obs[i,4]) + ' ' + str(obs[i,5]) + '\n')
f.close()
