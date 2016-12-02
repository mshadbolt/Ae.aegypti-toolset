#!/usr/bin/env python

import simuOpt

# simuOpt.setOptions(numThreads=x)

import simuPOP as sim

pop = sim.Population(size=[100]*100, loci=1)

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.PyOutput('gen:  mean freq   mean Ht (expected Ht)\n')
    ],
    preOps=[
        sim.Stat(alleleFreq=0, heteroFreq=0,
                 vars=['alleleFreq_sp', 'heteroFreq_sp'], step=20),
        sim.PyEval(r'"%2d:    %.4f        %.4f (%.4f)\n" % (gen,'
                   'sum([subPop[x]["alleleFreq"][0][1] for x in range(100)])/100.,'
                   'sum([subPop[x]["heteroFreq"][0] for x in range(100)])/100.,'
                   '0.5*(1-1/200.)**gen)', step=10)
    ],
    matingScheme=sim.RandomMating(),
    gen=100
)

