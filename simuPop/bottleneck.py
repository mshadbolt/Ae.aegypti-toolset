#!/usr/bin/env python

import simuOpt

# simuOpt.setOptions(numThreads=x)

import simuPOP as sim
from simuPOP.plotter import VarPlotter

def demo(gen):
    if gen >= 4 and gen <5:
        return 20
    else:
        return 1000

pop = sim.Population(size=1000, loci=1)
simu = sim.Simulator(pop, rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    postOps=[
        sim.Stat(alleleFreq=0),
        VarPlotter('alleleFreq[0][0]', update=1, saveAs='Figures/bottleneck.pdf')
    ],
    gen=11
)
