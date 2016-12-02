#!/usr/bin/env python

import simuOpt
simuOpt.setOptions(optimized=False, debug='DBG_WARNING', gui=False, alleleType='binary')

import simuPOP as sim
from simuPOP.utils import Exporter
from simuPOP.demography import AdmixtureEvent
import random

def mozParentsChooser(pop):
    '''Choose males from wild pop, females from Au pop'''
    males = [x for x in pop.individuals(subPop=[1]) if x.sex() == 1]
    print(len(males))
    females = [x for x in pop.individuals(subPop=[0]) if x.sex() == 2]
    print(len(females))
    while True:
        male = males[random.randint(0, len(males) - 1)]
        female = females[random.randint(0, len(females) - 1)]
        yield(male, female)

def printFreq(pop, loci):
    sim.stat(pop, alleleFreq=loci)
    print(', '.join(['{:.3f}'.format(pop.dvars().alleleFreq[x][0]) for x in loci]))

pop = sim.Population([100, 10000],
                     loci=[2, 2],
                     lociPos=[
                         0.1, 10.5,
                         3.5, 23.5],
                     ancGen=5,
                     alleleNames=['1', '2'],
                     subPopNames=['colony', 'wild'],
                     infoFields=['ind_id', 'father_id', 'mother_id']
                     )
sim.initSex(pop)

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
        sim.InitGenotype(subPops=[0], loci=[0], freq=[0,1]),
        sim.InitGenotype(subPops=[0], loci=[1], freq=[0, 1]),
        sim.InitGenotype(subPops=[0], loci=[2], freq=[0, 1]),
        sim.InitGenotype(subPops=[0], loci=[3], freq=[0, 1]),
        sim.InitGenotype(subPops=[1], loci=[0], freq=[1, 0]),
        sim.InitGenotype(subPops=[1], loci=[1], freq=[1, 0]),
        sim.InitGenotype(subPops=[1], loci=[2], freq=[1, 0]),
        sim.InitGenotype(subPops=[1], loci=[3], freq=[1, 0]),
    ],
    #preOps=sim.Migrator(rate=[100], mode=sim.BY_COUNTS, subPops=[0], toSubPops=[1]),
    matingScheme=sim.HomoMating(
        sim.PyParentsChooser(mozParentsChooser),
        sim.OffspringGenerator(ops=[
            sim.Recombinator(intensity=0.01),
            sim.IdTagger(),
            sim.PedigreeTagger()
        ]),
        subPopSize=sim.demography.AdmixtureEvent()
    ),
    gen=5
)
printFreq(pop, [0, 1])

sim.utils.export(pop, format='PED', output='test.ped')

#sim.initGenotype(pop, freq=[.4, .6])
#sim.dump(pop, max=6, structure=False)
#printFreq(pop, range(5))


#
# pop.evolve(
#     initOps = [
#         sim.InitSex(),
#         sim.InitGenotype()
#     ]
# )
#
# pop = sim.Population(size=[10000]*2, loci=10, lociPos=range(5) + range(10, 15))
# pop.evolve(
#     initOps=[
#         sim.InitSex(),
#         sim.InitGenotype(haplotypes=[[0]*10, [1]*10]),
#     ],
#     matingScheme=sim.RandomMating(ops=sim.Recombinator(intensity=0.0005)),
#     postOps=[
#         sim.Stat(LD=[[1,2], [4,5], [8,9], [0,9]], step=10),
#         sim.PyEval(r"'gen=%d\tLD12=%.3f (%.3f)\tLD45=%.3f (%.3f)\tLD09=%.3f\n'%"
#                    "(gen, LD[1][2], 0.25*0.9995**(gen+1), LD[4][5],"
#                    "0.25*0.9975**(gen+1),LD[0][9])", step=10)
#     ],
#     gen=100
# )
