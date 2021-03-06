import simuPOP as sim
from simuPOP.utils import Exporter

Inari=sim.Population(
   size=[1000,1000],
   ploidy=2,
   loci=[17,2,11,5,7,5,3,3,14,15,5,8,16,13,6,6,2,8,4,6,6,5,3,9,5,8,4,4],
   lociPos=[
      4.9,21.4,23.9,26.0,34.1,34.4,48.9,50.8,78.6,82.4,118.2,119.9,120.1,122.3,131.1,136.1,151.3,
      24.7,43.3,
      28.1,28.6,32.1,37.9,40.3,59.6,63.6,71.9,74.8,77.4,91.8,
      9.0,20.0,47.5,47.7,55.8,
      22.3,22.7,23.7,34.7,34.8,39.9,54.5,
      35.9,38.9,55.9,59.8,68.0,
      5.5,18.0,47.0,
      18.5,19.8,22.0,
      18.2,18.3,31.5,47.8,70.2,84.4,94.4,96.3,103.9,110.5,111.5,120.1,128.7,139.0,
      9.3,17.8,19.3,20.8,24.2,29.1,31.7,37.0,64.1,67.9,74.6,74.9,75.0,75.4,113.0,
      5.7,14.5,19.4,23.2,57.6,
      26.4,34.2,34.3,36.9,61.5,62.5,66.6,78.4,
      9.3,12.5,20.1,22.2,23.2,27.6,44.8,53.1,53.9,55.8,71.0,75.0,75.1,75.2,76.5,85.5,
      14.8,15.2,23.9,34.1,37.0,54.3,54.6,55.4,60.7,61.1,64.8,80.1,82.3,
      5.5,40.9,41.9,48.7,52.8,60.0,
      16.9,42.0,68.1,69.0,71.0,81.9,
      27.0,50.8,
      13.9,23.3,24.2,27.1,28.2,53.6,56.3,66.4,
      11.4,51.2,51.4,72.5,
      33.0,33.8,47.4,50.5,63.1,69.1,
      14.4,20.0,20.1,38.8,43.6,51.6,
      6.5,7.5,18.4,20.1,44.2,
      30.8,33.6,40.2,
      1.3,1.4,11.1,12.1,12.4,17.0,19.8,34.4,47.0,
      11.9,14.9,34.1,34.5,47.1,
      13.5,19.8,21.0,23.9,25.8,26.4,27.0,28.9,
      16.1,21.4,22.6,26.2,
      8.6,13.2,28.1,35.8
   ],

   ancGen=3,
   alleleNames=['1','2'],
   subPopNames=['escapee','wild'],
   infoFields=['migrate_to','ind_id','father_id', 'mother_id']
)
sim.initSex(Inari)

Inari.evolve(
   initOps=[
      
      sim.InitSex(),
      sim.IdTagger(),
      sim.InitGenotype(subPops=[0] , loci=[0], freq=[0.515,0.485]),
      sim.InitGenotype(subPops=[0] , loci=[1], freq=[0.617,0.383]),
      sim.InitGenotype(subPops=[0] , loci=[2], freq=[0.37,0.63]),
      sim.InitGenotype(subPops=[0] , loci=[3], freq=[0.447,0.553]),
      sim.InitGenotype(subPops=[0] , loci=[4], freq=[0.485,0.515]),
      sim.InitGenotype(subPops=[0] , loci=[5], freq=[0.599,0.401]),
      sim.InitGenotype(subPops=[0] , loci=[6], freq=[0.432,0.568]),
      sim.InitGenotype(subPops=[0] , loci=[7], freq=[0.495,0.505]),
      sim.InitGenotype(subPops=[0] , loci=[8], freq=[0.385,0.615]),
      sim.InitGenotype(subPops=[0] , loci=[9], freq=[0.576,0.424]),
      sim.InitGenotype(subPops=[0] , loci=[10], freq=[0.528,0.472]),
      sim.InitGenotype(subPops=[0] , loci=[11], freq=[0.152,0.848]),
      sim.InitGenotype(subPops=[0] , loci=[12], freq=[0.443,0.557]),
      sim.InitGenotype(subPops=[0] , loci=[13], freq=[0.532,0.468]),
      sim.InitGenotype(subPops=[0] , loci=[14], freq=[0.548,0.452]),
      sim.InitGenotype(subPops=[0] , loci=[15], freq=[0.544,0.456]),
      sim.InitGenotype(subPops=[0] , loci=[16], freq=[0.545,0.455]),
      sim.InitGenotype(subPops=[0] , loci=[17], freq=[0.564,0.436]),
      sim.InitGenotype(subPops=[0] , loci=[18], freq=[0.367,0.633]),
      sim.InitGenotype(subPops=[0] , loci=[19], freq=[0.273,0.727]),
      sim.InitGenotype(subPops=[0] , loci=[20], freq=[0.473,0.527]),
      sim.InitGenotype(subPops=[0] , loci=[21], freq=[0.535,0.465]),
      sim.InitGenotype(subPops=[0] , loci=[22], freq=[0.597,0.403]),
      sim.InitGenotype(subPops=[0] , loci=[23], freq=[0.552,0.448]),
      sim.InitGenotype(subPops=[0] , loci=[24], freq=[0.516,0.484]),
      sim.InitGenotype(subPops=[0] , loci=[25], freq=[0.544,0.456]),
      sim.InitGenotype(subPops=[0] , loci=[26], freq=[0.583,0.417]),
      sim.InitGenotype(subPops=[0] , loci=[27], freq=[0.495,0.505]),
      sim.InitGenotype(subPops=[0] , loci=[28], freq=[0.473,0.527]),
      sim.InitGenotype(subPops=[0] , loci=[29], freq=[0.286,0.714]),
      sim.InitGenotype(subPops=[0] , loci=[30], freq=[0.552,0.448]),
      sim.InitGenotype(subPops=[0] , loci=[31], freq=[0.667,0.333]),
      sim.InitGenotype(subPops=[0] , loci=[32], freq=[0.333,0.667]),
      sim.InitGenotype(subPops=[0] , loci=[33], freq=[0.466,0.534]),
      sim.InitGenotype(subPops=[0] , loci=[34], freq=[0.492,0.508]),
      sim.InitGenotype(subPops=[0] , loci=[35], freq=[0.496,0.504]),
      sim.InitGenotype(subPops=[0] , loci=[36], freq=[0.604,0.396]),
      sim.InitGenotype(subPops=[0] , loci=[37], freq=[0.325,0.675]),
      sim.InitGenotype(subPops=[0] , loci=[38], freq=[0.465,0.535]),
      sim.InitGenotype(subPops=[0] , loci=[39], freq=[0.459,0.541]),
      sim.InitGenotype(subPops=[0] , loci=[40], freq=[0.511,0.489]),
      sim.InitGenotype(subPops=[0] , loci=[41], freq=[0.605,0.395]),
      sim.InitGenotype(subPops=[0] , loci=[42], freq=[0.478,0.522]),
      sim.InitGenotype(subPops=[0] , loci=[43], freq=[0.445,0.555]),
      sim.InitGenotype(subPops=[0] , loci=[44], freq=[0.53,0.47]),
      sim.InitGenotype(subPops=[0] , loci=[45], freq=[0.522,0.478]),
      sim.InitGenotype(subPops=[0] , loci=[46], freq=[0.721,0.279]),
      sim.InitGenotype(subPops=[0] , loci=[47], freq=[0.59,0.41]),
      sim.InitGenotype(subPops=[0] , loci=[48], freq=[0.528,0.472]),
      sim.InitGenotype(subPops=[0] , loci=[49], freq=[0.407,0.593]),
      sim.InitGenotype(subPops=[0] , loci=[50], freq=[0.554,0.446]),
      sim.InitGenotype(subPops=[0] , loci=[51], freq=[0.524,0.476]),
      sim.InitGenotype(subPops=[0] , loci=[52], freq=[0.056,0.944]),
      sim.InitGenotype(subPops=[0] , loci=[53], freq=[0.435,0.565]),
      sim.InitGenotype(subPops=[0] , loci=[54], freq=[0.46,0.54]),
      sim.InitGenotype(subPops=[0] , loci=[55], freq=[0.587,0.413]),
      sim.InitGenotype(subPops=[0] , loci=[56], freq=[0.331,0.669]),
      sim.InitGenotype(subPops=[0] , loci=[57], freq=[0.47,0.53]),
      sim.InitGenotype(subPops=[0] , loci=[58], freq=[0.366,0.634]),
      sim.InitGenotype(subPops=[0] , loci=[59], freq=[0.367,0.633]),
      sim.InitGenotype(subPops=[0] , loci=[60], freq=[0.506,0.494]),
      sim.InitGenotype(subPops=[0] , loci=[61], freq=[0.678,0.322]),
      sim.InitGenotype(subPops=[0] , loci=[62], freq=[0.489,0.511]),
      sim.InitGenotype(subPops=[0] , loci=[63], freq=[0.524,0.476]),
      sim.InitGenotype(subPops=[0] , loci=[64], freq=[0.66,0.34]),
      sim.InitGenotype(subPops=[0] , loci=[65], freq=[0.735,0.265]),
      sim.InitGenotype(subPops=[0] , loci=[66], freq=[0.59,0.41]),
      sim.InitGenotype(subPops=[0] , loci=[67], freq=[0.47,0.53]),
      sim.InitGenotype(subPops=[0] , loci=[68], freq=[0.442,0.558]),
      sim.InitGenotype(subPops=[0] , loci=[69], freq=[0.461,0.539]),
      sim.InitGenotype(subPops=[0] , loci=[70], freq=[0.486,0.514]),
      sim.InitGenotype(subPops=[0] , loci=[71], freq=[0.544,0.456]),
      sim.InitGenotype(subPops=[0] , loci=[72], freq=[0.125,0.875]),
      sim.InitGenotype(subPops=[0] , loci=[73], freq=[0.405,0.595]),
      sim.InitGenotype(subPops=[0] , loci=[74], freq=[0.432,0.568]),
      sim.InitGenotype(subPops=[0] , loci=[75], freq=[0.526,0.474]),
      sim.InitGenotype(subPops=[0] , loci=[76], freq=[0.523,0.477]),
      sim.InitGenotype(subPops=[0] , loci=[77], freq=[0.337,0.663]),
      sim.InitGenotype(subPops=[0] , loci=[78], freq=[0.558,0.442]),
      sim.InitGenotype(subPops=[0] , loci=[79], freq=[0.45,0.55]),
      sim.InitGenotype(subPops=[0] , loci=[80], freq=[0.442,0.558]),
      sim.InitGenotype(subPops=[0] , loci=[81], freq=[0.489,0.511]),
      sim.InitGenotype(subPops=[0] , loci=[82], freq=[0.723,0.277]),
      sim.InitGenotype(subPops=[0] , loci=[83], freq=[0.465,0.535]),
      sim.InitGenotype(subPops=[0] , loci=[84], freq=[0.518,0.482]),
      sim.InitGenotype(subPops=[0] , loci=[85], freq=[0.676,0.324]),
      sim.InitGenotype(subPops=[0] , loci=[86], freq=[0.591,0.409]),
      sim.InitGenotype(subPops=[0] , loci=[87], freq=[0.556,0.444]),
      sim.InitGenotype(subPops=[0] , loci=[88], freq=[0.499,0.501]),
      sim.InitGenotype(subPops=[0] , loci=[89], freq=[0.5,0.5]),
      sim.InitGenotype(subPops=[0] , loci=[90], freq=[0.487,0.513]),
      sim.InitGenotype(subPops=[0] , loci=[91], freq=[0.682,0.318]),
      sim.InitGenotype(subPops=[0] , loci=[92], freq=[0.698,0.302]),
      sim.InitGenotype(subPops=[0] , loci=[93], freq=[0.506,0.494]),
      sim.InitGenotype(subPops=[0] , loci=[94], freq=[0.656,0.344]),
      sim.InitGenotype(subPops=[0] , loci=[95], freq=[0.474,0.526]),
      sim.InitGenotype(subPops=[0] , loci=[96], freq=[0.303,0.697]),
      sim.InitGenotype(subPops=[0] , loci=[97], freq=[0.451,0.549]),
      sim.InitGenotype(subPops=[0] , loci=[98], freq=[0.362,0.638]),
      sim.InitGenotype(subPops=[0] , loci=[99], freq=[0.527,0.473]),
      sim.InitGenotype(subPops=[0] , loci=[100], freq=[0.547,0.453]),
      sim.InitGenotype(subPops=[0] , loci=[101], freq=[0.508,0.492]),
      sim.InitGenotype(subPops=[0] , loci=[102], freq=[0.481,0.519]),
      sim.InitGenotype(subPops=[0] , loci=[103], freq=[0.455,0.545]),
      sim.InitGenotype(subPops=[0] , loci=[104], freq=[0.381,0.619]),
      sim.InitGenotype(subPops=[0] , loci=[105], freq=[0.524,0.476]),
      sim.InitGenotype(subPops=[0] , loci=[106], freq=[0.479,0.521]),
      sim.InitGenotype(subPops=[0] , loci=[107], freq=[0.463,0.537]),
      sim.InitGenotype(subPops=[0] , loci=[108], freq=[0.813,0.187]),
      sim.InitGenotype(subPops=[0] , loci=[109], freq=[0.52,0.48]),
      sim.InitGenotype(subPops=[0] , loci=[110], freq=[0.558,0.442]),
      sim.InitGenotype(subPops=[0] , loci=[111], freq=[0.559,0.441]),
      sim.InitGenotype(subPops=[0] , loci=[112], freq=[0.336,0.664]),
      sim.InitGenotype(subPops=[0] , loci=[113], freq=[0.628,0.372]),
      sim.InitGenotype(subPops=[0] , loci=[114], freq=[0.741,0.259]),
      sim.InitGenotype(subPops=[0] , loci=[115], freq=[0.931,0.069]),
      sim.InitGenotype(subPops=[0] , loci=[116], freq=[0.533,0.467]),
      sim.InitGenotype(subPops=[0] , loci=[117], freq=[0.6,0.4]),
      sim.InitGenotype(subPops=[0] , loci=[118], freq=[0.519,0.481]),
      sim.InitGenotype(subPops=[0] , loci=[119], freq=[0.537,0.463]),
      sim.InitGenotype(subPops=[0] , loci=[120], freq=[0.541,0.459]),
      sim.InitGenotype(subPops=[0] , loci=[121], freq=[0.495,0.505]),
      sim.InitGenotype(subPops=[0] , loci=[122], freq=[0.232,0.768]),
      sim.InitGenotype(subPops=[0] , loci=[123], freq=[0.69,0.31]),
      sim.InitGenotype(subPops=[0] , loci=[124], freq=[0.467,0.533]),
      sim.InitGenotype(subPops=[0] , loci=[125], freq=[0.504,0.496]),
      sim.InitGenotype(subPops=[0] , loci=[126], freq=[0.483,0.517]),
      sim.InitGenotype(subPops=[0] , loci=[127], freq=[0.465,0.535]),
      sim.InitGenotype(subPops=[0] , loci=[128], freq=[0.543,0.457]),
      sim.InitGenotype(subPops=[0] , loci=[129], freq=[0.769,0.231]),
      sim.InitGenotype(subPops=[0] , loci=[130], freq=[0.492,0.508]),
      sim.InitGenotype(subPops=[0] , loci=[131], freq=[0.479,0.521]),
      sim.InitGenotype(subPops=[0] , loci=[132], freq=[0.658,0.342]),
      sim.InitGenotype(subPops=[0] , loci=[133], freq=[0.559,0.441]),
      sim.InitGenotype(subPops=[0] , loci=[134], freq=[0.469,0.531]),
      sim.InitGenotype(subPops=[0] , loci=[135], freq=[0.41,0.59]),
      sim.InitGenotype(subPops=[0] , loci=[136], freq=[0.462,0.538]),
      sim.InitGenotype(subPops=[0] , loci=[137], freq=[0.503,0.497]),
      sim.InitGenotype(subPops=[0] , loci=[138], freq=[0.415,0.585]),
      sim.InitGenotype(subPops=[0] , loci=[139], freq=[0.687,0.313]),
      sim.InitGenotype(subPops=[0] , loci=[140], freq=[0.949,0.051]),
      sim.InitGenotype(subPops=[0] , loci=[141], freq=[0.558,0.442]),
      sim.InitGenotype(subPops=[0] , loci=[142], freq=[0.501,0.499]),
      sim.InitGenotype(subPops=[0] , loci=[143], freq=[0.553,0.447]),
      sim.InitGenotype(subPops=[0] , loci=[144], freq=[0.522,0.478]),
      sim.InitGenotype(subPops=[0] , loci=[145], freq=[0.57,0.43]),
      sim.InitGenotype(subPops=[0] , loci=[146], freq=[0.297,0.703]),
      sim.InitGenotype(subPops=[0] , loci=[147], freq=[0.549,0.451]),
      sim.InitGenotype(subPops=[0] , loci=[148], freq=[0.439,0.561]),
      sim.InitGenotype(subPops=[0] , loci=[149], freq=[0.614,0.386]),
      sim.InitGenotype(subPops=[0] , loci=[150], freq=[0.476,0.524]),
      sim.InitGenotype(subPops=[0] , loci=[151], freq=[0.286,0.714]),
      sim.InitGenotype(subPops=[0] , loci=[152], freq=[0.442,0.558]),
      sim.InitGenotype(subPops=[0] , loci=[153], freq=[0.694,0.306]),
      sim.InitGenotype(subPops=[0] , loci=[154], freq=[0.48,0.52]),
      sim.InitGenotype(subPops=[0] , loci=[155], freq=[0.574,0.426]),
      sim.InitGenotype(subPops=[0] , loci=[156], freq=[0.495,0.505]),
      sim.InitGenotype(subPops=[0] , loci=[157], freq=[0.536,0.464]),
      sim.InitGenotype(subPops=[0] , loci=[158], freq=[0.344,0.656]),
      sim.InitGenotype(subPops=[0] , loci=[159], freq=[0.476,0.524]),
      sim.InitGenotype(subPops=[0] , loci=[160], freq=[0.494,0.506]),
      sim.InitGenotype(subPops=[0] , loci=[161], freq=[0.759,0.241]),
      sim.InitGenotype(subPops=[0] , loci=[162], freq=[0.475,0.525]),
      sim.InitGenotype(subPops=[0] , loci=[163], freq=[0.678,0.322]),
      sim.InitGenotype(subPops=[0] , loci=[164], freq=[0.518,0.482]),
      sim.InitGenotype(subPops=[0] , loci=[165], freq=[0.632,0.368]),
      sim.InitGenotype(subPops=[0] , loci=[166], freq=[0.479,0.521]),
      sim.InitGenotype(subPops=[0] , loci=[167], freq=[0.353,0.647]),
      sim.InitGenotype(subPops=[0] , loci=[168], freq=[0.965,0.035]),
      sim.InitGenotype(subPops=[0] , loci=[169], freq=[0.483,0.517]),
      sim.InitGenotype(subPops=[0] , loci=[170], freq=[0.7,0.3]),
      sim.InitGenotype(subPops=[0] , loci=[171], freq=[0.519,0.481]),
      sim.InitGenotype(subPops=[0] , loci=[172], freq=[0.323,0.677]),
      sim.InitGenotype(subPops=[0] , loci=[173], freq=[0.463,0.537]),
      sim.InitGenotype(subPops=[0] , loci=[174], freq=[0.56,0.44]),
      sim.InitGenotype(subPops=[0] , loci=[175], freq=[0.391,0.609]),
      sim.InitGenotype(subPops=[0] , loci=[176], freq=[0.51,0.49]),
      sim.InitGenotype(subPops=[0] , loci=[177], freq=[0.466,0.534]),
      sim.InitGenotype(subPops=[0] , loci=[178], freq=[0.599,0.401]),
      sim.InitGenotype(subPops=[0] , loci=[179], freq=[0.528,0.472]),
      sim.InitGenotype(subPops=[0] , loci=[180], freq=[0.687,0.313]),
      sim.InitGenotype(subPops=[0] , loci=[181], freq=[0.341,0.659]),
      sim.InitGenotype(subPops=[0] , loci=[182], freq=[0.579,0.421]),
      sim.InitGenotype(subPops=[0] , loci=[183], freq=[0.474,0.526]),
      sim.InitGenotype(subPops=[0] , loci=[184], freq=[0.615,0.385]),
      sim.InitGenotype(subPops=[0] , loci=[185], freq=[0.818,0.182]),
      sim.InitGenotype(subPops=[0] , loci=[186], freq=[0.404,0.596]),
      sim.InitGenotype(subPops=[0] , loci=[187], freq=[0.567,0.433]),
      sim.InitGenotype(subPops=[0] , loci=[188], freq=[0.402,0.598]),
      sim.InitGenotype(subPops=[0] , loci=[189], freq=[0.489,0.511]),
      sim.InitGenotype(subPops=[0] , loci=[190], freq=[0.537,0.463]),
      sim.InitGenotype(subPops=[0] , loci=[191], freq=[0.514,0.486]),
      sim.InitGenotype(subPops=[0] , loci=[192], freq=[0.486,0.514]),
      sim.InitGenotype(subPops=[0] , loci=[193], freq=[0.565,0.435]),
      sim.InitGenotype(subPops=[0] , loci=[194], freq=[0.543,0.457]),
      sim.InitGenotype(subPops=[0] , loci=[195], freq=[0.633,0.367]),
      sim.InitGenotype(subPops=[0] , loci=[196], freq=[0.548,0.452]),
      sim.InitGenotype(subPops=[0] , loci=[197], freq=[0.458,0.542]),
      sim.InitGenotype(subPops=[0] , loci=[198], freq=[0.612,0.388]),
      sim.InitGenotype(subPops=[0] , loci=[199], freq=[0.572,0.428]),
      sim.InitGenotype(subPops=[1] , loci=[0], freq=[0.054,0.946]),
      sim.InitGenotype(subPops=[1] , loci=[1], freq=[0.983,0.017]),
      sim.InitGenotype(subPops=[1] , loci=[2], freq=[0.018,0.982]),
      sim.InitGenotype(subPops=[1] , loci=[3], freq=[0.013,0.987]),
      sim.InitGenotype(subPops=[1] , loci=[4], freq=[0.923,0.077]),
      sim.InitGenotype(subPops=[1] , loci=[5], freq=[0.976,0.024]),
      sim.InitGenotype(subPops=[1] , loci=[6], freq=[0.982,0.018]),
      sim.InitGenotype(subPops=[1] , loci=[7], freq=[0.023,0.977]),
      sim.InitGenotype(subPops=[1] , loci=[8], freq=[0.973,0.027]),
      sim.InitGenotype(subPops=[1] , loci=[9], freq=[0.991,0.009]),
      sim.InitGenotype(subPops=[1] , loci=[10], freq=[0.879,0.121]),
      sim.InitGenotype(subPops=[1] , loci=[11], freq=[0.819,0.181]),
      sim.InitGenotype(subPops=[1] , loci=[12], freq=[0.022,0.978]),
      sim.InitGenotype(subPops=[1] , loci=[13], freq=[0.953,0.047]),
      sim.InitGenotype(subPops=[1] , loci=[14], freq=[0.03,0.97]),
      sim.InitGenotype(subPops=[1] , loci=[15], freq=[0.98,0.02]),
      sim.InitGenotype(subPops=[1] , loci=[16], freq=[0.94,0.06]),
      sim.InitGenotype(subPops=[1] , loci=[17], freq=[0.992,0.008]),
      sim.InitGenotype(subPops=[1] , loci=[18], freq=[0.76,0.24]),
      sim.InitGenotype(subPops=[1] , loci=[19], freq=[0.039,0.961]),
      sim.InitGenotype(subPops=[1] , loci=[20], freq=[0.929,0.071]),
      sim.InitGenotype(subPops=[1] , loci=[21], freq=[0.146,0.854]),
      sim.InitGenotype(subPops=[1] , loci=[22], freq=[0.97,0.03]),
      sim.InitGenotype(subPops=[1] , loci=[23], freq=[0.01,0.99]),
      sim.InitGenotype(subPops=[1] , loci=[24], freq=[0.907,0.093]),
      sim.InitGenotype(subPops=[1] , loci=[25], freq=[0.978,0.022]),
      sim.InitGenotype(subPops=[1] , loci=[26], freq=[0.986,0.014]),
      sim.InitGenotype(subPops=[1] , loci=[27], freq=[0.948,0.052]),
      sim.InitGenotype(subPops=[1] , loci=[28], freq=[0.021,0.979]),
      sim.InitGenotype(subPops=[1] , loci=[29], freq=[0.018,0.982]),
      sim.InitGenotype(subPops=[1] , loci=[30], freq=[0.984,0.016]),
      sim.InitGenotype(subPops=[1] , loci=[31], freq=[0.984,0.016]),
      sim.InitGenotype(subPops=[1] , loci=[32], freq=[0.034,0.966]),
      sim.InitGenotype(subPops=[1] , loci=[33], freq=[0.895,0.105]),
      sim.InitGenotype(subPops=[1] , loci=[34], freq=[0.97,0.03]),
      sim.InitGenotype(subPops=[1] , loci=[35], freq=[0.978,0.022]),
      sim.InitGenotype(subPops=[1] , loci=[36], freq=[0.022,0.978]),
      sim.InitGenotype(subPops=[1] , loci=[37], freq=[0.981,0.019]),
      sim.InitGenotype(subPops=[1] , loci=[38], freq=[0.815,0.185]),
      sim.InitGenotype(subPops=[1] , loci=[39], freq=[0.106,0.894]),
      sim.InitGenotype(subPops=[1] , loci=[40], freq=[0.131,0.869]),
      sim.InitGenotype(subPops=[1] , loci=[41], freq=[0.976,0.024]),
      sim.InitGenotype(subPops=[1] , loci=[42], freq=[0.966,0.034]),
      sim.InitGenotype(subPops=[1] , loci=[43], freq=[0.043,0.957]),
      sim.InitGenotype(subPops=[1] , loci=[44], freq=[0.977,0.023]),
      sim.InitGenotype(subPops=[1] , loci=[45], freq=[0.895,0.105]),
      sim.InitGenotype(subPops=[1] , loci=[46], freq=[0.019,0.981]),
      sim.InitGenotype(subPops=[1] , loci=[47], freq=[0.979,0.021]),
      sim.InitGenotype(subPops=[1] , loci=[48], freq=[0.986,0.014]),
      sim.InitGenotype(subPops=[1] , loci=[49], freq=[0.017,0.983]),
      sim.InitGenotype(subPops=[1] , loci=[50], freq=[0.878,0.122]),
      sim.InitGenotype(subPops=[1] , loci=[51], freq=[0.979,0.021]),
      sim.InitGenotype(subPops=[1] , loci=[52], freq=[0.407,0.593]),
      sim.InitGenotype(subPops=[1] , loci=[53], freq=[0.027,0.973]),
      sim.InitGenotype(subPops=[1] , loci=[54], freq=[0.889,0.111]),
      sim.InitGenotype(subPops=[1] , loci=[55], freq=[0.986,0.014]),
      sim.InitGenotype(subPops=[1] , loci=[56], freq=[0.01,0.99]),
      sim.InitGenotype(subPops=[1] , loci=[57], freq=[0.962,0.038]),
      sim.InitGenotype(subPops=[1] , loci=[58], freq=[0.039,0.961]),
      sim.InitGenotype(subPops=[1] , loci=[59], freq=[0.046,0.954]),
      sim.InitGenotype(subPops=[1] , loci=[60], freq=[0.128,0.872]),
      sim.InitGenotype(subPops=[1] , loci=[61], freq=[0.964,0.036]),
      sim.InitGenotype(subPops=[1] , loci=[62], freq=[0.022,0.978]),
      sim.InitGenotype(subPops=[1] , loci=[63], freq=[0.959,0.041]),
      sim.InitGenotype(subPops=[1] , loci=[64], freq=[0.215,0.785]),
      sim.InitGenotype(subPops=[1] , loci=[65], freq=[0.094,0.906]),
      sim.InitGenotype(subPops=[1] , loci=[66], freq=[0.994,0.006]),
      sim.InitGenotype(subPops=[1] , loci=[67], freq=[0.925,0.075]),
      sim.InitGenotype(subPops=[1] , loci=[68], freq=[0.868,0.132]),
      sim.InitGenotype(subPops=[1] , loci=[69], freq=[0.046,0.954]),
      sim.InitGenotype(subPops=[1] , loci=[70], freq=[0.815,0.185]),
      sim.InitGenotype(subPops=[1] , loci=[71], freq=[0.018,0.982]),
      sim.InitGenotype(subPops=[1] , loci=[72], freq=[0.016,0.984]),
      sim.InitGenotype(subPops=[1] , loci=[73], freq=[0.976,0.024]),
      sim.InitGenotype(subPops=[1] , loci=[74], freq=[0.961,0.039]),
      sim.InitGenotype(subPops=[1] , loci=[75], freq=[0.966,0.034]),
      sim.InitGenotype(subPops=[1] , loci=[76], freq=[0.945,0.055]),
      sim.InitGenotype(subPops=[1] , loci=[77], freq=[0.977,0.023]),
      sim.InitGenotype(subPops=[1] , loci=[78], freq=[0.965,0.035]),
      sim.InitGenotype(subPops=[1] , loci=[79], freq=[0.087,0.913]),
      sim.InitGenotype(subPops=[1] , loci=[80], freq=[0.984,0.016]),
      sim.InitGenotype(subPops=[1] , loci=[81], freq=[0.106,0.894]),
      sim.InitGenotype(subPops=[1] , loci=[82], freq=[0.97,0.03]),
      sim.InitGenotype(subPops=[1] , loci=[83], freq=[0.932,0.068]),
      sim.InitGenotype(subPops=[1] , loci=[84], freq=[0.97,0.03]),
      sim.InitGenotype(subPops=[1] , loci=[85], freq=[0.989,0.011]),
      sim.InitGenotype(subPops=[1] , loci=[86], freq=[0.961,0.039]),
      sim.InitGenotype(subPops=[1] , loci=[87], freq=[0.981,0.019]),
      sim.InitGenotype(subPops=[1] , loci=[88], freq=[0.063,0.937]),
      sim.InitGenotype(subPops=[1] , loci=[89], freq=[0.95,0.05]),
      sim.InitGenotype(subPops=[1] , loci=[90], freq=[0.933,0.067]),
      sim.InitGenotype(subPops=[1] , loci=[91], freq=[0.961,0.039]),
      sim.InitGenotype(subPops=[1] , loci=[92], freq=[0.98,0.02]),
      sim.InitGenotype(subPops=[1] , loci=[93], freq=[0.03,0.97]),
      sim.InitGenotype(subPops=[1] , loci=[94], freq=[0.067,0.933]),
      sim.InitGenotype(subPops=[1] , loci=[95], freq=[0.045,0.955]),
      sim.InitGenotype(subPops=[1] , loci=[96], freq=[0.905,0.095]),
      sim.InitGenotype(subPops=[1] , loci=[97], freq=[0.98,0.02]),
      sim.InitGenotype(subPops=[1] , loci=[98], freq=[0.931,0.069]),
      sim.InitGenotype(subPops=[1] , loci=[99], freq=[0.853,0.147]),
      sim.InitGenotype(subPops=[1] , loci=[100], freq=[0.989,0.011]),
      sim.InitGenotype(subPops=[1] , loci=[101], freq=[0.07,0.93]),
      sim.InitGenotype(subPops=[1] , loci=[102], freq=[0.981,0.019]),
      sim.InitGenotype(subPops=[1] , loci=[103], freq=[0.907,0.093]),
      sim.InitGenotype(subPops=[1] , loci=[104], freq=[0.024,0.976]),
      sim.InitGenotype(subPops=[1] , loci=[105], freq=[0.03,0.97]),
      sim.InitGenotype(subPops=[1] , loci=[106], freq=[0.952,0.048]),
      sim.InitGenotype(subPops=[1] , loci=[107], freq=[0.01,0.99]),
      sim.InitGenotype(subPops=[1] , loci=[108], freq=[0.636,0.364]),
      sim.InitGenotype(subPops=[1] , loci=[109], freq=[0.18,0.82]),
      sim.InitGenotype(subPops=[1] , loci=[110], freq=[0.032,0.968]),
      sim.InitGenotype(subPops=[1] , loci=[111], freq=[0.991,0.009]),
      sim.InitGenotype(subPops=[1] , loci=[112], freq=[0.981,0.019]),
      sim.InitGenotype(subPops=[1] , loci=[113], freq=[0.062,0.938]),
      sim.InitGenotype(subPops=[1] , loci=[114], freq=[0.99,0.01]),
      sim.InitGenotype(subPops=[1] , loci=[115], freq=[0.28,0.72]),
      sim.InitGenotype(subPops=[1] , loci=[116], freq=[0.969,0.031]),
      sim.InitGenotype(subPops=[1] , loci=[117], freq=[0.98,0.02]),
      sim.InitGenotype(subPops=[1] , loci=[118], freq=[0.938,0.062]),
      sim.InitGenotype(subPops=[1] , loci=[119], freq=[0.033,0.967]),
      sim.InitGenotype(subPops=[1] , loci=[120], freq=[0.947,0.053]),
      sim.InitGenotype(subPops=[1] , loci=[121], freq=[0.046,0.954]),
      sim.InitGenotype(subPops=[1] , loci=[122], freq=[0.855,0.145]),
      sim.InitGenotype(subPops=[1] , loci=[123], freq=[0.972,0.028]),
      sim.InitGenotype(subPops=[1] , loci=[124], freq=[0.094,0.906]),
      sim.InitGenotype(subPops=[1] , loci=[125], freq=[0.089,0.911]),
      sim.InitGenotype(subPops=[1] , loci=[126], freq=[0.939,0.061]),
      sim.InitGenotype(subPops=[1] , loci=[127], freq=[0.046,0.954]),
      sim.InitGenotype(subPops=[1] , loci=[128], freq=[0.982,0.018]),
      sim.InitGenotype(subPops=[1] , loci=[129], freq=[0.316,0.684]),
      sim.InitGenotype(subPops=[1] , loci=[130], freq=[0.818,0.182]),
      sim.InitGenotype(subPops=[1] , loci=[131], freq=[0.951,0.049]),
      sim.InitGenotype(subPops=[1] , loci=[132], freq=[0.979,0.021]),
      sim.InitGenotype(subPops=[1] , loci=[133], freq=[0.958,0.042]),
      sim.InitGenotype(subPops=[1] , loci=[134], freq=[0.889,0.111]),
      sim.InitGenotype(subPops=[1] , loci=[135], freq=[0.012,0.988]),
      sim.InitGenotype(subPops=[1] , loci=[136], freq=[0.039,0.961]),
      sim.InitGenotype(subPops=[1] , loci=[137], freq=[0.036,0.964]),
      sim.InitGenotype(subPops=[1] , loci=[138], freq=[0.988,0.012]),
      sim.InitGenotype(subPops=[1] , loci=[139], freq=[0.081,0.919]),
      sim.InitGenotype(subPops=[1] , loci=[140], freq=[0.392,0.608]),
      sim.InitGenotype(subPops=[1] , loci=[141], freq=[0.033,0.967]),
      sim.InitGenotype(subPops=[1] , loci=[142], freq=[0.864,0.136]),
      sim.InitGenotype(subPops=[1] , loci=[143], freq=[0.903,0.097]),
      sim.InitGenotype(subPops=[1] , loci=[144], freq=[0.047,0.953]),
      sim.InitGenotype(subPops=[1] , loci=[145], freq=[0.95,0.05]),
      sim.InitGenotype(subPops=[1] , loci=[146], freq=[0.019,0.981]),
      sim.InitGenotype(subPops=[1] , loci=[147], freq=[0.18,0.82]),
      sim.InitGenotype(subPops=[1] , loci=[148], freq=[0.029,0.971]),
      sim.InitGenotype(subPops=[1] , loci=[149], freq=[0.958,0.042]),
      sim.InitGenotype(subPops=[1] , loci=[150], freq=[0.023,0.977]),
      sim.InitGenotype(subPops=[1] , loci=[151], freq=[0.676,0.324]),
      sim.InitGenotype(subPops=[1] , loci=[152], freq=[0.817,0.183]),
      sim.InitGenotype(subPops=[1] , loci=[153], freq=[0.979,0.021]),
      sim.InitGenotype(subPops=[1] , loci=[154], freq=[0.973,0.027]),
      sim.InitGenotype(subPops=[1] , loci=[155], freq=[0.967,0.033]),
      sim.InitGenotype(subPops=[1] , loci=[156], freq=[0.038,0.962]),
      sim.InitGenotype(subPops=[1] , loci=[157], freq=[0.065,0.935]),
      sim.InitGenotype(subPops=[1] , loci=[158], freq=[0.795,0.205]),
      sim.InitGenotype(subPops=[1] , loci=[159], freq=[0.988,0.012]),
      sim.InitGenotype(subPops=[1] , loci=[160], freq=[0.89,0.11]),
      sim.InitGenotype(subPops=[1] , loci=[161], freq=[0.028,0.972]),
      sim.InitGenotype(subPops=[1] , loci=[162], freq=[0.909,0.091]),
      sim.InitGenotype(subPops=[1] , loci=[163], freq=[0.968,0.032]),
      sim.InitGenotype(subPops=[1] , loci=[164], freq=[0.948,0.052]),
      sim.InitGenotype(subPops=[1] , loci=[165], freq=[0.039,0.961]),
      sim.InitGenotype(subPops=[1] , loci=[166], freq=[0.969,0.031]),
      sim.InitGenotype(subPops=[1] , loci=[167], freq=[0.037,0.963]),
      sim.InitGenotype(subPops=[1] , loci=[168], freq=[0.485,0.515]),
      sim.InitGenotype(subPops=[1] , loci=[169], freq=[0.972,0.028]),
      sim.InitGenotype(subPops=[1] , loci=[170], freq=[0.25,0.75]),
      sim.InitGenotype(subPops=[1] , loci=[171], freq=[0.936,0.064]),
      sim.InitGenotype(subPops=[1] , loci=[172], freq=[0.014,0.986]),
      sim.InitGenotype(subPops=[1] , loci=[173], freq=[0.954,0.046]),
      sim.InitGenotype(subPops=[1] , loci=[174], freq=[0.316,0.684]),
      sim.InitGenotype(subPops=[1] , loci=[175], freq=[0.034,0.966]),
      sim.InitGenotype(subPops=[1] , loci=[176], freq=[0.978,0.022]),
      sim.InitGenotype(subPops=[1] , loci=[177], freq=[0.012,0.988]),
      sim.InitGenotype(subPops=[1] , loci=[178], freq=[0.975,0.025]),
      sim.InitGenotype(subPops=[1] , loci=[179], freq=[0.915,0.085]),
      sim.InitGenotype(subPops=[1] , loci=[180], freq=[0.966,0.034]),
      sim.InitGenotype(subPops=[1] , loci=[181], freq=[0.944,0.056]),
      sim.InitGenotype(subPops=[1] , loci=[182], freq=[0.976,0.024]),
      sim.InitGenotype(subPops=[1] , loci=[183], freq=[0.012,0.988]),
      sim.InitGenotype(subPops=[1] , loci=[184], freq=[0.04,0.96]),
      sim.InitGenotype(subPops=[1] , loci=[185], freq=[0.036,0.964]),
      sim.InitGenotype(subPops=[1] , loci=[186], freq=[0.044,0.956]),
      sim.InitGenotype(subPops=[1] , loci=[187], freq=[0.948,0.052]),
      sim.InitGenotype(subPops=[1] , loci=[188], freq=[0.02,0.98]),
      sim.InitGenotype(subPops=[1] , loci=[189], freq=[0.034,0.966]),
      sim.InitGenotype(subPops=[1] , loci=[190], freq=[0.989,0.011]),
      sim.InitGenotype(subPops=[1] , loci=[191], freq=[0.022,0.978]),
      sim.InitGenotype(subPops=[1] , loci=[192], freq=[0.107,0.893]),
      sim.InitGenotype(subPops=[1] , loci=[193], freq=[0.97,0.03]),
      sim.InitGenotype(subPops=[1] , loci=[194], freq=[0.963,0.037]),
      sim.InitGenotype(subPops=[1] , loci=[195], freq=[0.044,0.956]),
      sim.InitGenotype(subPops=[1] , loci=[196], freq=[0.025,0.975]),
      sim.InitGenotype(subPops=[1] , loci=[197], freq=[0.914,0.086]),
      sim.InitGenotype(subPops=[1] , loci=[198], freq=[0.982,0.018]),
      sim.InitGenotype(subPops=[1] , loci=[199], freq=[0.983,0.017])

   ],

   preOps=sim.Migrator(rate=[100], mode=sim.BY_COUNTS, subPops=[0],toSubPops=[1]),
   matingScheme=sim.HomoMating(
      sim.RandomParentsChooser(),
      sim.OffspringGenerator(ops=[
         sim.Recombinator(intensity=0.01),
         sim.IdTagger(),
         sim.PedigreeTagger()
      ]),
      subPopSize=[1000,1000]
   ),
   gen=3
)

sim.utils.export(Inari, format='PED', output='Inarijoki_beta.ped')