import numpy as np
from tqdm import tqdm
import msprime
import demesdraw
from matplotlib import pyplot as plt
import tskit

# add demography

def repeat_simulations(mut, length, reco, num_simulations, seed=None):
    results = []
    Fst =[]
    dxy = []
    segsites1 = []
    segsites2 = []
    pi1 = []
    pi2 = []
    tajima1 = []
    tajima2 = []
    for i in tqdm(range(num_simulations), desc="Running simulations"):
        if seed is not None:
          np.random.seed(seed + i)
          demography = msprime.Demography()
          demography.add_population(name="N1", initial_size=150_000)
          demography.add_population(name="N2", initial_size=5_000)
          demography.add_population(name="ANC", initial_size=7_000_000)

        
          mu, sigma = 4000, 800
          Time = np.random.normal(mu, sigma)
          if Time > 8000 : 
            Time = 8000
          if Time < 1000 :
            Time = 1000
        
          demography.add_population_split(time=Time, derived=["N1", "N2"], ancestral="ANC") # split 2k generations ago
        
        # Simulate 50 diploid samples under the coalescent with recombination on a 1kb region.
          ts = msprime.sim_ancestry(
          {"N1": 50, "N2": 50},
          demography=demography,
          recombination_rate=reco, # as in humans
          sequence_length=1_000,
          random_seed=1234)
        # we can add mutations
          mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))
          diversity = mutated_ts.diversity()
          tajimas_d = mutated_ts.Tajimas_D()
          allele_frequency_spectrum = mutated_ts.allele_frequency_spectrum(polarised=True)
        #results.append((mutated_ts, None, diversity, tajimas_d, allele_frequency_spectrum))
        
          pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
          sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]
        
          Fst.append(mts.Fst(sample_sets))
        
    return results


mut= 3.5e-9
length = 1_000
seed = 4711
reco = 8.4e-9
num_simulations = 100

results = repeat_simulations(mut, sample_sizes, length, reco, pop_size, num_simulations, seed=seed)


# plot
diversities = [result[2] for result in results]
tajimas_ds = [result[3] for result in results]
allele_frequency_spectra = [result[4] for result in results]
plt.clf()
plt.figure(figsize=(10, 5))
plt.hist(diversities, bins=10, color='skyblue', edgecolor='black', alpha=0.7)
plt.xlabel("Nucleotide Diversity (π)")
plt.ylabel("Frequency")
plt.title("Histogram of Nucleotide Diversity Across Simulations")
plt.show()


