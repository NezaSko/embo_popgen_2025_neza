import msprime
import demesdraw
from matplotlib import pyplot as plt

# add demography
demography = msprime.Demography()
demography.add_population(name="N1", initial_size=50_000)
demography.add_population(name="N2", initial_size=1_667)

# instantaneous reduction of size 
#demography.add_population_parameters_change(population="EUR", time=100, initial_size=7_000_000)

demography.add_population(name="ANC", initial_size=7_000_000)
demography.add_population_split(time=1000, derived=["N1", "N2"], ancestral="ANC") # split 1k generations ago

# instatanous growth
#demography.add_population_parameters_change(time=4000, population="ANC", initial_size=100_000)

# instantaneous bottleneck
print(demography)

# Plot a schematic of the model
plt.clf()
demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.show()


def repeat_simulations(mut, length, reco, demography, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 
        # Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
        ts = msprime.sim_ancestry(
            {"N1" : 10, "N2" : 10},
            recombination_rate=reco,
            sequence_length=length,
            demography=demography,
            random_seed=np.random.randint(99999999))
        
        # we can add mutations
        mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))

        diversity = mutated_ts.diversity()
        tajimas_d = mutated_ts.Tajimas_D()
        allele_frequency_spectrum = mutated_ts.allele_frequency_spectrum(polarised=True)
        results.append((mutated_ts, None, diversity, tajimas_d, allele_frequency_spectrum))
    return results

mut= 3.5e-9
#sample_sizes = [20]
length = 1_000
seed = 4711
reco = 8.4e-9
num_simulations = 100

results = repeat_simulations(mut, length, reco, demography, num_simulations, seed=seed)

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


plt.clf()
plt.figure(figsize=(10, 5))
plt.hist(tajimas_ds, bins=10, color='pink', edgecolor='black', alpha=0.7)
plt.xlabel("Tajima's D")
plt.ylabel("Frequency")
plt.title("Distribution of Tajima's D Across Simulations")
plt.show()


plt.clf()
plt.figure(figsize=(10, 5))
bar_width = 0.8 / num_simulations 
colors = plt.cm.tab20(np.linspace(0, 1, num_simulations)) 
for i, afs in enumerate(allele_frequency_spectra):
    x_positions = np.arange(len(afs)) + i * bar_width  
    plt.bar(x_positions, afs, width=bar_width, color=colors[i], label=f'Simulation {i+1}')

plt.xlabel("Frequency")
plt.ylabel("Number of Sites")
plt.title("Allele Frequency Spectrum Across Simulations")
plt.show()

combined_afs = np.sum(allele_frequency_spectra, axis=0)
normalized_afs = combined_afs / np.sum(combined_afs)

plt.clf()
plt.figure(figsize=(10, 5))
plt.bar(range(len(normalized_afs)), normalized_afs, color='green', edgecolor='black', alpha=0.7)
plt.xlabel("Derived Allele Frequency")
plt.ylabel("Proportion of Sites")
plt.title("Allele Frequency Spectrum Across Simulations")
plt.show()

# Summary stats
print(f"Nucleotide Diversity: mean = {np.mean(diversities)}, std = {np.std(diversities)}")
print(f"Tajima's D: mean = {np.mean(tajimas_ds)}, std = {np.std(tajimas_ds)}")
