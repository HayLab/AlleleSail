alleleSail_sim is a python commandline program that runs discrete time & generation stochastic simulations, used for modeling the behavior of allele sails for both modification and suppression.

alleleSail_sim was built in Python3.7, and has been tested to run on Python3.8 and 3.10. Module requirements can be found in requirements.txt

To analyze the data, we wrote our own R scripts - the data we analyzed can be found in the data folder, and figures were generated using the R markdown script figures_4-4.Rmd.  

## How to use

More information on how to use this simulation can be found in the demo: [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/HayLab/AlleleSail/blob/main/AlleleSail_demo.ipynb)

As a quick tl;dr you can clone this github repository (which should take a few minutes), and simply run ```python alleleSail_sim.py -h``` to see the command line arguments that the simulation takes in. Data is written out to five files, ```_genotype``` which includes the number of individuals of each genotype, ```_allele``` which includes the number of individuals that carry each allele, ```_NEWallele``` which includes the absolute number of each allele in the popuolation, ```_total``` which includes the total population including transgenic additions, and ```_total_pop``` which includes total population WITHOUT including the transgenic additions.

## Data Files

The data files here are used to generate the figures in the paper. They are stored by date of generation, and the contents are as follows.

| folder | files | data | figure in paper |
| :----- | :--- | :--- | :-------------- |
| 1-11 | noSail_modification_MC ... | release of individuals with the "edit" but no "editor", for introduction frequencies between 0 and 0.5, and for editor costs between -0.05 (a benefit) per allele and 0.1 (cost) per allele. | used in Figure 2, Figure 3, Supplemental Figure 3, and Supplemental Figure 4.
| 1-11 | EvI/introFreq_{value} ... | release of individuals carrying both edit and editor, for the given introduction frequency, for fitness costs on the Edit, ranging between -0.05 and 0.2 per allele | used in figure 2
| 1-11 | EvI/last_gen_E_carrier | information pulled from multiple introduction frequency _allele files that gives the number of carriers of the edit "E" at generation 50 (the last generation of the simulation). For various fitness costs on the Edit | used in heatmaps of figure 2
| 1-11 | EvI/last_gen_O_carrier | information pulled from multiple introduction frequency _allele files that gives the number of carriers of the NOT-edit "O" at generation 50 (the last generation of the simulation). For various fitness costs on the Edit | used in heatmaps of figure 2 (homozygotes = total population - individuals that carry the wildtype)
| 1-11 | SvI/last_gen_E_carrier | information pulled from multiple introduction frequency _allele files that gives the number of carriers of the edit "E" at generation 50 (the last generation of the simulation). For various fitness costs on the Editor | used in heatmaps of figure 2
| 1-11 | SvI/last_gen_O_carrier | information pulled from multiple introduction frequency _allele files that gives the number of carriers of the NOT-edit "O" at generation 50 (the last generation of the simulation). For various fitness costs on the Editor | used in heatmaps of figure 2 (homozygotes = total population - individuals that carry the wildtype)
| 1-19 | sex_equilibriums... | 500 simulations of either XY or ZW sex systems, with edit/edit individuals released at various introduction frequencies | figure 6
| 2-5 | low_efficiencies_MC ... | population modification simulations where our editor has lower efficiencies, and maternal carryover of editing occurs. | Figure 3
| 2-5 | low_efficiencies_noMC ... | population modification simulations where our editor has lower efficiencies, and maternal carryover does not occur. | Figure 3
| 2-5 | rd_0-clvg_0.5_last_gen_E | modification simulations with 50% cleavage, and a recombination distance of zero (linked editor and edit). For various introduction frequencies and various fitness costs on the edit. | Figure 4
| 2-5 | rd_50-clvg_0.5_last_gen_E | modification simulations with 50% cleavage, and a maximal recombination distance (unlinked editor and edit). For various introduction frequencies and various fitness costs on the edit. | Figure 4
| 2-5 | svi_rd_0-clvg_0.5_last_gen_E | modification simulations with 50% cleavage, and a recombination distance of zero (linked editor and edit). For various introduction frequencies and various fitness costs on the editor. | Supplemental Figure 5
| 2-5 | svi_rd_50-clvg_0.5_last_gen_E | modification simulations with 50% cleavage, and a maximal recombination distance (unlinked editor and edit). For various introduction frequencies and various fitness costs on the editor. | Supplemental Figure 5
| 2-23 | ais_51releases_noMC... | continuous releases of aromatase cleavers / sex skew individuals, for various systems (XY, ZW, with and without sterililty, with different introductions of XX or XY or ZZ or ZW). | Figure 7, Supplemental Figure 9, Supplemental Figure 7.
| 2-23 | ais_ZW_ZWind_noMC ... | continuous releases of aromatase cleavers / sex skew individuals, where the cleaver/editor is on the Z chromosome | Supplemental Figure 7
| 2-23 | low_efficiencies_int_{value}_MC_...| population modification simulations where the editor has editing efficiencies of 15%, 50%, or 100%, maternal carryover occurs, and the release occurs at an introduction of {value}%. | Supplemental Figure 4 
| 2-23 | low_efficiencies_int_{value}_noMC_...| population modification simulations where the editor has editing efficiencies of 15%, 50%, or 100%, maternal carryover does not occur, and the release occurs at an introduction of {value}%. | Supplemental Figure 4
| 2-23 | 51releases_highEfficiencies_2... | name is slightly misleading - these files contain continuous releases of aromatase cleavers that have cleavage efficiency of less than 100 - the name compares to "lowEfficieny" modification runs that go as low as 15% efficiency, where these range from 80% to 100%. For XY, ZW, and fsRIDL systems | Supplemental Figure 9
| 3-5 | full_five_noMC_... | suppression simulations with only a single release, for XY releasing XX, XY releasing XY, ZW WW non-Viable, ZW WW-viable releasing ZW, ZW WW-viable releasing ZZ, and fsRIDL. | Used in Figure 5, Supplemental Figure 6, Supplemental Figure 7, and Supplemental Figure 8
| 3-11 | Sail_modification_MC_{fitness type}... | For various introduction frequencies and fitness costs on the editor, with maternal carryover occurring, population modification simulations where fitness costs are applied using the listed {fitness type} | Supplemental Figure 1
| 3-11 | noSail_modification_MC_{fitness type}... | The same as above, but individuals released do not carry the editor | Supplemental Figure 1
| 3-11 | ZW_Zwind_SingleRelease_noMC... | single releases of an aromatase cleaver attached to the Z chromosome | Supplemental Figure 7
| 3-28 | DAR_{modification type}_... | simulations for dominant, additive, and recessive costs on the editor, for the given {modification type} (somatic or germline editing). No Maternal carryover occurs. | Supplemental Figure 2
| 3-28 | DAR_{modification type}MC_... | simulations for dominant, additive, and recessive costs on the editor, for the given {modification type} (somatic or germline editing). Maternal carryover occurs. | Supplemental Figure 2

## Assumptions

**1) A panmictic population.** All individuals have an equal chance of mating with all other individuals

**2) All females mate.** If there is a sufficient number of males, each female will reproduce. Population is therefore dependent on the number of females

**3) Mating is monogamous.** If a female mates with a male, both individuals are removed from the breeding population.

**4) Population growth is logistic, dependent on adult population density.** Specifically, expected population follows the Beverton-Holt Model. When the population size is at carrying capacity, each individual mating will have two offspring that survive. When the poulation size is low, more or all offspring will survive. A very small population will experience logistic growth. *note: species such as Anopheles are more dependent on larval-density. As such, the suppression trends modeled here may not accurately represent behavior in some species, but comparisons should hold between suppression systems. More on this can be read below*

**5) Generations do not overlap.** Adults will mate, and produce offspring. All adults are then removed from the population, and offspring mate only within their generation. 

**6) Released individuals do not affect density-dependence.** We assume individuals are released immediately before mating. If we assume density-dependence occurs because individuals use up resources as they grow to adulthood, then the release of lab-grown adults right before mating should not drain resources, and therefore do not affect our measure of density, in terms of density-dependent growth.

## Simulation Structure and Details

This is an agent-based simulation using discrete generations to model the behavior of an allele sail over time. The ```main``` function takes in command-line parser arguments, and formats them so they can be passed onto the function ```run_stochastic_sim```.

```run_stochastic_sim ``` does some more processing of the inputs, and then runs a simulation. The simulation starts at time/generation 0, with a pool of females and a pool of males. For each female, we randomly select a male to father her children. Their genotypes are crossed, and probabilities are assigned to the likelihood of each offspring, depending on recombination rates and modifications that occur. Modifications and their probabilites are handled by the function ```drive_modification_stochastic```. The resulting possible offspring and their probabilities are stored in a dictionary, for faster access later. The number of offspring to create is pulled from a poisson distribution, where the expected value is the pre-defined number of offspring per female, multiplied by the mother and father's fertilities (which are between 0 and 1, 0 meaning sterile and 1 meaning completely fertile). We then randomly choose that number of offspring from the given possible offspring, weighted by the probabilities. The amount of individuals is then culled based on a survival_modifier (more on this below) to generate density-dependent growth. If sex determination is not specified (i.e., specified XY or ZW), then surviving offspring are randomly assigned to be male or female. Once all surviving offspring from all females have been generated, they are stored in separate male and female pools, and are used as the parents for the next generation.

### Density-Dependent Growth

According to [Wikipedia](https://en.wikipedia.org/wiki/Beverton%E2%80%93Holt_model) *(Accessed April 5, 2024)*:

 > The Beverton–Holt model is a classic discrete-time population model which gives the expected number $n_{t+1}$ (or density) of individuals in generation $t + 1$ as a function of the number of individuals in the previous generation, 
 > $$n_{t+1}={\frac {R_{0}n_{t}}{1+n_{t}/M}}$$
 > Here $R_0$ is interpreted as the proliferation rate per generation and $K = (R_0 − 1) M$ is the carrying capacity of the environment. 

We want an equation that can modify the number of offspring *per mating*, depending on the number of individuals in the population and $K$. As such, we divide by $n_t$ and rewrite in terms of $K$, and get some survival modifier $m$ where $$m = \frac{R_0}{1 + (R_0 - 1)\frac{n_t}{K}}$$ Here again $R_0$ represents some proliferation rate, $n_t$ represents the size of the previous generation, and $K$ is the carrying capacity.

When the population $n_t$ is at or near carrying capacity $K$, this modifier will be equal to 1. When the population is near zero, then this modifier will be equal to $R_0$. However, at carrying capacity we do not want every offspring of an individual to survive; for replacement, we want $2/n_t$ offspring per mating. We must also consider an offspring's viability $\omega_i$. As such, the chance of survival for some offspring $i$ is as follows.

$$P_{survival}(i) = \omega_i \cdot m \cdot \frac{2}{number\_ offspring}$$

Of note, the density-dependent modifier depends on the number of adults in the previous generation. This matches with the biological assumption that adult competition for resources affects the number of offspring that a mating can produce and raise to adulthood. This assumption matches well for capital breeders (organisms that used stored resources to support reproduction) or organisms with high levels or partental care. However, in species such as mosquitos which have no parental care, and are largely limited by larval competition, this assumption is a poor match. That said, we expect comparisons between different sytems (such as fsRIDL vs. Allele Sail) to hold for various population densities.
