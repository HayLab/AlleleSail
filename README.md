alleleSail_sim is a discrete time & generation stochastic simulations, used for modeling the behavior of allele sails for both modification and suppression.

alleleSail_sim was built in Python3.7, and requirements can be found in requirements.txt

## How to use

More information on how to use this simulation can be found in the demo: [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/HayLab/AlleleSail/blob/main/AlleleSail_demo.ipynb)


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