alleleSail_sim is a discrete time & generation stochastic simulations, used for modeling the behavior of allele sails for both modification and suppression.

alleleSail_sim was built in Python3.7, and requirements can be found in requirements.txt

## How to use

More information on how to use this simulation can be found in the demo: [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/HayLab/AlleleSail/blob/main/AlleleSail_demo.ipynb)


## Assumptions

**1) A panmictic population.** All individuals have an equal chance of mating with all other individuals

**2) All females mate.** If there is a sufficient number of males, each female will reproduce. Population is therefore dependent on the number of females

**3) Mating is monogamous.** If a female mates with a male, both individuals are removed from the breeding population.

**4) Population growth is logistic, dependent on adult population density.** Specifically, expected population follows the Beverton-Holt Model. When the population size is at carrying capacity, each individual mating will have two offspring that survive. When the poulation size is low, more or all offspring will survive. A very small population will experience logistic growth. *note: species such as Anopheles are more dependent on larval-density. As such, the suppression trends modeled here may not accurately represent behavior in some species, but comparisons should hold between suppression systems.*

**5) Generations do not overlap.** Adults will mate, and produce offspring. All adults are then removed from the population, and offspring mate only within their generation. 

**6) Released individuals do not affect density-dependence.** We assume individuals are released immediately before mating. If we assume density-dependence occurs because individuals use up resources as they grow to adulthood, then the release of lab-grown adults right before mating should not drain resources, and therefore do not affect our measure of density, in terms of density-dependent growth.

## Simulation Structure and Details
