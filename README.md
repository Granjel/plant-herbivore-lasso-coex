Computational analyses for the article "Simple interaction structure governs niche and fitness differences in a multitrophic grassland community."

## Data

### `data/raw/empirical-dataset.txt`

This dataset contains plant community composition and cover data from a grassland experiment manipulating herbivore assemblages and plant communities. Each row represents a sampling unit within a given block, treatment, and time point. Metadata columns describe the experimental design (time, date, block, treatment, datapoint, focal species, and total cover), while the majority of columns correspond to individual plant species recorded as percent cover. Additional columns summarise functional groups (legumes, grasses, other) and the composition of the herbivore community (grasshopper species abundances or presence).

The data are structured to support the estimation of plant–plant and plant–herbivore interactions in a multitrophic context. Species-level cover values form the community matrix used for modelling, while herbivore variables encode the experimentally controlled treatments that allow inference of both direct effects and higher-order interaction modifications. The dataset includes repeated observations across treatments and time, enabling analyses of community dynamics and coexistence mechanisms.

### `data/raw/full-plant-community-cover.txt`

This dataset contains species-level plant cover observations in long format. Each row represents a single observation of a species (“Spcode”) and its corresponding percent cover value (“Cover”) recorded within a sampling unit. It is used to produce Table S1.

## Code

The scripts in `code` represent the main analytical pipeline and can be run using the `code/00-pipeline.R` script, which is self-explanatory. The scripts in `suppl-mat/suppl-code` can be run in any order and reproduce the results shown in the supplementary materials.
