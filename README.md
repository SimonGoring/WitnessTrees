# WitnessTrees
This is the repository for Goring *et al*., but also contains much of the core processing and analytic code for the Upper Midwestern Public Land Survey data.

The code supports the publication of the Goring *et al*. and supports secondary processing of other datasets, including the gridding and aggregation of other partially complete datasets in southern Michigan, Indiana and Illinois.

Key data includes:
 * data/input
   * [NaturalEarth](https://github.com/SimonGoring/WitnessTrees/tree/master/data/input/NaturalEarth) - This data is used for mapping, but is used in mapping for the paper.  Data is obtained from http://naturalearthdata.org
   * [rasters](https://github.com/SimonGoring/WitnessTrees/tree/master/data/input/rasters) - includes all FIA base data (from Sydne Record) gridded to an 8x8km Albers grid used in analysis, as well as a base unit map.
   * [relation_tables](https://github.com/SimonGoring/WitnessTrees/tree/master/data/input/relation_tables) - These are key tables for this analysis, they include tables to convert FIA data to PalEON taxa, conversions from original settlement taxon assignments to PalEON taxa, and PalEON assignments to PFTs.  This also includes the correction table, to help build better quality assessments of stem density, as presented in Goring et al.
   * [shapes](https://github.com/SimonGoring/WitnessTrees/tree/master/data/input/shapes) - shapefiles for Canada and the USA for use in mapping.
  * data/output
    * [aggregated_midwest](https://github.com/SimonGoring/WitnessTrees/tree/master/data/output/aggregated_midwest) - This should be the files output by running either the southern Michigan data, or the Goring *et al*. paper code.
    * [gridded](https://github.com/SimonGoring/WitnessTrees/tree/master/data/output/gridded) - This data is used by the Paciorek et al. composition model, it is the data generated from gridding each of the constituent domains.
    * [southern_MI](https://github.com/SimonGoring/WitnessTrees/tree/master/data/output/southern_MI) - 
    * [tests](https://github.com/SimonGoring/WitnessTrees/tree/master/data/output/tests) - clean_bind_test in this folder checks the species conversions to make sure they're capturing the proper taxa in the conversion table.
    * [wiki_outputs](https://github.com/SimonGoring/WitnessTrees/tree/master/data/output/wiki_outputs) - outputs for the PalEON wiki, includes CSVs in Albers projection for density, basal area biomass, point numbers and plot numbers.
  * R
    * [composition_raster](https://github.com/SimonGoring/WitnessTrees/tree/master/R/composition_raster) - The code used to generate the composition raster used in Paciorek *et al*. from individual datasets in southern Michigan, Illinois, Indiana and the upper Midwest.
    * [paper](https://github.com/SimonGoring/WitnessTrees/tree/master/R/paper) - The full code for Goring *et al*.
    * [process_raw](https://github.com/SimonGoring/WitnessTrees/tree/master/R/process_raw) - The code to join Michigan, Minnesota and Wisconsin, and to fill the hole in the upper Peninsula of Michigan using the southern Michigan dataset.
