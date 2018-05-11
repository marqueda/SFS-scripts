# SFS-scripts

A collection of tools useful for preparing or manipulating site-frequency spectrum (SFS) files in the format as used by fastsimcoal2 (see http://cmpg.unibe.ch/software/fastsimcoal2/).

### sampleKgenotypesPerPop.py and choosekgenotypes.py

These scripts are useful for the preparation of SFS files from genotypes in a VCF file, for which no missing data is allowed. They randomly draw k non-missing genotypes per population from a VCF file. Sites with < k non-missing genotypes per population are discarded.
sampleKgenotypesPerPop.py: Uses a population file with format <IndividualTABPopulation> and allows to specify different k's for different popoulations.
choosekgenotypes.py: Assumes that the individual names in the VCF file header contain the population code in the format "individual.popcode[.more]" separated by a period from the individual name.

### SFStools.R

This script can:
- convert multi-dimensional SFS to marginal 2D-SFSs
- convert 2D-SFS to the marginal 1D-SFSs
- visualize 2D-SFSs and compare observed and simulated 2D-SFSs present in the same folder
- visualize 1D-SFSs and compare observed and simulated 1D-SFSs present in the same folder
SFS must be in the correct format (see above) and it currenlty supports only SFS files with a single SFS stored.

### fold1DSFS.R & fold2DSFS.R

These scripts fold unfolded 1D or 2D SFS. Currently supports only SFS files with a single SFS stored.
