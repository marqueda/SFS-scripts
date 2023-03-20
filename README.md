# SFS-scripts

A collection of tools useful for preparing or manipulating site-frequency spectrum (SFS) files in the format as used by fastsimcoal2 (see http://cmpg.unibe.ch/software/fastsimcoal2/).

### sampleKgenotypesPerPop.py and choosekgenotypes.py

These scripts are useful for the preparation of SFS files from genotypes in a VCF file, for which no missing data is allowed. They randomly draw k non-missing genotypes per population from a VCF file. Sites with < k non-missing genotypes per population are discarded.

sampleKgenotypesPerPop.py: Uses a population file with format IndividualTABPopulation and allows to specify different k's for different popoulations.

choosekgenotypes.py: Assumes that the individual names in the VCF file header contain the population code in the format "individual.popcode[.more]" separated by a period from the individual name.

### vcf2sfs.py

This script converts VCF files into site frequency spectra (SFS) files of the respective dimensions in fastsimcoal2 format (1D, 2D or multi-D, depending on the number of populations in the VCF file) by simply counting alleles from genotypes for sites without missing data. It can also compute SFS for non-overlapping windows in the genome (defined by distance along the chromsome or number of sequenced sites in the VCF file), compute bootstrap replicates of the SFS (by resampling single sites) or block-bootstrap repliates of the SFS (by resampling non-overlapping windows in the genome).
IMPORTANT: Sites with missing data are discarded. Use the script choosekgenotypes.py (see above) first to subsample genotypes per population to get rid of missing data.

### SFStools.R

This script can:
- convert multi-dimensional SFS to marginal 2D-SFSs
- convert 2D-SFS to the marginal 1D-SFSs
- visualize 2D-SFSs and compare observed and simulated 2D-SFSs present in the same folder
- visualize 1D-SFSs and compare observed and simulated 1D-SFSs present in the same folder
SFS must be in the correct format (see above) and it currenlty supports only SFS files with a single SFS stored.

### fold1DSFS.R, fold2DSFS.R, foldSFS.py

These scripts fold unfolded 1D, 2D SFS or SFS of any dimension. Currently supports only SFS files with a single SFS stored.
