# TOPSMODE

This application calculates the spectral moments descriptors based on the TOPological Sub-Structural MOlecular DEsign (TOPS-MODE) approach and the contribution of this along the molecule. The application is capable to calculate the spectral moments of the bond and atom matrix.
### Installation (choose one)
1. git clone https://github.com/sytoperez/TOPSMODE.git
2. git clone git@github.com:sytoperez/TOPSMODE.git
3. gh repo clone sytoperez/TOPSMODE
4. Download the .zip and unzip it in the supercomputing centers you are going to use
### Available arguments
-i input file (allowed formats: sdf) with standardized structures
-o output file name with calculated descriptors
-t descriptors type: ato for atomic calculation and bond for bond calculation
-d descriptors to build. Default all
-n Integer value. Number of descriptors of each type. Default: 15
-w Name of unique ID for compounds (sdf). gen - auto-generated names and titles - sdf titles will be used. If omitted for sdf molecule titles will be used or auto-generated names
-a Nname of field with activity values for compounds (sdf). none - activity values is ommited
-v Integer value. 0 - print no details. 1 and more - verbose output. Default: 0.

### CHANGELOG
**v1.1.0** (22/4/2024)
- Fix errors in atomic contribution calculation.
**v1.0.0** (7/11/2023)
Initial version
