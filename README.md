# TOPSMODE

This application calculates the spectral moments descriptors based on the TOPological Sub-Structural MOlecular DEsign (TOPS-MODE) approach and the contribution of this along the molecule. The application is capable to calculate the spectral moments of the bond and atom matrix.
### Installation (choose one)
1. git clone https://github.com/sytoperez/TOPSMODE.git
2. git clone git@github.com:sytoperez/TOPSMODE.git
3. gh repo clone sytoperez/TOPSMODE
4. Download the .zip and unzip it in the supercomputing centers you are going to use
### Usage
Returns 2D TOPSMODE descriptors for compounds
```text TOPS -i input.sdf -o output (args) ```


Returns molecular contributions based on models coefficients/importance of descriptors
```text TOPS_Ctbr -i input.rdf -o output (args) ```

### Available arguments for TOPS.py (descriptors calculation)
-i input file (allowed formats: sdf) with standardized structures
-o output file name with calculated descriptors
-t descriptors type: ato for atomic calculation and bond for bond calculation
-d descriptors to build. Default all
-n Integer value. Number of descriptors of each type. Default: 15
-w Name of unique ID for compounds (sdf). gen - auto-generated names and titles - sdf titles will be used. If omitted for sdf molecule titles will be used or auto-generated names
-a Nname of field with activity values for compounds (sdf). none - activity values is ommited
-v Integer value. 0 - print no details. 1 and more - verbose output. Default: 0.

### Available arguments for TOPS_Ctbr.py (structure atom/bond contribution calculation)
-i input file (allowed formats: sdf) with standardized structures
-m File name (with full path) for contributions. Should contain at least these columns (named): "n", "variables" (separated with |), "coeff" (separated with |). If the coefficient are of linear model the first one must be the value of intercept
 -o output file with calculated contributions.
-w Name of unique ID for compounds (sdf). gen - auto-generated names and titles - sdf titles will be used
-t Descriptors type: ato for atomic calculation and bond for bond calculation.
-l Model technique used (lin-linear, non-nonlinear). Default linear
-d Only output .csv with data without structures in .png
-v Integer value. 0 - print no details. 1 and more - verbose output. Default: 0.

### CHANGELOG
**v1.1.0** (22/4/2024)
- Fix errors in atomic contribution calculation.

**v1.0.0** (7/11/2023)
Initial version
