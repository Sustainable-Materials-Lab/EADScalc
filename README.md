# EADScalc

## Description
EADScalc is a python script that can be used to determine the degree of substitution (DS) of modified cellulose (nanomaterials) using organic elemental analysis results, inlcuding correction for water content of the typically hygroscopic cellulose samples.

## Installation
Requires lmfit and numoy. In windows, install the free anaconda python 3 distribution. Then install lmfit using pip:

`python -m pip install lmfit`

## Usage
To get help with the script, run:

`python EADScalc.py -h`

The minimal amount of information you need to enter to run the calculation is:

`python EADScalc.py <chain_ratio> <water_content> <carbon_content> <hydrogen_content> <nitrogen_content> <sulfur_content> <result_filename_without_ext> --nmod <number_of_modificationsf_(1)> --C1 <number_C_first_mod> --O1 <number_O_first_mod> --S1 <number_S_first_mod> --H1 <number_H_first_mod> (--C2...)`

For example:

`python EADScalc.py 1 3.0 41 6 0 0.5 'sulfate' --nmod 1 --C1 0 --01 4 --S1 1 --H1 1`

When considering the empirical formula of the modification, consider all atoms that REPLACE the hydroxyl group that you are modifying. E.g. if you create cellulose acetate, the modification would be C2H3O2.

## Support
If there are issues, contact Samuel Eyley

## Roadmap
No updates are planned for this project and it is in a dormant state.

## Contributing
Contributors are welcome.

## Authors and acknowledgment
Originally authored by Samuel Eyley.

## License
CC-SA 4.0

## Project status
Dormant - I no longer use this in my day to day work, so I am not actively working on it, however, I am willing to fix problems.
