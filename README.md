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

`python EADScalc.py <chain_ratio> <water_content> <carbon_content> <hydrogen_content> <nitrogen_content> <sulfur_content> <result_filename_without_ext> --nmod <number_of_modificationsf_(1)> --C1 0 --O1 4 --S1 1 --H1 1 (--C2...)`

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
