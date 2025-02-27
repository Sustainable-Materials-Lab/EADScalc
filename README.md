# EADScalc

## Description
EADScalc is a python script that can be used to determine the degree of substitution (DS) of modified cellulose (nanomaterials) using organic elemental analysis results, inlcuding correction for water content of the typically hygroscopic cellulose samples.

## Installation
Use a python virtual environment and install via
`pip install EADScalc --extra-index-url https://gitlab.kuleuven.be/api/v4/projects/9529/packages/pypi/simple`

## Usage
To get help with the script, run:

`EADScalc -h`

The minimal amount of information you need to enter to run the calculation is:

`EADScalc <carbon_content> <hydrogen_content> <nitrogen_content> --nmod <number_of_modifications> --O1 <number_O_first_mod> --H1 <number_H_first_mod>`

However, the water content of the sample (humidity) and likely the sulfur content will be needed in most calculations.

For example (sulfate):

`EADScalc 41 6 0 --humidity 5 --sulfur 0.5 --nmod 1 --01 4 --S1 1 --H1 1`

When considering the empirical formula of the modification, consider all atoms that REPLACE the hydroxyl group that you are modifying. E.g. if you create cellulose acetate, the modification would be C2H3O2.

### Output
Results are output to the terminal and to the file named in the command line as a text file.

## Support
If there are issues, contact @u0092172

## Roadmap
No updates are planned for this project and it is in a dormant state.

## Contributing
Contributors are welcome.

## Authors and acknowledgment
Originally authored by @u0092172

## License
MIT License

## Project status
Dormant - I no longer use this in my day to day work, so I am not actively working on it, however, I am willing to fix problems.
