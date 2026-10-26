# EADScalc

## Description

EADScalc calculates the degree of substitution (DS) of modified cellulose from
organic elemental analysis results. It can account for water in hygroscopic
cellulose samples.

## Installation

Install the package in your Python environment:

```console
pip install EADScalc
```

## Usage

Show the command-line options with:

```console
EADScalc --help
```

The required positional values are detected carbon, hydrogen, and nitrogen
percentages. You must also set `--nmod` and provide the hydrogen and oxygen
counts for the first modification (`--H1` and `--O1`):

```console
EADScalc <carbon_percent> <hydrogen_percent> <nitrogen_percent> --nmod <number_of_modifications> --H1 <modification_hydrogen_count> --O1 <modification_oxygen_count>
```

For example, this estimates one sulfate modification from elemental analysis
percentages. The detected sulfur option is `--sulf`; modification options are
case-sensitive (`--O1`, not `--01`):

```console
EADScalc 41 6 0 --humidity 5 --sulf 0.5 --nmod 1 --O1 4 --S1 1 --H1 1
```

The elemental-analysis values and `--humidity` are percentages. Modification
counts describe the atoms in the group replacing a cellulose hydroxyl. For
example, use `C2H3O2` for an acetate group (`--C1 2 --H1 3 --O1 2`).

The optimizer uses `scipy.optimize.least_squares`. `--method` accepts `trf`
(the default) or `dogbox`; both support the DS bounds. Set `--G1`, `--G1min`,
and `--G1max` to adjust the first modification's initial guess and bounds. The
corresponding `G2` and `G3` options apply when fitting multiple modifications.

### Output

Results are printed to the terminal and written to `result.txt` in the current
working directory.

## Support
If there are issues, report via https://github.com/Sustainable-Materials-Lab/EADScalc/issues

## Roadmap
No updates are planned for this project and it is in a dormant state.

## Recent Updates
**2026**: Replaced lmfit with SciPy's bounded `least_squares` optimizer. The
available methods are `trf` and `dogbox`; the default is `trf`.

**2025**: Major code refactoring completed
- Reduced code duplication by consolidating repetitive patterns into data-driven structures
- Implemented decorator factories for dynamic Click option generation
- Replaced 39+ manual option decorators with 2 reusable functions
- Consolidated 13+ atomic mass variables into a single dictionary
- Refactored residual calculation from 13 if-blocks to a single loop
- Improved code maintainability through data-driven structures and reusable decorators

## Contributing
Contributors are welcome.

## Authors and acknowledgment
Originally authored by Samuel Eyley. 

Code refactored with the help of GitHub Copilot.

## License
MIT License

## Project status
Dormant - I no longer use this in my day to day work, so I am not actively working on it, however, I am willing to fix problems.
