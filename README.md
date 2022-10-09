# EADScalc



## Getting started

- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Automatically merge when pipeline succeeds](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

## Name
EADScalc

## Description
EADScalc is a python script that can be used to determine the degree of substitution (DS) of modified cellulose (nanomaterials) using organic elemental analysis results, inlcuding correction for water content of the typically hygroscopic cellulose samples.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

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
