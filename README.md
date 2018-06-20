# MTseeker

[![Build Status](https://travis-ci.org/trichelab/MTseeker/MTseeker.png?branch=master)](https://travis-ci.org/trichelab/MTseeker/MTseeker)  [![codecov](https://codecov.io/gh/trichelab/MTseeker/MTseeker/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/MTseeker/MTseeker)

## How to install

devtools::install_github("trichelab/MTseeker", username="yourGithubUsername", auth_token="yourGitHubUserTokenFromGithub.com/settings/tokens") 

## ToDo items:

### 1. Travis

Now you can go to [Travis](https://travis-ci.org/profile/trichelab/MTseeker) and turn on continuous integration for your new package. You may need to click the "Sync account" button to get your new package to show up in the list.

If you have a codecov.io account, running your tests on Travis will trigger the code coverage job. No additional configuration is necessary

### 2. Appveyor

Go to [Appveyor's new project page](https://ci.appveyor.com/projects/new) and select your new repository from the list. Then you can go to the [badges](https://ci.appveyor.com/project/trichelab/MTseeker/MTseeker/settings/badges) page, copy the markdown code it provides, and paste it up with the other badges above. (Their badge API has a random token in it, so `skeletor` can't include it in the template for you.)

## Installing

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/hadley/devtools) package:

```
# install.packages("devtools")
devtools::install_github("trichelab/MTseeker/MTseeker",
                         username="yourGithubUsername", 
                         auth_token="yourTokenFromGithub.com/settings/tokens",                           build_vignettes=TRUE)
```

Note that the access token and username are required as the repo is private.

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
