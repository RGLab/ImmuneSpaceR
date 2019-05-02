# Contributing to `ImmuneSpaceR`

## Opening issues

Please briefly describe your problem and what output you expect and include a minimal reprex.

The goal of a reprex is to make it as easy as possible for me to recreate your problem so that I can fix it. If you've never heard of a reprex before, start by reading [this](https://github.com/jennybc/reprex#what-is-a-reprex).

---

Brief description of the problem

```r
# insert reprex here
```

---

### Issue labels

* bug: found something that's giving you a headache? use this flag!
* question: not sure how to use or what is going on? Ask away
* enhancement: new feature? refactor? optimization?
* infrastructure: tests, Travis, codecov, bioconductor, pkgdown
* help wanted: maybe you are stuck? need some help? use this flag to get attention



## General development guidelines

If you'd like to contribute changes to `ImmuneSpaceR`, we use [the GitHub flow](https://guides.github.com/introduction/flow/index.html) for proposing, submitting, reviewing, and accepting changes. If you haven't done this before, Hadley Wickham provides [a nice overview of git](http://r-pkgs.had.co.nz/git.html), as well as [best practices](http://r-pkgs.had.co.nz/git.html#pr-make) for submitting pull requests. We also recommend using [his style guide](http://adv-r.had.co.nz/Style.html) when writing code.


---


# Development Flow

## New Feature

`fb_newFeature` -> Code Review -> `dev` -> `master` -> BioConductor Submission

For example:

a. Run into a problem or come up with an idea for a new feature or optimization/enhancement
a. Create an issue to address the problem or propose the new feature
a. Create a new branch from `dev`
a. Make changes in the new branch
a. Build and install
a. Manually test out your changes
a. Make a commit
a. If you have more changes to make, repeat steps from d to g
a. Run `R CMD check` and `BiocCheck`
a. If you get an error or warning, repeat steps from d to i
a. Push your commit(s) to GitHub. This will kick off Travis build
a. Wait for Travis build to be done and repeat steps from d to k if there are errors or you'd like to make more changes
a. When the new branch is ready to merge, create a Pull Request. The maintainer will assign reviewers for your PR
a. Reviewers will go through the PR checklist and give you feedback. Repeat steps from d to l accordingly
a. Once the reviewers approve your PR, the maintainer will merge your branch to `dev` branch

## Hot Fix

`dev` -> `master` -> BioConductor Submission


---


# Package Structure (WIP)

## ImmuneSpaceR.R

* Package description
* Declare global variables for `data.table` scoping

## CreateConnection.R

* `CreateConnection()`: constructor for `ISCon` class

## ISCon.R

* `ISCon` definition and documentation
  * `self$study`
  * `self$config`
  * `self$availableDatasets`
  * `self$cache`
  * `private$.constants`
* `self$initialize()`
* `.check_internet()`
* `.get_host()`
* `.check_portal()`
* `.get_url_base()`
* `.get_user_email()`
* `.get_url_path()`
* `.set_curl_options()`
* `.check_credential()`

## ISCon-cytometry.R

* `self$listWorkspaces()`
* `self$listGatingSets()`
* `self$summarizeCyto()`
* `self$summarizeGatingSet()`
* `self$loadGatingSet()`
* `private$.openWorkspace()`
* `private$.downloadCytoData()`
* `private$.mergePD()`
* `.isRstudioDocker()`
* `.assertRstudio()`
* `.buildGSPath()`

## ISCon-dataset.R

* `self$listDatasets()`
* `self$getDataset()`
* `private()$.setAvailableDatasets()`
* `private$.getColnames()`
* `private$.checkFilter()`
* `.extractNames()`
* `.fixColFilter()`
* `.transformData()`

## ISCon-geneExpression.R

* `self$listGEMatrices()`
* `self$listGEAnalysis()`
* `self$getGEMatrix()`
* `self$getGEAnalysis()`
* `self$getGEInputs()`
* `self$getGEFiles()`
* `self$downloadGEFiles()`
* `self$addTreatment()`
* `self$mapSampleNames()`
* `private$.downloadMatrix()`
* `private$.getGEFeatures()`
* `private$.constructExpressionSet()`
* `private$.mungeFeatureId()`
* `.setCacheName()`
* `.combineEMs()`

## ISCon-participantGroup.R

* `self$listParticipantGroups()`
* `self$getParticipantData()`
* `private.assertAllStudyConnection()`

## ISCon-plot.R

* `self$plot()`
* `.plot()`
* `.standardize_time()`
* `.qpHeatmap2()`
* `.qpBoxplotViolin()`
* `.qpLineplot()`
* `.getDataToPlot()`
* `.heatmapAnnotations()`
* `.format_lab()`

## ISCon-utils.R

* `self$print()`
* `self$clearCache()`
* `private$.checkStudy()`
* `private$.munge()`
* `private$.isProject()`
* `private$.isRunningLocally()`
* `private$.localStudyPath()`
* `private$.listISFiles()`
* `.getLKtbl()`
* `.mixedsort()`

## netrc.R

* `interactive_netrc()`
* `write_netrc()`
* `check_netrc()`
* `.get_env_netrc()`
* `.get_env_url()`
* `.get_path()`
* `.check_con()`

## template.R

* `template_IS()`

## theme.R

* `theme_IS()`
* `ISpalette()`
* `.override_scale()`

## utils.R

* `loadConnection()`
* `saveConnection()`

## zzz.R

* `.onAttach()`


---


# `ISCon-*.R` Structure (WIP)

```R
#' @include ISCon.R
NULL



# PUBLIC -----------------------------------------------------------------------
# (fields and methods that will be public)
# (they should be documented in ISCon.R file)

# (comment what/how this method is doing)
ISCon$set(
  which = "public",
  name = "aPublicMethod",
  value = function() {
  
  }
)


# (comment about the method)
ISCon$set(
  which = "public",
  name = "anotherPublicMethod",
  value = function() {

  }
)



# PRIVATE -----------------------------------------------------------------------

# (comment about this method)
ISCon$set(
  which = "private",
  name = ".aPrivateMethod",
  value = function() {
  
  }
)


# (comment)
ISCon$set(
  which = "private",
  name = ".anotherPrivateMethod",
  value = function() {
  
  }
)



# HELPERS -----------------------------------------------------------------------

# (comment)
.aHelper <- function(x) {

}


# (comment)
.anotherHelper <- function() {

}
```

* `ISCon-*.R` files are divided into four sections: PURPOSE, PUBLIC, PRIVATE, HELPERS
* PURPOSE section
    * Describes the purpose of this .R file and functions/methods as group
* PUBLIC section
    * For functions that require fields or methods in self
    * These methods should be documented in `ISCon.R` file
* PRIVATE section
    * Like public methods, this section are for functions that require fields or methods in self, not accessible to the user
    * The name of a PRIVATE method should start with a dot (`.`)
* HELPER section
    * Contains helper functions exclusively used in public/private methods in this .R file
    * The name of a HELPER function should start with a dot (`.`)
* Spacing
    * Two new lines between methods/functions in the same section
    * Three new lines between the sections


---


# Package Development Guide

## Functional Code Guide

* Use `data.table::data.table()` instead of `data.frame()`
* Use `data.table::fread()` instead of `read.csv()` or `read.table()`
* Use `TRUE` or `FALSE` instead of `T` or `F`
* Avoid `sapply()`, use `vapply()` instead
* Use `message()` and `warning()` to communicate with the user in your functions
* Use `self$selectRows()` or `self$executeSql()` instead of `labkey.selectRows()` or `labkey.executeSql()` (WIP)
* Use `fieldname` for `colNameOpt` argument if you decide to use `labkey.selectRows()` or `labkey.executeSql()`
* Use [the `profvis` package](https://rstudio.github.io/profvis/) to test the performance of your code

## Style Guide

* We are mainly following [the tidyverse style guide](https://style.tidyverse.org/) except for the object naming convention
* Run [`styler::style_pkg()`](https://styler.r-lib.org/reference/style_pkg.html) for convenience
* The maximum line length is **80 characters**
* Use `<-` not `=` for assignment
* When indenting your code, use **two spaces**. Never use tabs or mix tabs and spaces
* Use common sense and **BE CONSISTENT**

## Naming Things

* Use `lowerCamelCase` for variables, functions, methods, arguments
* Use `snake_case` for column names in tables
* Use specific words or look up shorthands in the glossary table
* Avoid generic names
* Use verbs for functions and methods

## Commenting Things

* Comment every functions and methods
* Comments should explain the why, not the what

## Notes

* Even a smallest change can break the package and the portal, so be mindful about the changes you make

## Glossary (WIP)

* WIP: Work in Progress


---


# Code Review for Pull Requests

* Reviewing code is a lot of work for both code owner and reviewer. The purpose of the following checklists is to set expectations for the code review and reduce time spent on reading the code.
* Code owners, please respect reviewer's time. Reviewers, please be respectful and kind to the code owners in your reviews

## Checklist for Reviewer

* [ ]  Build it yourself and test the new features or changes
* [ ]  Is it functional, and does it meet the requirements?
* [ ]  Is it readable?
  * Easy to understand
* [ ]  Is it following the coding standards/guidelines?
* [ ]  Functional coding guide
* [ ]  Style guide (No changes in `styler::style_pkg()`)
* [ ]  Is it well documented (clear and sufficient)?
  * Comments make sense
  * Updated or new vignette or manual
* [ ]  Does it have new/updated unit tests?
  * Use the `covr` package to examine the extent of test coverage
* [ ]  Does it have new dependencies? Is it justifiable?
* [ ]  Is there code duplication that should be reduced?
* [ ]  Are there user interface improvements that could be made?
* [ ]  Are there performance improvements that could be made?
* [ ]  Is it passing `R CMD check`?
* [ ]  Is it passing `BiocCheck::BiocCheck()`?


---


# Maintainer Guide

## Merge to `master` branch Checklist

* [ ]  Is the version bumped?
* [ ]  Is the documentation up to date?
* [ ]  Is the pkgdown site up to date?
* [ ]  Is it passing `BiocCheck::BiocCheck()`?
* [ ]  Is it passing against the production?

## BioConductor Submission Checklist

* [ ]  Is the version valid?
* [ ]  Is it passing `BiocCheck`?
* [ ]  Does it need to be merged to the release branch of BiocConductor?


---


# FAQ (WIP)

### What is the point of this?

* The point of having these guidelines is "to have **a common vocabulary** of coding so people can **concentrate on what you are saying**, rather than on how you are saying it" ([Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml)).

### How does `ImmuneSpaceR` work?

* Well

### How do I setup my machine for `ImmuneSpaceR` development

* See [Setup local for travis-like testing](https://www.notion.so/0d597ef4-4b78-49bd-a49f-a4d81c1b0f98)


---


# Further Readings

## R related

* [R Packages](http://r-pkgs.had.co.nz/style.html) by Hadley Wickham
* [The tidyverse style guide](https://style.tidyverse.org/)
* [Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml)
* [rOpenSci Guide](https://ropensci.github.io/dev_guide/)
* [BioConductor Package Guidelines](https://www.bioconductor.org/developers/package-guidelines/)
* [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)

## General coding best practices

* [The Art of Readable Code](https://www.oreilly.com/library/view/the-art-of/9781449318482/) by Trevor Foucher and Dustin Boswell
* [Clean Code](https://www.oreilly.com/library/view/clean-code/9780136083238/) by Robert C. Martin
