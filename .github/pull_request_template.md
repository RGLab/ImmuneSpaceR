<!-- DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->


## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4" or similar - or if just relates to an issue make sure to mention it like "#4" -->


## Example
<!--- if introducing a new feature or changing behavior of existing methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing grammar, please include new tests for your change -->


## Checklist for Reviewer

<!---  Reviewing code is a lot of work for both code owner and reviewer. The purpose of the following checklists is to set expectations for the code review and reduce time spent on reading the code --->
<!--- Code owners, please respect reviewer's time. Reviewers, please be respectful and kind to the code owners in your reviews --->

This Pull Request should follow [ImmuneSpaceR Package Development Guide](
https://github.com/RGLab/ImmuneSpaceR/blob/main/CONTRIBUTING.md#package-development-guide)

* [ ]  Can you build it yourself and test the new features or changes?
* [ ]  Is it functional, and does it meet the requirements?
* [ ]  Is it readable?
  * Easy to understand
* [ ]  Is it following the functional coding guide?
* [ ]  Is it following the style guide (No changes in `styler::style_pkg()`)?
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
