# ImmuneSpaceR 1.7.3

## Major changes

* Switched away from R's standard reference classes to `R6` classes (#57).
* `httr` is now the new dependency replacing `RCurl` (#54).

### `ImmuneSpaceConnection` class

* The connection class is now based on the `R6` class system (#57).
* `getGEMatrix` method can now access the raw matrix for a cohort using the 'outputType = raw' argument and return probe-level (microarray) or gene-counts (RNAseq), non-normalized or summarized data (#63).
* `getGEmatrix` method now has a boolean 'verbose' argument that when set to TRUE will print information on how the matrix was constructed, such as normalization procedure and annotation version (#66).
* Updated `print` method to show a prettier summary of the object (#62).
* `getDataset` method will now return an empty data frame for invalid dataset names (#68).
* Restructured fields and methods into public and private (#59, #60).
* Some fields and methods were renamed to be consistent with the existing naming convention (verb + camelCase), and others are documented for users (#59, #60).
* These are currently available fields and methods. Please note the name changes and deprecated fields and methods in the new version. Check the class documentation for detailed description of fields and methods.

#### Fields

* `study`
* `availableDatasets` (renamed from `available_datasets`)
* `cache` (renamed from `data_cache`)
* `config`

#### Methods

* `print` (renamed from `show`)
* `listDatasets`
* `listGEMatrices` (renamed from `GeneExpressionMatrices`)
* `listGEAnalysis`
* `listParticipantGroups`
* `getDataset`
* `getGEMatrix`
* `getGEAnalysis`
* `getGEFiles`
* `getGEInputs` (renamed from `GeneExpressionInputs`)
* `getParticipantData`
* `addTreatment`
* `mapSampleNames` (renamed from `EMNames`)
* `plot` (renamed from `quick_plot`)
* `clearCache`  (renamed from `clear_cache`)

## Minor bug fixes and improvements

* Reorganized the file structure in `/R` (#60).
* Fixed `.onLoad()` note on R CMD CHECK (#61).
* Removed setting of RCurlOptions cainfo and CA certificates as it is no longer needed for `httr` (#61).
* Modified user agent to capture more information about user's R environment (i.e., "R/3.4.3 (Linux x86_64) Rlabkey/2.2.0 ImmuneSpaceR/1.7.2") (#61).
* Built a pkgdown site (#67).
* Added an argument (`mc.cores`) to `checkRawFiles` method as a temporary fix for `httr`'s bug (not working well with `mclapply`) (#64).
* Updated `checkStudyCompliance` (#65).
