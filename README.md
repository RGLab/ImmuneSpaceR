ImmuneSpaceR
============

A thin wrapper around Rlabkey to access the ImmuneSpace database from R
This package simplifies access to the HIPC ImmuneSpace database, for R programmers.

It takes advantage of the standardization of the ImmuneSpace database to hide all the `Rlabkey` specific code away from the user.
Study-specific datasets can be accessed via an object-oriented paradigm.

## Instalation
The package requires a new version of `Rlabkey` available [here](https://github.com/RGLab/Rlabkey).

The package can be downloaded here and installed like any other R packages or installed directly from github using devtools.

    library(devtools)
    install_github("RGLab/Rlabkey")
    install_github("RGLab/ImmuneSpaceR")


The database is accessed with the user's credentials. A `.netrc` file storing login and password information is required.

Create a file named .netrc in the home directory of the computer running R.
The file should contain 3 rows with key/value separated by tabs, as follow:  

    machine	www.immunespace.org  
    login   <user.login>  
    password        <user.password>  

See `man netrc` for the official documentation.

## Usage
The general idea is that the user creates an instance of an `ImmuneSpaceConnection` class. 
The instance configures itself to connect to a specific study, and datasets and gene expression matrices can be retrieved by name.

For example:

```
study<-CreateConnection("SDY269)
```

will create an instance of study 269.  The user needs credentials stored in a `.netrc` file to access the database. 
Datasets can be listed by:

```
study$listDatasets()
```

which will print names of available datasets and gene expression matrices.

Gene expression matrices or data sets can be retreived by:

```
study$getGEMarix("name")
#or
study$getDataset("name")
```

The study object *caches* data, so once it is retrieved, the next time you access it, it will use the local cached copy. 

The package uses a simple S5 reference class to represent the connection to a study and get around some of R's copy-on-change behaviour.

## Quick plots

The `quick_plot` function uses ggplot2's qplot function to generate quick plots of data sets, leveraging the standardized data set tables. 

```
quick_plot(study$getDataset("hai"))
```


