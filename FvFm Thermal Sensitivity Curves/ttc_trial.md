Thermal Sensitivity Curves
================
German Vargas G.

Fv/Fm Thermal Sensitivity Curves:
-----------------------

These are measurements done in order to obtain the temperature at which Fv/Fm drops 50% (T50, $^{\circ}C$), (Tcrit, $^{\circ}C$) and the maximum Fv/Fm.

This file is a tutorial on how to analyze thermal sensitivity curves data. For more details on how to perform pressure volume curves and theoretical framework see:

-list of citations

Load the functions
------------------

``` r
source(file = "functions/T_tol_func_V1.R")#directory to the scrpt files containing the function
```

Load the data
-------------

``` r
trial <- read.csv(file = "data/Thermal_tolerance_sample.csv")

head(trial,15)
```

Fiting the TT curves
--------------------

### Just for one individual of a given species

``` r
Thermal_output <- Thermal_tol(fv_fm = fv.fm_aft,
                              temp = tmp,
                              species = code,
                              sp_levels = unique("SIDCAP"),
                              id = ind,
                              data = trial,
                              print_graph = TRUE,
                              graph_dir = "graphs/thermal_sp_id/")
```
### Leaf area provided

In this case the output will contain data for each individual of each species in the data

``` r

```
