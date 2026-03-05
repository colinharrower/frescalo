frescalo
================
Colin Harrower
2025-05-20

The frescalo package contains an R implementation of Frescalo (FREquency
SCAling LOcal) modeling methodology invented by Dr Mark O. Hill to aid
with analysis of biological record data. The methodology attempts to
correct for spatial and temporal variation in recorder effort that is
common in biological recording data (i.e. species occurrence data).

## Installation

You can install the development version of frescalo from
[GitHub](https://github.com/colinharrower/frescalo) with:

``` r
# install.packages("pak")
pak::pak("colinharrower/frescalo")
```

## Example

The package contains simulated data which can be used as an example and
or to test the package. The simulated data comrpises two data objects
`s` containing occurrence data and `d` containing neighborhood
weightings. The occurrence data in `s` has already been summarized to
the format required by frescalo, specifically a `data.frame` containing
all unique combinations of time period (`time`), sites or locations
(`location`) and species ids or names (`species`). The neighbourhood
weights data specifies for each focal location (`location1`) which other
sites form its “neighbourhood” (`location2`) and the relevant weighting
score between the focal site and each respective neighbourhood site
(`w`).

``` r
# Load the package
library(frescalo)

# View the first 10 rows of the test dataset s (occurrence data)
head(s,10)
```

    ##    time location    species
    ## 1     1     SO34 Species 11
    ## 2     1     SO37 Species 11
    ## 3     2     SO43 Species 11
    ## 4     1     SO34 Species 11
    ## 5     1     SN96 Species 11
    ## 6     2     SO43 Species 11
    ## 7     2     SN96 Species 11
    ## 8     1     SO05 Species 11
    ## 9     1     SO05 Species 11
    ## 10    2     SN62 Species 11

``` r
# View the first 10 rows of the weights datasets (d)
head(d,10)
```

    ##    location1 location2      w
    ## 1       TV69      SU40 0.0004
    ## 2       TV69      SU41 0.0000
    ## 3       TV69      SU42 0.0000
    ## 4       TV69      SU50 0.0137
    ## 5       TV69      SU60 0.0639
    ## 6       TV69      SU65 0.0005
    ## 7       TV69      SU70 0.0410
    ## 8       TV69      SU76 0.0001
    ## 9       TV69      SU77 0.0006
    ## 10      TV69      SU83 0.0088

In this example dataset the species are identified using simplistic
species IDs and the locations are identified using Ordnance Survey of
Great Britain 10 km grid references, however the site and species
identifiers can be any appropriate names or coded identifiers. The time
periods can be identified as simple coded or numerical identifiers, as
in the example data `d`, or using appropriate numerical values
(e.g. mid-point of the time periods).

The main wrapper function `frescalo()` is used to fit apply the frescalo
model to the data, returning a list comprised of 4`data.frames` `locs`
with the location metrics, `freq` with original and re-scaled species
frequencies, `trend` with tFactors for each species, `site_time` with
the site by time recording effort metrics used for the adjustments. The
outputs are equivalent to the `samples.txt`, `frequencies.out` and
`trend.out` output files produced by Mark Hill’s original fortran
frescalo program.

``` r
# Use test dataset (s) and weights data (d) included with the package
out_fres = frescalo(s,d,in_parallel = FALSE,filter_wts = TRUE)

# View outputs
## Location
head(out_fres[["locs"]])
```

    ##   location nSpecies    phi_in    alpha phi_out spnum_in spnum_out iter
    ## 1     TR36        2 0.1367205 20.08779    0.74 1.952163  16.64654    7
    ## 2     TR34        1 0.1190483 34.61611    0.74 1.100281  13.62739    5
    ## 3     TR26        3 0.1283845 30.75372    0.74 1.282218  13.81837    6
    ## 4     TR16       10 0.1869391 10.88084    0.74 2.980106  15.73302    7
    ## 5     TR15        1 0.1855399 10.81524    0.74 2.970722  15.77083    7
    ## 6     TR14        2 0.1410191 14.16595    0.74 2.017961  13.92664    6

``` r
## Frequency
head(out_fres[["freq"]])
```

    ##   location    species pres       freq    freq_1 rank     rank_1 benchmark
    ## 1     TR36  Species 1    1 0.27165159 0.9982833    1 0.06007253         1
    ## 2     TR36 Species 13    1 0.23737656 0.9956761    2 0.12014507         1
    ## 3     TR36  Species 2    0 0.21097121 0.9914331    3 0.18021760         1
    ## 4     TR36 Species 32    0 0.13732783 0.9485621    4 0.24029013         1
    ## 5     TR36  Species 8    0 0.10259340 0.8863268    5 0.30036266         0
    ## 6     TR36  Species 5    0 0.09602601 0.8683959    6 0.36043520         0

``` r
## Trend
head(out_fres[["trend"]])
```

    ##      species time   tFactor      StDev
    ## 1  Species 1    1 0.6291580 0.05149657
    ## 2  Species 1    2 0.4358530 0.04694728
    ## 3 Species 10    1 0.2470118 0.04133537
    ## 4 Species 10    2 0.2549233 0.04915949
    ## 5 Species 11    1 0.3490210 0.14781598
    ## 6 Species 11    2 0.2216005 0.13160626

``` r
## Site Time
head(out_fres[["site_time"]])
```

    ##   location time s_it     w
    ## 1     NC02    1    1 1.000
    ## 2     NC02    2    0 0.005
    ## 3     NC11    1    1 1.000
    ## 4     NC11    2    0 0.005
    ## 5     NC13    1    1 1.000
    ## 6     NC13    2    0 0.005
