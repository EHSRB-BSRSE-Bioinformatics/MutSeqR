# Find the appropriate BS genome for the specified organism and genome.

This function will browse available BSgenomes, indicating which one
should be installed for the specified organism and genome assembly
version. If you cannot specify both organism and genome, the function
can return a list of available genomes for a specified species.

## Usage

``` r
find_BS_genome(organism, genome, masked = FALSE)
```

## Arguments

- organism:

  the name of the organism for which to install the reference genome.
  This can be the scientific name or a common name. For example Homo
  Sapiens, H. sapiens, or human

- genome:

  The reference genome assembly version. Ex. hg18, mm10, rn6.

- masked:

  Logical value. Whether to search for the 'masked' BSgenome. Default is
  FALSE.

## Value

a BSgenome object

## Examples

``` r
# Find the reference genome for Mouse, mm10 assembly:
mouse_mm10 <- find_BS_genome("mouse", "mm10")
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cran.rstudio.com
#> Selected reference genome: BSgenome.Mmusculus.UCSC.mm10
#> Reference genome is already installed.
#> Once installed, supply: BSgenome.Mmusculus.UCSC.mm10 as the BS_genome parameter in import_mut/vcf_data()
# Find all possible mouse BS genomesL
mouse_all <- find_BS_genome("mouse")
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cran.rstudio.com
#> Possible BS genomes for organism = 'mouse', masked = FALSE: BSgenome.Mmusculus.UCSC.mm10BSgenome.Mmusculus.UCSC.mm39BSgenome.Mmusculus.UCSC.mm8BSgenome.Mmusculus.UCSC.mm9. Please install one of the possible BS genomes using BiocManager::install('pkgname') and provide the pkgname to import_mut/vcf_data().
```
