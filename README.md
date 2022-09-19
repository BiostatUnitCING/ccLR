# Case-control likelihood ratio (LR) #
**ccLR** (Case-Control Likelihood Ratio)

This repository contains freely accessible scripts for the calculation of case-control likelihood ratios (LRs). This method compares the likelihood of the distribution of the variant of interest between cases and controls, under the hypothesis that the variant has similar age-specific relative risks (RRs) as the “average” *BRCA1/2* truncating variant, compared to the hypothesis that it is not associated with increased BC risk. This method can be also applied for the calculation of case-control LRs in variants residing in other genes, other than *BRCA1* and *BRCA2*.

## Installation ##
All of the components behind this analysis are freely available.

The following programms and packages are required:
* [R](https://www.r-project.org/) programming language

* R package [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
   ```
   install.packages("devtools")
   devtools::install_github("tidyverse/tidyverse")
   ```
* R package [R.utils](https://cran.r-project.org/web/packages/R.utils/index.html)
   ```
   install.packages("R.utils")
   ```

**Once all the prerequisite packages are installed, you can manually install the script from this GitHub repository.**

## Inputs ##

**To run the analysis you will use the script `ccLR.R`, for which you need to provide:**

1. The gene of interest

    According to the gene of interest, the user must specify the parameter `--gene` as:
    
   `--gene BRCA1` or `--gene BRCA2` or `--gene Custom`
   
______________________________________________________________________________________________________________________________________________________________________________
   
2. A phenotype file at a tab-delimited format (.txt) with a header line. It can only accept files including the following columns:
   1. `sample_ids` - Sample identifier
   2. `status` - Disease status (1 for cases, 0 for controls)
   3. `ageInt` - Age at interview
   4. `AgeDiagIndex` - Age at breast cancer diagnosis
   5. `Country` - Country of origin

    *The phenotype file should be like:*
   ```
   sample_ids  status    ageInt          AgeDiagIndex          Country
   64302       1         68              68                    Greece
   64303       1	      65	      68                    UK
   64304       0	      65	      68                    USA
   64305       0	      65	      68                    Australia
   64306       1	      61	      68                    UK
   64307       1	      62	      68                    USA
   64308       0	      61	      68                    Belgium
   ```
   **We emphasize that for *BRCA1* or *BRCA2* variants, the script will terminate if there exist samples in the input phenotype file for which age at diagnosis or interview are below 21 or above 80 years (this happens because age-specific breast cancer rates and risks are not available for ages below 21 or above 80 years).**

______________________________________________________________________________________________________________________________________________________________________________

3. A directory including genotype files in comma-delimited format (.csv) (one for each variant). The script can only accept files including the following two columns:
   1. `sample_ids` - Sample identifier
   2. `gt` - Sample genotype (0 for homozygous major (AA), 1 for heterozygous (Aa), 2 for homozygous minor (aa) and NA for samples with unknown genotype)
   
    *The genotype file should be like:*
   ```
   sample_ids  gt 
   64302       0  
   64303       0	
   64305       0	 
   64306       1	 
   64307       2	  
   64308       NA	   
   ```
______________________________________________________________________________________________________________________________________________________________________________

4. A comma-delimited file (.csv), which includes age-specific disease relative risks, cumulative incidence rates, penetrance and cumulative risk for the general population (non-carriers) and mutation-carriers. Breast cancer rates are published by [Dorling et al. 2021](https://www.nejm.org/doi/10.1056/NEJMoa1913948?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed), [Kuchenbaecker et al. 2017](https://jamanetwork.com/journals/jama/fullarticle/2632503), [Antoniou et al. 2003](https://www.cell.com/ajhg/fulltext/S0002-9297(07)60640-5). The GitHub archive contains the directory `rates` with the files required for *BRCA1* and *BRCA2* case-control LR calculations. These files must be supplied as inputs using the parameters: `--dorling` , `--kuchenbaecker` or `--antoniou`. 
 
    If you want to perform the case-control likelihood ratio calculations for variants in other genes then you need to create the file following the format shown below and supply it as input in the parameter `--customrates`:

  ```
  Age	Relative_risk_carriers	Cumulative_Incidence_Rates_carriers	Cumulative_Incidence_Rates_non-carriers	 Penetrance_carriers	Penetrance_non-carriers	    Cumulative_Risk_carriers	Cumulative_Risk_non-carriers
  21 	18	                19	                                0.007069091	                         0.007461818	        0.000392727	            0.007069091	                0.007461818
  22      18	                19	                                0.007069091	                         0.007461818	        0.000392727	            0.007069091	                0.007461818
  23      18	                19	                                0.007069091	                         0.007461818	        0.000392727	            0.007069091	                0.007461818
  24	18	                19	                                0.007069091	                         0.007461818	        0.000392727	            0.007069091	                0.007461818
  25	18	                19	                                0.007069091	                         0.007461818	        0.000392727	            0.007069091	                0.007461818
  ```
      
   So, according to the rates used, the user must also specify the parameters `--rates` as:
   
   `--rates Dorling` or `--rates Kuchenbaecker` or `--rates Antoniou` or `--rates Custom`

______________________________________________________________________________________________________________________________________________________________________________

5.  The analysis output directory and prefix.
    
    This can be set using the parameter `--output`. In case this parameter is not set, the output file names will have prefix `ccLR` and will be saved in the genotypes directory set with the parameter `--genotypes`.
    
## Running the script ##

**Finally to perform the analysis, use the following script:**

```
Rscript ccLR.R --gene BRCA2 --phenotype /ccLR/example_data/phenotypes.txt --genotypes /ccLR/example_data/genotypes --rates Dorling --dorling /ccLR/rates/Dorling_etal.2021.csv --output /ccLR/example_data/example
```
|Parameter      |Input       |
|:---    |:---   |
|`--gene`       |`BRCA1`, `BRCA2` or `Custom` |
|`--phenotype`  |phenotype file (.txt) with sample_ids, status, age at interview and age at disease diagnosis |
|`--genotypes`  |directory including genotype files (.csv), containing sample_ids and genotype |
|`--rates`      |`Dorling`, `Kuchenbaecker`, `Antoniou`,  `Custom` |
|`--dorling`    |file (.csv) of Dorling BC rates (parameter required only if `--rates Dorling` |
|`--kuchenbaecker`|file (.csv) of Kuchenbaecker BC rates (parameter required only if `--rates Kuchenbaecker` |
|`--antoniou`|file (.csv) of Antoniou BC rates (parameter required only if `--rates Antoniou` |
|`--customrates`|file (.csv) of custom disease rates (parameter required only if `--rates Custom` |
|`--output`|directory and prefix of the output file (Optional)|


## Excel Calculations ##
Detailed calculations can also be performed by the user, using the excel spreadsheet "ccLR.xlsx" included in the GitHub repository.

## General information ##
This method was developed by the [ENIGMA](https://enigmaconsortium.org/) (Evidence-based Network for the Interpretation of Germline Mutant Alleles) consortium for the calculation of the case-control likelihood ratio (LR) of *BRCA1* and *BRCA2* variants.
However, scripts can also be applied for other genes too. 
These scripts have been implemented in the context of the following "Under review" manuscript:

Zanti M et al. "Case-control likelihood ratio calculation for clinical classification of variants of uncertain significance in the *BRCA1* and *BRCA2* genes".

## Author ##
Zanti M.

## Contact ##
For questions and comments, please open a Github issue (preferred) or contact Maria Zanti at mariaz@cing.ac.cy




