# Case-control likelihood ratio (LR) #
**ccLR** (Case-Control Likelihood Ratio)

This repository contains freely accessible R scripts and pre-formatted excel calculators for the calculation of case-control likelihood ratios (LRs). 

This method compares the likelihood of the distribution of the variant of interest among cases and controls, under the hypothesis that the variant is associated with similar risks of the disease in question, as the “average” pathogenic variant, compared to the likelihood under the hypothesis that it is a benign variant not associated with increased risk. 

These resources may be readily applied for the calculation of LRs to be used in classification of VUS in the *BRCA1* and *BRCA2*, and other disease susceptibility genes with known penetrance values. 

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
   5. `StudyCountry` - Country of origin

    *The phenotype file should be like:*
   ```
   sample_ids  status    ageInt          AgeDiagIndex          StudyCountry
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

4. A comma-delimited file (.csv), which includes age-specific disease relative risks, and disease incidence for the general population. Breast cancer rates are published by [Dorling et al. 2021](https://www.nejm.org/doi/10.1056/NEJMoa1913948?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed), [Kuchenbaecker et al. 2017](https://jamanetwork.com/journals/jama/fullarticle/2632503), [Antoniou et al. 2003](https://www.cell.com/ajhg/fulltext/S0002-9297(07)60640-5). The GitHub archive contains the directory `rates` with the files required for *BRCA1* and *BRCA2* case-control LR calculations. These files must be supplied as inputs using the parameters: `--dorling` , `--kuchenbaecker` or `--antoniou`. 
 
    If you want to perform the case-control likelihood ratio calculations for variants in other genes then you need to create the file following the format shown below and supply it as input in the parameter `--customrates`:

  ```
  Age	Relative_risk_carriers		Disease_Incidence	    	
  21 	18	                        0.0000392727	            
  22      18	                        0.0000392727	           
  23      18	                        0.0000392727	       
  24	18	                        0.0000392727	           
  25	18	                        0.0000392727	           
  ```
      
   So, according to the rates used, the user must also specify the parameters `--rates` as:
   
   `--rates Dorling` or `--rates Kuchenbaecker` or `--rates Antoniou` or `--rates Custom`

______________________________________________________________________________________________________________________________________________________________________________

5.  The analysis output directory and prefix.
    
    This can be set using the parameter `--output`. In case this parameter is not set, the output file names will have prefix `ccLR` and will be saved in the genotypes directory set with the parameter `--genotypes`.
    
## Running the script ##

**Finally to perform the analysis, use the following script:**

```
Rscript ccLR.R --gene BRCA2 --phenotype ~/ccLR-main/example_data/phenotypes.txt --genotypes ~/ccLR-main/example_data/genotypes --rates Dorling --dorling ~/ccLR-main/rates/Dorling_etal.2021.csv --output ~/ccLR-main/example_data/example
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
Case-control likelihood ratio calculations can also be performed, using the pre-formatted excel calculators. User can input either individual-level data ("ccLRCalculator.xlsx") or tabulated age groups ("ccLRCalculator_TabulatedFormat.xlsx").

## General information ##
This method was developed on behalf of the [Breast Cancer Association Consortium (BCAC)](https://bcac.ccge.medschl.cam.ac.uk/) and the [Evidence-based Network for the Intepretation of Germline Mutant Alleles (ENIGMA)](https://enigmaconsortium.org/) for the calculation of the case-control likelihood ratio (LR) for *BRCA1*, *BRCA2* and other disease-genes with well-established penetrance estimates.

These scripts have been implemented in the context of the following manuscripts:

Zanti M et al. A likelihood ratio approach for utilizing case-control data in the clinical classification of rare sequence variants: application to *BRCA1* and *BRCA2* https://onlinelibrary.wiley.com/doi/10.1155/2023/9961341

## Author ##
Zanti M.

## Contact ##
For questions and comments, please open a Github issue (preferred) or contact Maria Zanti at mariaz at cing.ac.cy

## License ##
The method is distributed under the GPL either version 3 or later license.






