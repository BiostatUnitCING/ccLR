# This code can be used to replicate the case-control likelihood ratio calculation
# This analysis implements age-specific relative risks or odds ratios of breast cancer published either by Dorling L et al. 2021, Kuchenbaecker KB et al. 2017 or Antoniou A et al. 2003

#Required R packages
library(tidyverse)
library(R.utils)

########## Process inputs ########## 
args=commandArgs(asValue=TRUE)

dorling=args$dorling                                                                        # location and name of input file of published data from Dorling et al. 2021 (These are supplied in the GitHub repository)
kuchenbaecker=args$kuchenbaecker                                                          # location and name of input file of published data from Kuchenbaecker et al. 2017 (These are supplied in the GitHub repository)
antoniou=args$antoniou                                                                  # location and name of input file of published data from Antoniou et al. 2003 (These are supplied in the GitHub repository)
phenotype=args$phenotype                                                              # location and name of phenotype file with 4 columns (sample_ids, status, age at interview and age at diagnosis)
output=args$output                                                                  # location and name of file to save the output to; will be saved as 'output.csv'
rates=args$rates                                                                  # disease rates to be used (select between Dorling, Kuchenbaecker and Antoniou)
genotypes=args$genotypes                                                        # directory of input genotype files in csv format (2 columns needed: sample_ids and genotype) 
gene=args$gene                                                                # gene of interest (BRCA1 or BRCA2) 
customrates=args$customrates                                                # gene of interest (BRCA1 or BRCA2) 

if(is.null(phenotype)){
  print("ERROR phenotype file must be provided. Execution halted")
  quit()
}

if(is.null(genotypes)){
  print("ERROR genotype files must be provided. Execution halted")
  quit()
}

if(is.null(output)){output="ccLR"}

pheno <- read.table(paste(phenotype), header = TRUE)
MAF <- nrow(pheno)/100

pheno$ageInt <- floor(pheno$ageInt)
pheno$AgeDiagIndex <- floor(pheno$AgeDiagIndex)

results <- NULL                        # creates a data frame for the analysis results  
results <- data.frame(results) 

ldf <- list()                          # creates a list
listcsv <- dir(paste(genotypes), pattern = "*.csv")  # creates the list of all the csv files in the genotypes directory
setwd(paste(genotypes))
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.csv(listcsv[k], header = T)
  colnames(ldf[[k]]) = c("sample_ids", listcsv[k])
  colnames(ldf[[k]]) <-  sub(".csv", "", colnames(ldf[[k]]))
}

ldf1 <- ldf %>% reduce(inner_join, by = "sample_ids")

DF <- merge(ldf1,pheno,by="sample_ids")
df <- DF
df$AgeDiagIndex[is.na(df$AgeDiagIndex)] <- 888
df$ageInt[is.na(df$ageInt)] <- 888

list <- list(colnames(df))
colnum <- ncol(df) - 3
for (i in 2:colnum){
  snp_file <- paste('variant',list[[1]][i],".csv",sep="")
  df1 <- df[,c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex")]
  df2 <- df1 %>% drop_na(list[[1]][i])
  ncon <- nrow(df2[df2$status == "0",])
  ncas <- nrow(df2[df2$status == "1",])
  cases <- df2[which(df2$status==1),]
  controls <- df2[which(df2$status==0),]
  carriers1 <- df2[which(df2[,2]==1),]
  ncar1 <- nrow(carriers1)
  carriers2 <- df2[which(df2[,2]==2),]
  ncar2 <- nrow(carriers2)
  
  if(ncar1<MAF & ncar1!=0 & (ncar1 > ncar2)){
    
    tca1 <- cases[,c("sample_ids",list[[1]][i])]
    #tca <- table(cases$sample,cases[,list[[1]][i]])
    #tca1 <- as.data.frame.matrix(tca)
    hetca <- tca1[which(tca1[,2]>0),]
    hetca1 <- hetca[which(hetca[,2]==1),]
    hetca2 <- hetca[which(hetca[,2]==2),]
    freq_ca <- nrow(hetca1)/ncas
    
    tco1 <- controls[,c("sample_ids",list[[1]][i])]
    #tco <- table(controls$sample,controls[,list[[1]][i]])
    #tco1 <- as.data.frame.matrix(tco)
    hetco <- tco1[which(tco1[,2]>0),]
    hetco1 <- hetco[which(hetco[,2]==1),]
    hetco2 <- hetco[which(hetco[,2]==2),]
    freq_con <- nrow(hetco1)/ncon
    
    #####
    
    prop_cases <- ncas/(ncas+ncon)
    
    m1 <- df2[which(df2[,list[[1]][i]]==1),]
    m1$age_pen[m1$status=="0"] <- floor(m1$ageInt[m1$status=="0"])
    m1$age_pen[m1$status=="1"] <- floor(m1$AgeDiagIndex[m1$status=="1"])
    
    ##remove missing ages and filter for age
    m2 <- m1[which(m1$age_pen!=888),]
    if(any((m2$age_pen<21) | (m2$age_pen>80))){
      print("ERROR This analysis is only applicable for samples diagnosed or interviewed between the ages of 21-80. Please remove samples not following these criteria. Execution halted")
      quit()
    }
    rownames(m2) = NULL
    
    if(gene=="BRCA1"){
      
    if(rates=="Dorling"){
      
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      #Breast Cancer Odds ratios taken from Dorling L et al.  N. Engl. J. Med. 2021;384:428-439.
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      if(is.null(dorling)){
        print("ERROR Dorling_etal.2021.csv file (located in the script folder) must be provided. Execution halted")
        quit()
      }
      
      Dorling <- read.csv(paste(dorling), header = TRUE)
    
      for (j in 1:nrow(m2)){
      
        m2$D[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
        m2$I[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
        m2$J[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
        m2$K[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
      
      }
    } else if (rates=="Kuchenbaecker"){
      
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      #Breast Cancer Relative Risk taken from Kuchenbaecker KB et al. JAMA. 2017;317(23):2402-2416.
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      if(is.null(kuchenbaecker)){
        print("ERROR Kuchenbaecker_etal.2017.csv file (located in the script folder) must be provided. Execution halted")
        quit()
      }
      
      Kuchenbaecker <- read.csv(paste(kuchenbaecker), header = TRUE)
      
      for (j in 1:nrow(m2)){
        
        m2$D[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
        m2$I[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
        m2$J[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
        m2$K[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
        
      }
      
    } else if(rates=="Antoniou"){
      
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      #Breast Cancer Relative Risk taken from Antoniou A et al. Am J Hum Genet. 2003;72(5):1117-1130. 
      #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      if(is.null(antoniou)){
        print("ERROR Antoniou_etal.2003.csv file (located in the script folder) must be provided. Execution halted")
        quit()
      }
      
      Antoniou <- read.csv(paste(antoniou), header = TRUE)
      
      for (j in 1:nrow(m2)){
        
        m2$D[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
        m2$I[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
        m2$J[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
        m2$K[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
        
      }
    } else if (is.null(rates)){
      print("ERROR disease rates required for LR calculation must be provided. Execution halted")
      quit()
    }
    } else if (gene=="BRCA2"){
      if(rates=="Dorling"){
        
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #Breast Cancer Odds ratios taken from Dorling L et al.  N. Engl. J. Med. 2021;384:428-439.
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        if(is.null(dorling)){
          print("ERROR Dorling_etal.2021.csv file (located in the script folder) must be provided. Execution halted")
          quit()
        }
        
        Dorling <- read.csv(paste(dorling), header = TRUE)
        
        for (j in 1:nrow(m2)){
          
          m2$D[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
          m2$I[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
          m2$J[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
          m2$K[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
          
        }
      } else if (rates=="Kuchenbaecker"){
        
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #Breast Cancer Relative Risk taken from Kuchenbaecker KB et al. JAMA. 2017;317(23):2402-2416.
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        if(is.null(kuchenbaecker)){
          print("ERROR Kuchenbaecker_etal.2017.csv file (located in the script folder) must be provided. Execution halted")
          quit()
        }
        
        Kuchenbaecker <- read.csv(paste(kuchenbaecker), header = TRUE)
        
        for (j in 1:nrow(m2)){
          
          m2$D[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
          m2$I[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
          m2$J[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
          m2$K[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
          
        }
        
      } else if(rates=="Antoniou"){
        
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #Breast Cancer Relative Risk taken from Antoniou A et al. Am J Hum Genet. 2003;72(5):1117-1130. 
        #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        if(is.null(antoniou)){
          print("ERROR Antoniou_etal.2003.csv file (located in the script folder) must be provided. Execution halted")
          quit()
        }
        
        Antoniou <- read.csv(paste(antoniou), header = TRUE)
        
        for (j in 1:nrow(m2)){
          
          m2$D[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Cumulative_Risk_BRCA1.carriers"]
          m2$I[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Cumulative_Risk_non.carriers"]
          m2$J[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Penetrance_BRCA1.carriers"]
          m2$K[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Penetrance_non.carriers"]
          
        }
      }
      }else if (gene=="Custom"){
        if (rates=="Custom"){
          if(is.null(customrates)){
            print("ERROR custom rates file for LR calculation must be provided. Execution halted")
            quit()
          }
          Customrates <- read.csv(paste(customrates), header = TRUE)
          for (j in 1:nrow(m2)){
            
            m2$D[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Cumulative_Risk_carriers"]
            m2$I[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Cumulative_Risk_non-carriers"]
            m2$J[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Penetrance_carriers"]
            m2$K[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Penetrance_non-carriers"]
          }
        }else if (is.null(rates)){
          print("ERROR rates for LR calculation must be provided. Execution halted")
          quit()
        } 
      }else if (is.null(gene)){
      print("ERROR Gene for LR calculation must be provided. Execution halted")
      quit()
      }  
    
    if(freq_con==0){
        freq_con <- (0 + 1) / (ncon + 2)
      }
    
      m2$E <- freq_con * (1-m2$D) / (freq_con * (1-m2$D) + (1-freq_con)*(1-m2$I))
      m2$F <- freq_con * m2$J / (freq_con * m2$J + (1 - freq_con) * m2$K)
      m2$G <- m2$F / m2$E
      m2$H <- prop_cases * m2$G / (prop_cases * m2$G + 1 - prop_cases)
      m2$L <- ifelse(m2$status==1,m2$H,1-m2$H)
      m2$M <- log10(m2$L)
  
      F16 <- (prop_cases^nrow(hetca1))*((1-prop_cases)^nrow(hetco1))
      H14 <- sum(m2$M)
      H16 <- log10(F16)
      H18 <- H14-H16
      I18 <- 10^(H18)
      J18 <- 10^(-H18)
    
      SNP <- paste(list[[1]][i])
    
      results[SNP,"Total_n_cases"] <- ncas
      results[SNP,"Total_n_controls"] <- ncon
      results[SNP,"n_carriers_cases"] <- nrow(hetca1)
      results[SNP,"n_carriers_controls"] <- nrow(hetco1)
      results[SNP,"freq_cases"] <- freq_ca
      results[SNP,"freq_controls"] <- freq_con
      results[SNP,"oddsfavour"] <- I18
      results[SNP,"oddsagainst"] <- J18
      print(paste("Case-control likelihood ratio calculated with success for variant",list[[1]][i], sep=" "))
      }
  else{
  i<-i+1}}

write.csv(results, paste(output,".csv",sep=""))


