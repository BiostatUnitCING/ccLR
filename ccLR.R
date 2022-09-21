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
phenotype=args$phenotype                                                              # location and name of phenotype file with 5 columns (sample_ids, status, age at interview, age at diagnosis and country of origin)
output=args$output                                                                  # location and name of file to save the output to; will be saved as 'output.csv'
rates=args$rates                                                                  # disease rates to be used (select between Dorling, Kuchenbaecker and Antoniou)
genotypes=args$genotypes                                                        # directory of input genotype files in csv format (2 columns needed: sample_ids and genotype) 
gene=args$gene                                                                # gene of interest (BRCA1, BRCA2 or Custom) 
customrates=args$customrates                                                #  directory of input files with penetrance and relative risks; required if we submit custom rates

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
colnum <- ncol(df) - 4
country_result <- NULL
country_result <- data.frame(country_result)
results <- NULL
results <- data.frame(results)
for (i in 2:colnum){
  paste(i)
  paste(list[[1]][i])
  snp_file <- paste('variant',list[[1]][i],".csv",sep="")
  df1 <- df[,c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex","StudyCountry")]
  df1 <- df1 %>% drop_na(list[[1]][i])
  df1$age_pen[df1$status=="0"] <- floor(df1$ageInt[df1$status=="0"])
  df1$age_pen[df1$status=="1"] <- floor(df1$AgeDiagIndex[df1$status=="1"])
  ##remove missing ages and filter for age
  df1 <- df1[which(df1$age_pen!=888),]
  if(any((df1$age_pen<21) | (df1$age_pen>80))){
    print("ERROR This analysis is only applicable for samples diagnosed or interviewed between the ages of 21-80. Please remove samples not following these criteria. Execution halted")
    quit()
  }
  non_carriers <- df1[which(df1[,list[[1]][i]]==0),]
  hets <- df1[which(df1[,list[[1]][i]]==1),]
  hets_cas <- hets[which(hets$status ==1),]
  if(nrow(hets_cas)!=0){
  meanagehets_cas <- round(mean(hets_cas$age_pen),digits=2)
  sdagehets_cas <- round(sd(hets_cas$age_pen),digits=2)
  minagehets_cas <- min(hets_cas$age_pen)
  maxagehets_cas <- max(hets_cas$age_pen)}else{
  meanagehets_cas <-NULL
  sdagehets_cas<-NULL
  minagehets_cas<-NULL
  maxagehets_cas<-NULL
  }
  hets_con <- hets[which(hets$status ==0),]
  if(nrow(hets_con)!=0){
  meanagehets_con <- round(mean(hets_con$age_pen),digits=2)
  sdagehets_con <- round(sd(hets_con$age_pen),digits=2)
  minagehets_con <- min(hets_con$age_pen)
  maxagehets_con <- max(hets_con$age_pen)}else{
  meanagehets_con <-NULL
  sdagehets_con<-NULL
  minagehets_con<-NULL
  maxagehets_con<-NULL
  }
  homo <- df1[which(df1[,list[[1]][i]]==2),]
  cas <- df1[which(df1$status == 1),]
  con <- df1[which(df1$status == 0),]
  freq_cases <- nrow(hets_cas)/nrow(cas)
  freq_controls <- nrow(hets_con)/nrow(con)

  if(nrow(hets)<MAF & nrow(hets)!=0 & (nrow(hets) > nrow(homo))){
    df1$StudyCountry <- as.factor(df1$StudyCountry)
    country <- levels(df1$StudyCountry)
    for (h in 1:length(country)) { # Germany - 8 code
      df2 <- df1[which(df1$StudyCountry == country[h]),]
      N=nrow(df2)
      ncon <- nrow(df2[df2$status == "0",])
      ncas <- nrow(df2[df2$status == "1",]) #number of cases and controls when NA were removed
      cases <- df2[which(df2$status==1),]
      controls <- df2[which(df2$status==0),]
      car1 <- df2[which(df2[,list[[1]][i]]==1),]
      ncar1 <- nrow(car1)
      car2 <-  df2[which(df2[,list[[1]][i]]==2),]
      ncar2 <- nrow(car2)
      m2 <- df2
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
      
        m2$RR[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Relative_risk_BRCA1.carriers"]
        m2$Incidence[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Incidence"]

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
        
        m2$RR[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Relative_risk_BRCA1.carriers"]
        m2$Incidence[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Incidence"]

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
        
        m2$RR[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Relative_risk_BRCA1.carriers"]
        m2$Incidence[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Incidence"]

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
          
          m2$RR[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Relative_risk_BRCA2.carriers"]
          m2$Incidence[j] <- Dorling[Dorling$Age==m2$age_pen[j],"BC_Incidence"]

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
          
          m2$RR[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Relative_risk_BRCA2.carriers"]
          m2$Incidence[j] <- Kuchenbaecker[Kuchenbaecker$Age==m2$age_pen[j],"BC_Incidence"]

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
          
          m2$RR[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Relative_risk_BRCA2.carriers"]
          m2$Incidence[j] <- Antoniou[Antoniou$Age==m2$age_pen[j],"BC_Incidence"]

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
            
            m2$RR[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Relative_risk_carriers"]
            m2$Incidence[j] <- Customrates[Customrates$Age==m2$age_pen[j],"Disease_Incidence"]
          }
        }else if (is.null(rates)){
          print("ERROR rates for LR calculation must be provided. Execution halted")
          quit()
        } 
      }else if (is.null(gene)){
      print("ERROR Gene for LR calculation must be provided. Execution halted")
      quit()
      }  
    
      #Remove homozygotes
      m2 <- m2[which(m2[,2]!=2),]
      m2$S0 <- exp((-(m2$Incidence))*m2$age_pen)
      m2$S1 <- exp((-(m2$Incidence))*m2$age_pen*m2$RR)
      K=ncar1 #total number of variant carriers (N defined earlier-stratum specific)
      m2$Likelihood <- ifelse(m2$status==0, ((m2$S1)/(m2$S0))*(m2$RR^0),((m2$S1)/(m2$S0))*(m2$RR^1))
      m3 <- m2[which(m2[,2] == 1),]
      LR <- ifelse(K>0, (N^K*prod(m3$Likelihood))/((sum(m2$Likelihood))^K), 1)
      country_result <- rbind(country_result)
      country_result[country[h],"K"] <- K
      country_result[country[h],"N"] <- N
      country_result[country[h],"Prod"] <- prod(m3$Likelihood)
      country_result[country[h],"Sum"] <- sum(m2$Likelihood)
      country_result[country[h],"LR_stratum"] <- LR 
    }
      
      SNP <- paste(list[[1]][i])
    
      results[SNP,"Total_n_cases"] <- nrow(cas)  
      results[SNP,"Total_n_controls"] <- nrow(con) 
      results[SNP,"n_carriers_cases"] <- nrow(hets_cas)
      results[SNP,"n_carriers_controls"] <- nrow(hets_con)
      results[SNP,"freq_cases"] <- freq_cases
      results[SNP,"freq_controls"] <- freq_controls
      if(nrow(hets_cas)>1){
        results[SNP,"Case_Carriers_Age"] <- paste(meanagehets_cas,"±",sdagehets_cas," (",minagehets_cas,"-",maxagehets_cas,")",sep="")}else{results[SNP,"Case_Carriers_Age"] <- meanagehets_cas}
      if(nrow(hets_con)>1){
        results[SNP,"Control_Carriers_Age"] <- paste(meanagehets_con,"±",sdagehets_con," (",minagehets_con,"-",maxagehets_con,")",sep="")}else{results[SNP,"Control_Carriers_Age"]<-NA}
      country_result$LR_stratum[is.na(country_result$LR_stratum)] <- 1
      results[SNP,"LR"] <- prod(country_result$LR_stratum)
  }
  else{
    i<-i+1}}

write.csv(results, paste(output,".csv",sep=""))


