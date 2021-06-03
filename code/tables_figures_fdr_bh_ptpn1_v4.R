#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file




##################################################################################
####################### SCRIPT FOR PREPARING TABLES AND FIGURES ##################
##################################################################################

#Script where we create all the figures and tables of the manuscript




####################################################
##### CHANGES RESPECT TO PREVIOS VERSIONS ##########
####################################################

#Respect previous papers
    #Here, I changed the code for table 2 in lpl code to include some suggestions for reviewers of cntf

#Version 2
    #The version 2 is run in the David Enard laptop. I changed the titles of the supplementary figures and calculate the correlation between leptin and adiposity.
    #also remove the binding of supplementary pdfs, now this is not needed.

#Version 3
    #In the version 3 I have changed the name of the supplementary figures 2,3,4. Now it is called figure SX...

#Version 4
    #In the version 4, I have included the reviewer suggestions of ped obesity.
    #New versions of the 2 tables, and complete removal of all figures




########################################
############ DATA PREPARATION ##########
########################################

### load the environment with the analyses run
load("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/rdata/analysis_v2.RData")
require(SNPassoc)
require(genetics)


### set the working directory
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7")


### create a data.frame withe the names of phenotypes in dataset and the real name per the figure
pheno_names = cbind.data.frame(
    c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap","obesity","CRF_BMI","CRF_waist","waist_height","CRF_hip","waist_hip","CRF_Body_fat_PC","FMI", "TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","HOMA","QUICKI","Leptin_ng_ml","SBP","DBP"),
    c("Center (%)", "Age (y)", "Weight (kg)", "Height (cm)", "Triceps skinfold (mm)", "Subescapular skinfold (mm)", "% obesity/overweight", "BMI (kg/m^2)", "Waist circum. (cm)", "Waist/Height ratio", "Hip circum. (cm)", "Waist/Hip ratio", "Body fat (%)", "FMI (kg/m^2)", "Total cholesterol (mg/dL)","LDL-C (mg/dL)","HDL-C (mg/dL)","Total cholesterol/HDL-C","LDL-C/HDL-C","Triglycerides (mg/dL)","Triglycerides/HDL-C","ApoA1 (mg/dL)","ApoB (mg/dL)","ApoB/ApoA1","ApoB/LDL-C","Insulin (micro lU/mL)","HOMA","QUICKI","Leptin (ng/ml)","SBP (mm Hg)","DBP (mm Hg)"))
colnames(pheno_names) <- c("var_name", "final_name")

### binomial phenotypes included in this set of results
binomial_pheno = c("obesity", "CVi_BP")


### the same for SSB variable
ssb_names = cbind.data.frame(
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"),
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"))
colnames(ssb_names) <- c("var_name", "final_name")


### level names of factor variables
factor_variables_levels = list(c("Non-overweight", "Overweight"), c("Low physic. act. (<60 min/day)", "High physic. act. (≥60 min/day)"))
names(factor_variables_levels) <- c("obesity", "PA_factor")


### load allele and genes names of SNPs according to HELENA and ncbi
#gene names
gene_names = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/chromosome_snps.csv", sep=",", header=T)
#select snps from the studied gene
gene_names = gene_names[which(gene_names$selected_snp %in% labels(myData_ptpn1)),]
nrow(gene_names) == length(labels(myData_ptpn1))

#load allele names
#alleles = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/alleles_ptpn1_v2.csv", sep=",", header=T)
    #allele names have been updated in the main script of the analyses after checking allele frequencies with the 1KGP.
nrow(alleles) == length(labels(myData_ptpn1))


##check no individual has NA for all phenotypes and that no phenotype has NA for all individuals
#select phenotype columns    
only_pheno = myData_ptpn1[, which(!colnames(myData_ptpn1) %in% labels(myData_ptpn1))]
#which rows (individuals) have NA for all columns (phenotypes)
length(which(rowSums(is.na(only_pheno)) == ncol(only_pheno))) == 0 
length(which(colSums(is.na(only_pheno)) == nrow(only_pheno))) == 0 
    #first test for rows (individuals) with all column entries (phenotypes) as NA
    #second test for columns (phenotypes) with all rows entries (individuals) as NA
    #it should be zero in both cases




###################
##### TABLE 1 #####
###################

### start the table
#extract summary of snps
summary_snps = summary(myData_ptpn1)

#extract the major allele frequency and convert to minor allele frequency
MAF = 1-(summary_snps[, which(names(summary_snps) == "major.allele.freq")]/100)

#extract the HWE p.values
HWE = summary_snps[, which(names(summary_snps) == "HWE")]

#bind both metrics
table_1 = cbind.data.frame(MAF, HWE)
rownames(table_1) <- row.names(summary_snps)
table_1

#check that the HWE and MAF correspond to the same snps
table_1$MAF == 1-(summary_snps$major.allele.freq/100)
table_1$HWE == summary_snps$HWE


### add alleles to table 1

#add two columns for Major and minor alleles
table_1 = cbind.data.frame(NA, NA, table_1)

#add the corresponding names (first major, second minor)
colnames(table_1)[which(colnames(table_1) == "NA")] <- c("Major allele", "Minor allele")

#add the alleles
for(i in 1:nrow(table_1)){#for each row of table 1

    #select the [i] row
    selected_row = table_1[i,]

    #select the snp of the [i] row
    selected_snp = row.names(selected_row)

    #extract the ncbi snps from [i] row
    ncbi_name = alleles[which(alleles$snp == selected_snp),]$ncbi
    
    #select major and minor allele (they are ordered, major always first)
    ncbi_major = strsplit(as.character(ncbi_name), split="/")[[1]][1]
    ncbi_minor = strsplit(as.character(ncbi_name), split="/")[[1]][2]

    #save them in the corresponding columns
    table_1[i,which(colnames(table_1) == "Major allele")] <- ncbi_major
    table_1[i,which(colnames(table_1) == "Minor allele")] <- ncbi_minor
}

#extract the allele names from the table 1 (combined major and minor) along with snp names
alleles_from_table_1 = cbind.data.frame(row.names(table_1), paste(table_1[,which(colnames(table_1) == "Major allele")], "/", table_1[,which(colnames(table_1) == "Minor allele")], sep=""))
colnames(alleles_from_table_1) <- c("snp", "alleles_from_ncbi")

#merge these alleles names with the original allele names from alleles df
df_check_alleles = merge(alleles_from_table_1, alleles)

#check that the colums of alleles names are similar
identical(df_check_alleles$alleles_from_ncbi, df_check_alleles$ncbi)
    #CHECK THIS LINE, IT GIVES FALSE!!!

#reorder the table following the order in the chromosome 
table_1 = table_1[match(ptpn1_snps$snp, row.names(table_1)),]
row.names(table_1) == ptpn1_snps$snp


### add the gene name to the table. More sense when SNPs of multiple genes are included.
#loop for each gene
ptpn1s = unique(ptpn1_snps$gen)
#for each gene name
for(i in 1:length(ptpn1s)){

    #select the [i] gene
    ptpn1_selected = ptpn1s[i]

    #select the position of the first snp of the [i] gene
    position_to_add = which(row.names(table_1) == ptpn1_snps[ptpn1_snps$gen==ptpn1_selected,][1,]$snp)

    #add the gene name
    if(position_to_add == 1){# if the position is the 1

        #add in the position of MAF and HWE two NAs in a row
        row_to_add = cbind(NA, NA, NA, NA)
        colnames(row_to_add) <- c("Major allele", "Minor allele", "MAF", "HWE") #add the same names that in the data.frame

        #convert to upper case
        row.names(row_to_add) <- toupper(ptpn1_selected)

        #bind first the gene name
        table_1 = rbind(row_to_add, table_1) 
    } else{ #if the gene name have to be include in the middle of the table
        
        #add in the position of MAF and HWE two NAs in a row
        row_to_add = cbind(NA, NA, NA, NA)
        colnames(row_to_add) <- c("Major allele", "Minor allele", "MAF", "HWE")#add the same names that in the data.frame

        #convert to upper case        
        row.names(row_to_add) <- toupper(ptpn1_selected)

        #add the gene name in the middle of the table, after the prior gene and before the snps of this gene [i]        
        table_1 = rbind(table_1[1:(position_to_add-1),], row_to_add, table_1[(position_to_add):nrow(table_1),])        
    }
}


### add the minor alleles frequencies according to 1000 Genomes Project

#Minor allele frequencies obtained according to 1000 KGP obtained from the NCBI
    #rs6067472: T=0.3648
        #https://www.ncbi.nlm.nih.gov/snp/rs6067472#frequency_tab
    #rs10485614: C=0.0755
        #https://www.ncbi.nlm.nih.gov/snp/rs10485614#frequency_tab
    #rs2143511: C=0.4543
        #https://www.ncbi.nlm.nih.gov/snp/rs2143511#frequency_tab
    #rs6020608: T=0.2744
        #https://www.ncbi.nlm.nih.gov/snp/rs6020608#frequency_tab
    #rs968701: G=0.4861
        #https://www.ncbi.nlm.nih.gov/snp/rs968701#frequency_tab

#create a table with these minor allele frequencies
maf_1kgp = rbind.data.frame(NA, 0.3648, 0.0755, 0.4543, 0.2744, 0.4861)
row.names(maf_1kgp) <- c("PTPN1", "rs6067472", "rs10485614", "rs2143511", "rs6020608", "rs968701")
colnames(maf_1kgp) <- c("MAF 1KGP - Europe")

#merge with table 1 
table_1_1KGP = merge(table_1, maf_1kgp, by="row.names")

#reorder the table again following the order in the chromosome 
table_1_1KGP = table_1_1KGP[match(c("PTPN1", as.character(ptpn1_position$snp)), table_1_1KGP$Row.names),]
table_1_1KGP$Row.names == c("PTPN1", as.character(ptpn1_position$snp))

#remove the column name for row names
colnames(table_1_1KGP)[which(colnames(table_1_1KGP) == "Row.names")] <- ""




####################
##### TABLE 2 ######
####################

#This table 2 calculates the number of individuals in each center and then the proportion of overweight WITHIN each center. In that way, we can compare the prevalence of overweight while considering the differences in sample size between centers. 

#extract the list of centers
list_centers = unique(myData_ptpn1$center)

#create the empty table
table_2 = data.frame(cbind(NA, NA, NA, NA, NA, NA, NA))
names(table_2)[1] <- "Center"
names(table_2)[2] <- paste("All sample size", sep="")
names(table_2)[3] <- paste("All overweight (%)", sep="")
names(table_2)[4] <- paste("Male sample size", sep="")
names(table_2)[5] <- paste("Male overweight (%)", sep="")
names(table_2)[6] <- paste("Female sample size", sep="")
names(table_2)[7] <- paste("Female overweight (%)", sep="")

#for each phenotype
for (i in 1:length(list_centers)){

    #select the [i] phenotype
    center_selected = list_centers[i]


    ## data from both sexes
    #the total number of individuals of both sexes from the [i] center that have data for obesity
    all_sample_size = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$center == '", center_selected, "' & !is.na(myData_ptpn1$obesity)),]", sep=""))))
    #the total number of individuals of both sexes from the [i] center that are overweight
    all_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$center == '", center_selected, "' & myData_ptpn1$obesity == 1),]", sep=""))))


    ## data from males
    #the total number of males from the [i] center that have data for obesity
    male_sample_size = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$center == '", center_selected, "' & !is.na(myData_ptpn1$obesity)),]", sep=""))))
    #the total number of males from the [i] center that are overweight
    male_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$center == '", center_selected, "' & myData_ptpn1$obesity == 1),]", sep=""))))


    ## data from females
    #the total number of females from the [i] center that have data for obesity
    female_sample_size = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$center == '", center_selected, "' & !is.na(myData_ptpn1$obesity)),]", sep=""))))
    #the total number of females from the [i] center that are overweight    
    female_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$center == '", center_selected, "' & myData_ptpn1$obesity == 1),]", sep=""))))


    #bind results in one row
    results = cbind.data.frame(
        center_selected,
        all_sample_size,
        round((all_overweight / all_sample_size)*100, 2),
        male_sample_size,
        round((male_overweight / male_sample_size)*100, 2),
        female_sample_size,
        round((female_overweight / female_sample_size)*100, 2))


    #change columns names to match names table 2
    names(results)[1] <- "Center"
    names(results)[2] <- paste("All sample size", sep="")
    names(results)[3] <- paste("All overweight (%)", sep="")
    names(results)[4] <- paste("Male sample size", sep="")
    names(results)[5] <- paste("Male overweight (%)", sep="")
    names(results)[6] <- paste("Female sample size", sep="")
    names(results)[7] <- paste("Female overweight (%)", sep="")

    #add to the table
    table_2 = rbind(table_2, results)
}

#remove the first row
table_2 = table_2[-which(rowSums(is.na(table_2)) == ncol(table_2)),]

#remove the rows of obesity phenotype ("% obesity/overweight"; we now have columns separated by overweight) and two centers for which we don't have data: Modena and Birmingham
table_2 = table_2[which(!table_2$Center %in% c("Birmingham* in UK", "Modena (Italy)")),]

#see the table
table_2


##for the phenotypes with percentage (%) or squared (^2), make some changes to be acceptable for latex.

#We have to add slash (\\) to avoid problems in latex (two because one is an expression for R).

#pheno with slash
column_slash = which(grepl("%", colnames(table_2))) #rows with percentage as phenotype name
#change names of these phenotypes modifying "%" by "\\%"
colnames(table_2)[column_slash] <- gsub("%", "\\%", colnames(table_2)[column_slash], fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression




####################
##### TABLE 3 ######
####################

#This table 3 calculated the mean +- SD for different phenotypes in each category of sex*weight status. If you want scripts for calculating percentage of discrete factors inside each category, see previous versions of this script.

#select the phenotypes to show in table 3
response_pheno_table_3 = c("CRF_age", "CRF_weight", "CRF_height", "CRF_BMI","CRF_waist","waist_height","CRF_hip","waist_hip","CRF_Body_fat_PC","FMI", "TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","HOMA","QUICKI","Leptin_ng_ml","SBP","DBP")
length(response_pheno_table_3) == nrow(pheno_names[which(!pheno_names$var_name %in% c("center",  "CRF_trici", "CRF_subscap", "obesity")),])

#calculate n for both sexes
n_all = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("male", "female")),])
n_male = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex=="male"),])
n_female = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex=="female"),])

#calculate n for sex and obesity
n_all_healthy_weight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("male", "female") & myData_ptpn1$obesity == 0),])
n_all_overweight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("male", "female") & myData_ptpn1$obesity == 1),])
n_male_healthy_weight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("male") & myData_ptpn1$obesity == 0),])
n_male_overweight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("male") & myData_ptpn1$obesity == 1),])
n_female_healthy_weight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("female") & myData_ptpn1$obesity == 0),])
n_female_overweight = nrow(myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c("female") & myData_ptpn1$obesity == 1),])

#check that overwight and healthi weight in each group are correct
n_all_healthy_weight+n_all_overweight == n_all
n_male_healthy_weight+n_male_overweight == n_male
n_female_healthy_weight+n_female_overweight == n_female
    #WARNING!! THIS TOTAL NUMBER OF INDNVIDUALS DOES OT TAKE INTO ACCOUNT THE EXISTENCE OF NAs IN SPECIFIC PHENOTYPES. This is the number of individuals for which we have data for most of studied SNPs and phenotypes.

#create the empty table
table_3 = data.frame(cbind(NA, NA, NA, NA, NA, NA, NA))
names(table_3)[1] <- "Phenotype"
names(table_3)[2] <- paste("All non-overweight (n=", n_all_healthy_weight, ")", sep="")
names(table_3)[3] <- paste("All overweight (n=", n_all_overweight, ")", sep="")
names(table_3)[4] <- paste("Male non-overweight (n=", n_male_healthy_weight, ")", sep="")
names(table_3)[5] <- paste("Male overweight (n=", n_male_overweight, ")", sep="")
names(table_3)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
names(table_3)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

#for each phenotype
for (i in 1:length(response_pheno_table_3)){

    #select the [i] phenotype
    pheno_selected = response_pheno_table_3[i]

    #extract the complete name
    pheno_table = pheno_names$final_name[which(pheno_names$var_name == pheno_selected)]

    #set the number decimals
    #some phenotypes will have 1 decimal
    if(pheno_selected %in% c()){

        #1 decimal
        number_decimals = 1
    } else {#if not

        #and the phenotype is
        if(pheno_selected %in% c()){

            #0 decimals
            number_decimals = 0
        } else {#if the phenotype is none of the latter

            #2 decimals
            number_decimals = 2
        }
    }

    #calculate the mean and sd across the whole panel in healthy weight
    all_mean_healthy_weight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)
    all_sd_healthy_weight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)

    #calculate the mean and sd across the whole panel in overweight individuals
    all_mean_overweight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)
    all_sd_overweight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)

    #calculate the mean and sd across male with healthy weight
    male_mean_healthy_weight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)
    male_sd_healthy_weight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)

    #calculate the mean and sd across the whole panel in overweight individuals
    male_mean_overweight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)
    male_sd_overweight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)

    #calculate the mean and sd across female with healthy weight
    female_mean_healthy_weight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)
    female_sd_healthy_weight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))),number_decimals)

    #calculate the mean and sd across the whole panel in overweight individuals
    female_mean_overweight = round(mean(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)
    female_sd_overweight = round(sd(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))),number_decimals)

    #bind results into one row
    results = cbind.data.frame(
        pheno_table,
        paste(all_mean_healthy_weight, "$\\pm$", all_sd_healthy_weight, sep=""),
        paste(all_mean_overweight, "$\\pm$", all_sd_overweight, sep=""),
        paste(male_mean_healthy_weight, "$\\pm$", male_sd_healthy_weight, sep=""),
        paste(male_mean_overweight, "$\\pm$", male_sd_overweight, sep=""),
        paste(female_mean_healthy_weight, "$\\pm$", female_sd_healthy_weight, sep=""),
        paste(female_mean_overweight, "$\\pm$", female_sd_overweight, sep=""))

    #change columns names to match names table 2
    names(results)[1] <- "Phenotype"
    names(results)[2] <- paste("All non-overweight (n=", n_all_healthy_weight, ")", sep="")
    names(results)[3] <- paste("All overweight (n=", n_all_overweight, ")", sep="")
    names(results)[4] <- paste("Male non-overweight (n=", n_male_healthy_weight, ")", sep="")
    names(results)[5] <- paste("Male overweight (n=", n_male_overweight, ")", sep="")
    names(results)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
    names(results)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

    #bind to the table
    table_3 = rbind(table_3, results)
}

#remove the first row
table_3 = table_3[-which(rowSums(is.na(table_3)) == ncol(table_3)),]

#see the table
table_3


##for the phenotypes with percentage (%) or squared (^2), make some changes to be acceptable for latex.

#We have to add slash (\\) to avoid problems in latex (two because one is an expression for R).

#pheno with slash
pheno_slash = which(grepl("%", table_3$Phenotype)) #rows with percentage as phenotype name
#change names of these phenotypes modifying "%" by "\\%"
table_3[pheno_slash,]$Phenotype <- gsub("%", "\\%", table_3[pheno_slash,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression

#pheno with ^2
pheno_squared = which(grepl("\\^2", table_3$Phenotype)) #rows with squared as phenotype name
#change names of these phenotypes modifying "^2" by "\\textsuperscript{2}"
table_3[pheno_squared,]$Phenotype <- gsub("^2", "\\textsuperscript{2}", table_3[pheno_squared,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression




####################
##### TABLE 4 ######
####################

#This table will show the average values of phenotypes per each genotype along with the FDR of the additive and codominant model.


### check if considering additive/codominant models cover all significant associations ###

#copy the supple data to do some operations
crude_assocs = suppl_data_1
    #we use the supple dataset file because it has the final phenotypes and columns we need

#create a variable with the combination of phenotype and snp of each association
crude_assocs$pheno_snp_combination = interaction(crude_assocs$phenotype, crude_assocs$snp, sep="-")
#check
summary(crude_assocs$pheno_snp_combination == paste(crude_assocs$phenotype, crude_assocs$snp, sep="-"))

#select those associations with an FDR<0.1
assoc_fdr_less_01 = crude_assocs[which(crude_assocs$fdr<0.1),]

#from these associations, select those with the additive and codominant model
assoc_fdr_less_01_only_add_cod = crude_assocs[which(crude_assocs$fdr<0.1 & crude_assocs$heritage_model %in% c("additive", "codominant")),]
#check
nrow(assoc_fdr_less_01_only_add_cod[which(assoc_fdr_less_01_only_add_cod$fdr<0.1 & assoc_fdr_less_01_only_add_cod$heritage_model %in% c("additive", "codominant")),]) == nrow(assoc_fdr_less_01_only_add_cod)

#check whether all significant combinations are included when restricting to add/codominant
summary(unique(assoc_fdr_less_01$pheno_snp_combination) %in% unique(assoc_fdr_less_01_only_add_cod$pheno_snp_combination))

#select those not included
assoc_fdr_less_01[which(!assoc_fdr_less_01$pheno_snp_combination %in% assoc_fdr_less_01_only_add_cod$pheno_snp_combination),]
    #The only phenoype*snp combination with FDR<0.1 not covered by additive/codominant models is Leptin_ng_ml.rs10485614. This is only association for leptin with an FDR below 0.1 (FDR=0.075). This is not a relevant association for the manuscript, and in any case, it is mentioned in the discussion along with the genotype with higher leptin, corresponding FDR and R2. We can leave it that way. 


### make a loop to extract phenotype values per genotype ###

#open empty data.frame to save results
table_4 = rbind.data.frame(rep(NA, 11))
colnames(table_4)[1] <- "SNP"
colnames(table_4)[2] <- "Phenotype"
colnames(table_4)[3] <- "11"
colnames(table_4)[4] <- "12"
colnames(table_4)[5] <- "22"
colnames(table_4)[6] <- "P add" 
colnames(table_4)[7] <- "FDR add"
colnames(table_4)[8] <- "R2 add"
colnames(table_4)[9] <- "P cod"
colnames(table_4)[10] <- "FDR cod"
colnames(table_4)[11] <- "R2 cod"

#reorder the file with significant add/cod associations based on snps
assoc_fdr_less_01_only_add_cod_copy = assoc_fdr_less_01_only_add_cod[order(assoc_fdr_less_01_only_add_cod$snp),]

#select those pheno-snp combinations for associations with FDr<0.1 in the additive or codominant models
pheno_snp_combinations_table_4 = unique(assoc_fdr_less_01_only_add_cod_copy$pheno_snp_combination)
    #we need unique because an snp-pheno association that is significant for additive and codominant will be selected two times and included two times in the table. We only need it one time, if the combination is present two times in "assoc_fdr_less_01_only_add_cod_copy" (2 rows), the script will take the FDR of additive and codominant. If it is not present it will look for the non-significant model in "crude_assocs". 

#for each pheno-snp combination, considering only associations with FDR<0.1 and additive model
for(i in 1:length(pheno_snp_combinations_table_4)){

    #select the [i] combination
    selected_combination = pheno_snp_combinations_table_4[i]

    #check
    print(paste("##############################################"))
    print(paste("STARTING WITH TABLE 4: ", selected_combination, sep=""))
    print(paste("##############################################"))

    #select the rows of the significant associations for the [i] pheno-snp combination
    selected_rows = assoc_fdr_less_01_only_add_cod_copy[which(assoc_fdr_less_01_only_add_cod_copy$pheno_snp_combination == selected_combination),]
        #We can have the same phenotype-snp combination for two models, additive and codominant.

    #select the [i] phenotype and snp
    selected_pheno = unique(selected_rows$phenotype)
    selected_snp = unique(selected_rows$snp)
        #unique because we can have two rows (additive and codominant), but in both cases the snp and phenotype should be the same indicated in selected_combination
    #check
    print(paste("############################"))
    print(paste("WE SELECTED THE CORRECT PHENOTYPE AND SNP?"))
    print(paste(selected_pheno, "-", selected_snp, sep="") == selected_combination)
    print(paste("############################"))

    #extract the genotype data of the [i] snp
    geno_data = eval(parse(text=paste("na.omit(myData_ptpn1$", selected_snp, ")", sep="")))

    #extract genotype levels
    geno_data_levels = unique(geno_data)

    #for each genotype level calculate the number of individuals
    genotype_sample_size = data.frame(selected_genotype=NA, sample_size=NA) #open an empty data.frame to save resuls
    for(j in 1:length(geno_data_levels)){

        #select the [j] genotype
        selected_genotype = geno_data_levels[j]

        #calculate the number of individuals with the [j] snp
        sample_size = length(which(geno_data == selected_genotype))

        #save the results
        genotype_sample_size = rbind.data.frame(genotype_sample_size, cbind.data.frame(selected_genotype, sample_size))
    }

    #remove the first row with all NAs
    genotype_sample_size = genotype_sample_size[-which(rowSums(is.na(genotype_sample_size)) == ncol(genotype_sample_size)),]

    #select genotypes of the homozygotes
    genotype_sample_size_homo = genotype_sample_size[which(!genotype_sample_size$selected_genotype %in% c("1/2", "2/1")),]

    #select the homozygote genotype with the lowest sample size, that is, minor homozygote
    minor_homo = genotype_sample_size_homo[which(genotype_sample_size_homo$sample_size == min(genotype_sample_size_homo$sample_size)),]$selected_genotype

    #select the homozygote genotype with the highest sample size, that is, major homozygote
    major_homo = genotype_sample_size_homo[which(genotype_sample_size_homo$sample_size == max(genotype_sample_size_homo$sample_size)),]$selected_genotype

    #select those individuals with each type of genotype
    subset_minor_homo = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", selected_snp, " == '", minor_homo, "' ),]", sep="")))
    subset_hetero = eval(parse(text=paste("myData_ptpn1[which(!myData_ptpn1$", selected_snp, " %in% c('", minor_homo, "', '", major_homo, "') & !is.na(myData_ptpn1$", selected_snp, ")),]", sep="")))
    subset_major_homo = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", selected_snp, " == '", major_homo, "' ),]", sep="")))
        #"parse()" returns the parsed but unevaluated expressions in an "expression" object. Therefore, the code is not run yet.
        #"eval()" evaluates an R expression. It runs the expression. 
            #dummy example:
                #expression_sum = parse(text="1+1") #create a expression to sum 1+1
                #eval(expression_sum) #evaluate the expression, giving 2.

    #check the subset worked
    print(paste("############################"))
    print(paste("CHECK THE SUBSET WAS WELL FOR: ", selected_combination, sep=""))
    print(paste("############################"))
    #no other genotype rather than minor homo should exist in subset_minor_homo
    print(nrow(eval(parse(text=paste("subset_minor_homo[which(subset_minor_homo$", selected_snp, " != '", minor_homo, "'),]", sep="")))) == 0)
    #no other genotype rather than hetero should exist in subset_hetero
    print(nrow(eval(parse(text=paste("subset_hetero[which(subset_hetero$", selected_snp, " %in% c('", minor_homo, "', '", major_homo, "') | is.na(subset_hetero$", selected_snp, ")),]", sep="")))) == 0)
    #no other genotype rather than major homo should exist in subset_major_homo
    print(nrow(eval(parse(text=paste("subset_major_homo[which(subset_major_homo$", selected_snp, " != '", major_homo, "'),]", sep="")))) == 0)

    #set the number decimals
    #and the phenotype is
    if(selected_pheno %in% c()){

        #0 decimals
        number_decimals = 0
    } else {#if the phenotype is none of the latter

        #2 decimals
        number_decimals = 2
    }

    #extract the data of the selected phenotype to make a condition
    selected_pheno_data = eval(parse(text=paste("myData_ptpn1$", selected_pheno, sep="")))

    #if the selected phenotype is not a factor
    if(!is.factor(selected_pheno_data)){

        #extract the average of each genotype
        minor_homo_average = round(mean(na.omit(eval(parse(text=paste("subset_minor_homo$", selected_pheno, sep=""))))), number_decimals)
        hetero_average = round(mean(na.omit(eval(parse(text=paste("subset_hetero$", selected_pheno, sep=""))))), number_decimals)
        major_homo_average = round(mean(na.omit(eval(parse(text=paste("subset_major_homo$", selected_pheno, sep=""))))), number_decimals)

        #extract the SD of each genotype
        minor_homo_sd = round(sd(na.omit(eval(parse(text=paste("subset_minor_homo$", selected_pheno, sep=""))))), number_decimals)
        hetero_sd = round(sd(na.omit(eval(parse(text=paste("subset_hetero$", selected_pheno, sep=""))))), number_decimals)
        major_homo_sd = round(sd(na.omit(eval(parse(text=paste("subset_major_homo$", selected_pheno, sep=""))))), number_decimals)
    
        #create a single variable with each average and sd
        minor_homo_pheno = paste(minor_homo_average, "$\\pm$", minor_homo_sd, sep="")
        hetero_pheno = paste(hetero_average, "$\\pm$", hetero_sd, sep="")
        major_homo_pheno = paste(major_homo_average, "$\\pm$", major_homo_sd, sep="")
    } else { #if not, and hence the phenotype is a factor

        #calculate the number of minor homo with obesity and the total number of minor homo with data about obesity status
        cases_minor_homo = nrow(eval(parse(text=paste("subset_minor_homo[which(subset_minor_homo$", selected_pheno, "== 1),]", sep=""))))
        total_minor_homo = nrow(eval(parse(text=paste("subset_minor_homo[which(!is.na(subset_minor_homo$", selected_pheno, ")),]", sep=""))))
            #we need all rows with data for the selected phenotype included in the subset of the genotype
        
        #calculate the number of hetero with obesity and the total number of hetero with data about obesity status
        cases_hetero = nrow(eval(parse(text=paste("subset_hetero[which(subset_hetero$", selected_pheno, "== 1),]", sep=""))))
        total_hetero = nrow(eval(parse(text=paste("subset_hetero[which(!is.na(subset_hetero$", selected_pheno, ")),]", sep=""))))
            #we need all rows with data for the selected phenotype included in the subset of the genotype

        #calculate the number of major homo with obesity and the total number of major homo with data about obesity status
        cases_major_homo = nrow(eval(parse(text=paste("subset_major_homo[which(subset_major_homo$", selected_pheno, "== 1),]", sep=""))))
        total_major_homo = nrow(eval(parse(text=paste("subset_major_homo[which(!is.na(subset_major_homo$", selected_pheno, ")),]", sep=""))))
            #we need all rows with data for the selected phenotype included in the subset of the genotype

        #calculate the percentage of individuals with obesity respect the total number of individual within each genotype
        minor_homo_pheno=round((cases_minor_homo/total_minor_homo)*100, number_decimals)
        hetero_pheno=round((cases_hetero/total_hetero)*100, number_decimals)
        major_homo_pheno=round((cases_major_homo/total_major_homo)*100, number_decimals)

        #I have checked the figure of PTPN1 for obesity (rs2143511 - additive: pdf with FDR<0.1). I find it very strange to me. I calculated the number of individuals with a genotype and without obesity and then divide by the total number of cases without obesity across the three genotypes. I did this for each genotype within obesity and non-obesity. I have checked that the sample sizes and percentages are correct (see annotated lines below), BUT it is very strange way to present these results. You can se how TT (major homo) frequency decreases a lot in overweight, while the minor homo (CC) increases, suggesting the former is protective. But it is very difficult to see because the percentage of individuals with the TC is the highest, but the significant model is additive... It is better to calculate the percentage of obesity within a genotype, in that way you can see a decrease in the percentage of obesity from the minor homo to the major homo.

            #sample size in each category
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$obesity == 0),])

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$obesity == 1),])

            #percentage of individuals in each category

                #total_cases_non_obese = length(which(!is.na(myData_ptpn1$rs2143511) & myData_ptpn1$obesity == 0))
                #total_cases_obese = length(which(!is.na(myData_ptpn1$rs2143511) & myData_ptpn1$obesity == 1))
                    #number of cases with and without obesity for the three genotypes (no NA for the SNP)

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$obesity == 0),])/total_cases_non_obese
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$obesity == 0),])/total_cases_non_obese
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$obesity == 0),])/total_cases_non_obese

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$obesity == 1),])/total_cases_obese
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$obesity == 1),])/total_cases_obese
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$obesity == 1),])/total_cases_obese
    }

    #if the number of selected rows is 2, and hence the spn-phenotype association is significant for additive and codominant model
    if(nrow(selected_rows) == 2){

        #extract the association results of additive model from selected rows
        p_value_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$p_value, 3)
        fdr_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$fdr, 3)
        r2_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$r2_percentage, 3)

        #extract the association results of codominant model from selected rows
        p_value_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$p_value, 3)
        fdr_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$fdr, 3)
        r2_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$r2_percentage, 3)
    } else { #if not, only one heritage model is significant

        #thus we should have only 1 selected row
        if(nrow(selected_rows) == 1){
            
            #save the name of the significant
            model_significant = selected_rows$heritage_model
                #the significant model is in selected row, because this is the association selected as significant
            
            #save the name of the non-significant model
            model_no_significant = ifelse(model_significant == "additive", "codominant", "additive")
                #if the significant model is additive, then the non-significant model is codominant. If not, then the significant model is codominant and hence the non-significant is additive.

            #extract the results of the significant model from selected rows (there you have the significant associations) and then assign these results to objects named with the name of the significant model
            assign(paste("p_value_", model_significant, sep=""), round(selected_rows$p_value, 3))
            assign(paste("fdr_", model_significant, sep=""), round(selected_rows$fdr, 3))
            assign(paste("r2_", model_significant, sep=""), round(selected_rows$r2_percentage, 3))

            #extract the results of the non-significant model from crude_assocs (there you have all associations)
            results_model_no_significant = crude_assocs[which(crude_assocs$pheno_snp_combination == selected_combination & crude_assocs$heritage_model == model_no_significant),]
                #we need the whole supple (non-significant associations), but with the combination pheno-snp selected and the non-significant model

            #then assign these results to objects named with the name of the non-significant model
            assign(paste("p_value_", model_no_significant, sep=""), round(results_model_no_significant$p_value,3))
            assign(paste("fdr_", model_no_significant, sep=""), round(results_model_no_significant$fdr,3))
            assign(paste("r2_", model_no_significant, sep=""), round(results_model_no_significant$r2_percentage,3))
        }
    }

    #extract the final name of the selected phenotype for the table
    pheno_table = pheno_names[which(pheno_names$var_name == selected_pheno),]$final_name

    #bind results into one row
    results = cbind.data.frame(
        selected_snp,
        pheno_table,
        minor_homo_pheno,
        hetero_pheno,
        major_homo_pheno,
        p_value_additive,
        fdr_additive,
        r2_additive,
        p_value_codominant,
        fdr_codominant,
        r2_codominant)

    #change columns names to match names table 2
    names(results)[1] <- "SNP"
    names(results)[2] <- "Phenotype"
    names(results)[3] <- "11"
    names(results)[4] <- "12"
    names(results)[5] <- "22"
    names(results)[6] <- "P add"    
    names(results)[7] <- "FDR add"
    names(results)[8] <- "R2 add"
    names(results)[9] <- "P cod"    
    names(results)[10] <- "FDR cod"
    names(results)[11] <- "R2 cod"

    #add the results as a row in the final data.frame
    table_4 = rbind.data.frame(table_4, results)
}

#remove the first row with all NAs
table_4 = table_4[-which(rowSums(is.na(table_4)) == ncol(table_4)),]


## change some name entries
#add R squared and percentage to the column names with R2
colnames(table_4)[which(colnames(table_4) == "R2 add")] <- "R\\textsuperscript{2} add (\\%)"
colnames(table_4)[which(colnames(table_4) == "R2 cod")] <- "R\\textsuperscript{2} cod (\\%)"


##for the phenotypes with percentage (%) or squared (^2), make some changes to be acceptable for latex.

#We have to add slash (\\) to avoid problems in latex (two because one is an expression for R).

#pheno with slash
pheno_slash = which(grepl("%", table_4$Phenotype)) #rows with percentage as phenotype name
#change names of these phenotypes modifying "%" by "\\%"
table_4[pheno_slash,]$Phenotype <- gsub("%", "\\%", table_4[pheno_slash,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression

#pheno with ^2
pheno_squared = which(grepl("\\^2", table_4$Phenotype)) #rows with squared as phenotype name
#change names of these phenotypes modifying "^2" by "\\textsuperscript{2}"
table_4[pheno_squared,]$Phenotype <- gsub("^2", "\\textsuperscript{2}", table_4[pheno_squared,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression




####################
##### TABLE 5 ######
####################

#This table will show the average values of phenotypes per each genotype*PA level along with the FDR of the additive and codominant model.

#copy the supple data to do some operations
crude_interacts = suppl_data_2
    #we use the supple dataset file because it has the final phenotypes and columns we need


### check that what happens if an crude SNP-phenotype is not significant for one of the models (codominant or additive) ###

#No problem. For the interaction analyses, we considered those phenotypes and SNPs that were previously associated with FDR<0.1, independently of the heritage model. For example, if BMI was associated with rs2143511 under dominant model, we tested the interaction between that SNP and physical activity on BMI across the 5 heritage models. In other words, I the interaction between a SNP and PA on a phenotype is calculated under additive model, it is also calculated under the codominant model.

#Therefore, if a SNP-phenotype association has an FDR lower than 0.1 under additive model, but it is not significant under the codominant model, there is no problem. The script will take the non-significant FDR for the interaction under codominant model from crude_interacts (supple 2), where the results for the interactions under all heritage models are included.

#you can check that there is no problem. Unannotate the next line and you will see how setting as NA the results from codominant model gives NA for all entry of the table for codominant. 
#crude_interacts[which(crude_interacts$heritage_model == "codominant"), c("p_value", "fdr", "r2_percentage")] <- NA
    #select the p_value, fdr and r2 for all interaction results belonging to the codominant model and set them as NA.
    #The results is NA for the p_value, fdr and r2 under codominant model for all snp-phenotype combinations, while the additive model includes all results.


### check if considering additive/codominant models cover all significant associations ###

#create a variable with the combination of phenotype and snp of each association
crude_interacts$pheno_snp_combination = interaction(crude_interacts$phenotype, crude_interacts$snp, sep="-")
#check
summary(crude_interacts$pheno_snp_combination == paste(crude_interacts$phenotype, crude_interacts$snp, sep="-"))

#select those associations with an FDR<0.05
interact_fdr_less_005 = crude_interacts[which(crude_interacts$fdr<0.05),]

#from these associations, select those with the additive and codominant model
interact_fdr_less_005_only_add_cod = crude_interacts[which(crude_interacts$fdr<0.05 & crude_interacts$heritage_model %in% c("additive", "codominant")),]
#check
nrow(interact_fdr_less_005_only_add_cod[which(interact_fdr_less_005_only_add_cod$fdr<0.05 & interact_fdr_less_005_only_add_cod$heritage_model %in% c("additive", "codominant")),]) == nrow(interact_fdr_less_005_only_add_cod)

#check whether all significant combinations are included when restricting to add/codominant
summary(unique(interact_fdr_less_005$pheno_snp_combination) %in% unique(interact_fdr_less_005_only_add_cod$pheno_snp_combination))

#select those not included
interact_fdr_less_005[which(!interact_fdr_less_005$pheno_snp_combination %in% interact_fdr_less_005_only_add_cod$pheno_snp_combination),]
    #There are only 3 interactions not included: 
        #obesity overdominant rs2143511 0.01062509724 0.04250038896
        # FMI    recessive    rs6067472 0.02494014612 0.04500762131
        # FMI    recessive    rs6020608 0.03375571598 0.04500762131
    #All three are very close to FDR=0.05.
        #FMI is very associated with other two snps in other models.
        #Obesity is not associated with other snps or models, but I am not very confident with this association. It only shows significance for the overdominant model, similar phenotype for the homozygous.
    #There is no problem if we lose these associations in the table of main text. The main results and analyses remain similar in main text, and the these three can be found in supplementary 2.


### make a loop to extract phenotype values per genotype ###

#open empty data.frame to save results
table_5 = rbind.data.frame(rep(NA, 14))
colnames(table_5)[1] <- "SNP"
colnames(table_5)[2] <- "Phenotype"
colnames(table_5)[3] <- "11 PI"
colnames(table_5)[4] <- "12 PI"
colnames(table_5)[5] <- "22 PI"
colnames(table_5)[6] <- "11 PA"
colnames(table_5)[7] <- "12 PA"
colnames(table_5)[8] <- "22 PA"
colnames(table_5)[9] <- "P add" 
colnames(table_5)[10] <- "FDR add"
colnames(table_5)[11] <- "R2 add"
colnames(table_5)[12] <- "P cod"
colnames(table_5)[13] <- "FDR cod"
colnames(table_5)[14] <- "R2 cod"

#reorder the file with significant add/cod associations based on snps
interact_fdr_less_005_only_add_cod_copy = interact_fdr_less_005_only_add_cod[order(interact_fdr_less_005_only_add_cod$snp),]

#select those pheno-snp combinations for associations with FDr<0.1 in the additive or codominant models
pheno_snp_combinations_table_5 = unique(interact_fdr_less_005_only_add_cod_copy$pheno_snp_combination)
    #we need unique because an snp-pheno association that is significant for additive and codominant will be selected two times and included two times in the table. We only need it one time, if the combination is present two times in "interact_fdr_less_005_only_add_cod_copy" (2 rows), the script will take the FDR of additive and codominant. If it is not present it will look for the non-significant model in "crude_interacts". 

#set the name of the environmental variable
env_variable_table_5 = "PA_factor"

#for each pheno-snp combination, considering only associations with FDR<0.1 and additive model
for(i in 1:length(pheno_snp_combinations_table_5)){

    #select the [i] combination
    selected_combination = pheno_snp_combinations_table_5[i]

    #check
    print(paste("##############################################"))
    print(paste("STARTING WITH TABLE 5: ", selected_combination, sep=""))
    print(paste("##############################################"))

    #select the rows of the significant associations for the [i] pheno-snp combination
    selected_rows = interact_fdr_less_005_only_add_cod_copy[which(interact_fdr_less_005_only_add_cod_copy$pheno_snp_combination == selected_combination),]
        #We can have the same phenotype-snp combination for two models, additive and codominant.

    #select the [i] phenotype and snp
    selected_pheno = unique(selected_rows$phenotype)
    selected_snp = unique(selected_rows$snp)
        #unique because we can have two rows (additive and codominant), but in both cases the snp and phenotype should be the same indicated in selected_combination
    #check
    print(paste("############################"))
    print(paste("WE SELECTED THE CORRECT PHENOTYPE AND SNP?"))
    print(paste(selected_pheno, "-", selected_snp, sep="") == selected_combination)
    print(paste("############################"))

    #extract the genotype data of the [i] snp
    geno_data = eval(parse(text=paste("na.omit(myData_ptpn1$", selected_snp, ")", sep="")))

    #extract genotype levels
    geno_data_levels = unique(geno_data)

    #for each genotype level calculate the number of individuals
    genotype_sample_size = data.frame(selected_genotype=NA, sample_size=NA) #open an empty data.frame to save resuls
    for(j in 1:length(geno_data_levels)){

        #select the [j] genotype
        selected_genotype = geno_data_levels[j]

        #calculate the number of individuals with the [j] snp
        sample_size = length(which(geno_data == selected_genotype))

        #save the results
        genotype_sample_size = rbind.data.frame(genotype_sample_size, cbind.data.frame(selected_genotype, sample_size))
    }

    #remove the first row with all NAs
    genotype_sample_size = genotype_sample_size[-which(rowSums(is.na(genotype_sample_size)) == ncol(genotype_sample_size)),]

    #select genotypes of the homozygotes
    genotype_sample_size_homo = genotype_sample_size[which(!genotype_sample_size$selected_genotype %in% c("1/2", "2/1")),]

    #select the homozygote genotype with the lowest sample size, that is, minor homozygote
    minor_homo = genotype_sample_size_homo[which(genotype_sample_size_homo$sample_size == min(genotype_sample_size_homo$sample_size)),]$selected_genotype

    #select the homozygote genotype with the highest sample size, that is, major homozygote
    major_homo = genotype_sample_size_homo[which(genotype_sample_size_homo$sample_size == max(genotype_sample_size_homo$sample_size)),]$selected_genotype

    #select those individuals with each type of genotype
    subset_minor_homo = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", selected_snp, " == '", minor_homo, "' ),]", sep="")))
    subset_hetero = eval(parse(text=paste("myData_ptpn1[which(!myData_ptpn1$", selected_snp, " %in% c('", minor_homo, "', '", major_homo, "') & !is.na(myData_ptpn1$", selected_snp, ")),]", sep="")))
    subset_major_homo = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", selected_snp, " == '", major_homo, "' ),]", sep="")))
        #"parse()" returns the parsed but unevaluated expressions in an "expression" object. Therefore, the code is not run yet.
        #"eval()" evaluates an R expression. It runs the expression. 
            #dummy example:
                #expression_sum = parse(text="1+1") #create a expression to sum 1+1
                #eval(expression_sum) #evaluate the expression, giving 2.

    #check the subset worked
    print(paste("############################"))
    print(paste("CHECK THE SUBSET WAS WELL FOR: ", selected_combination, sep=""))
    print(paste("############################"))
    #no other genotype rather than minor homo should exist in subset_minor_homo
    print(nrow(eval(parse(text=paste("subset_minor_homo[which(subset_minor_homo$", selected_snp, " != '", minor_homo, "'),]", sep="")))) == 0)
    #no other genotype rather than hetero should exist in subset_hetero
    print(nrow(eval(parse(text=paste("subset_hetero[which(subset_hetero$", selected_snp, " %in% c('", minor_homo, "', '", major_homo, "') | is.na(subset_hetero$", selected_snp, ")),]", sep="")))) == 0)
    #no other genotype rather than major homo should exist in subset_major_homo
    print(nrow(eval(parse(text=paste("subset_major_homo[which(subset_major_homo$", selected_snp, " != '", major_homo, "'),]", sep="")))) == 0)

    #set the number decimals
    #if the phenotype is
    if(selected_pheno %in% c()){

        #0 decimals
        number_decimals = 0
    } else {#if the phenotype is none of the latter

        #2 decimals
        number_decimals = 2
    }

    #extract the data of the selected phenotype to make a condition
    selected_pheno_data = eval(parse(text=paste("myData_ptpn1$", selected_pheno, sep="")))

    #extract the environmental variable to do a loop inside the ifelse
    selected_env_var_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", env_variable_table_5, sep="")))))

    #if the selected phenotype is not a factor
    if(!is.factor(selected_pheno_data)){

        #for each level of the environmental variable
        for(j in 1:length(selected_env_var_levels)){

            #select the [j] level
            selected_env_level = selected_env_var_levels[j]

            #select those individuals with the [j] level of the environmental variable in each genotype
            subset_minor_homo_env = eval(parse(text=paste("subset_minor_homo[which(subset_minor_homo$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))
            subset_hetero_env = eval(parse(text=paste("subset_hetero[which(subset_hetero$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))
            subset_major_homo_env = eval(parse(text=paste("subset_major_homo[which(subset_major_homo$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))

            #check the subset worked
            print(paste("############################"))
            print(paste("CHECK THE ENV SUBSET WAS WELL FOR: ", selected_combination, sep=""))
            print(paste("############################"))
            #no other genotype rather than minor homo should exist in subset_minor_homo
            print(nrow(eval(parse(text=paste("subset_minor_homo_env[which(subset_minor_homo_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)
            print(nrow(eval(parse(text=paste("subset_hetero_env[which(subset_hetero_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)
            print(nrow(eval(parse(text=paste("subset_major_homo_env[which(subset_major_homo_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)

            #extract the average of each genotype
            minor_homo_average = round(mean(na.omit(eval(parse(text=paste("subset_minor_homo_env$", selected_pheno, sep=""))))), number_decimals)
            hetero_average = round(mean(na.omit(eval(parse(text=paste("subset_hetero_env$", selected_pheno, sep=""))))), number_decimals)
            major_homo_average = round(mean(na.omit(eval(parse(text=paste("subset_major_homo_env$", selected_pheno, sep=""))))), number_decimals)

            #extract the SD of each genotype
            minor_homo_sd = round(sd(na.omit(eval(parse(text=paste("subset_minor_homo_env$", selected_pheno, sep=""))))), number_decimals)
            hetero_sd = round(sd(na.omit(eval(parse(text=paste("subset_hetero_env$", selected_pheno, sep=""))))), number_decimals)
            major_homo_sd = round(sd(na.omit(eval(parse(text=paste("subset_major_homo_env$", selected_pheno, sep=""))))), number_decimals)
    
            #create a single variable with each average and sd and assign it to the corresponding object
            assign(paste("minor_homo_", env_variable_table_5, "_", selected_env_level, sep=""), paste(minor_homo_average, "$\\pm$", minor_homo_sd, sep=""))
            assign(paste("hetero_", env_variable_table_5, "_", selected_env_level, sep=""), paste(hetero_average, "$\\pm$", hetero_sd, sep=""))
            assign(paste("major_homo_", env_variable_table_5, "_", selected_env_level, sep=""), paste(major_homo_average, "$\\pm$", major_homo_sd, sep=""))
        }
    } else { #if not, and hence the phenotype is a factor
        
        #for each level of the environmental variable
        for(j in 1:length(selected_env_var_levels)){

            #select the [j] level
            selected_env_level = selected_env_var_levels[j]

            #select those individuals with the [j] level of the environmental variable in each genotype
            subset_minor_homo_env = eval(parse(text=paste("subset_minor_homo[which(subset_minor_homo$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))
            subset_hetero_env = eval(parse(text=paste("subset_hetero[which(subset_hetero$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))
            subset_major_homo_env = eval(parse(text=paste("subset_major_homo[which(subset_major_homo$", env_variable_table_5, "==", selected_env_level, "),]", sep="")))

            #check the subset worked
            print(paste("############################"))
            print(paste("CHECK THE ENV SUBSET WAS WELL FOR: ", selected_combination, sep=""))
            print(paste("############################"))
            #no other genotype rather than minor homo should exist in subset_minor_homo
            print(nrow(eval(parse(text=paste("subset_minor_homo_env[which(subset_minor_homo_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)
            print(nrow(eval(parse(text=paste("subset_hetero_env[which(subset_hetero_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)
            print(nrow(eval(parse(text=paste("subset_major_homo_env[which(subset_major_homo_env$", env_variable_table_5, " != '", selected_env_level, "'),]", sep="")))) == 0)

            #calculate the number of minor homo with obesity and the total number of minor homo with data about obesity status
            cases_minor_homo = nrow(eval(parse(text=paste("subset_minor_homo_env[which(subset_minor_homo_env$", selected_pheno, "== 1),]", sep=""))))
            total_minor_homo = nrow(eval(parse(text=paste("subset_minor_homo_env[which(!is.na(subset_minor_homo_env$", selected_pheno, ")),]", sep=""))))
                #we need all rows with data for the selected phenotype included in the subset of the genotype

            #calculate the number of hetero with obesity and the total number of hetero with data about obesity status
            cases_hetero = nrow(eval(parse(text=paste("subset_hetero_env[which(subset_hetero_env$", selected_pheno, "== 1),]", sep=""))))
            total_hetero = nrow(eval(parse(text=paste("subset_hetero_env[which(!is.na(subset_hetero_env$", selected_pheno, ")),]", sep=""))))
                #we need all rows with data for the selected phenotype included in the subset of the genotype

            #calculate the number of major homo with obesity and the total number of major homo with data about obesity status
            cases_major_homo = nrow(eval(parse(text=paste("subset_major_homo_env[which(subset_major_homo_env$", selected_pheno, "== 1),]", sep=""))))
            total_major_homo = nrow(eval(parse(text=paste("subset_major_homo_env[which(!is.na(subset_major_homo_env$", selected_pheno, ")),]", sep=""))))
                #we need all rows with data for the selected phenotype included in the subset of the genotype

            #calculate the percentage of individuals with obesity respect the total number of individual within each genotype and assign it to the corresponding object
            assign(paste("minor_homo_", env_variable_table_5, "_", selected_env_level, sep=""), round((cases_minor_homo/total_minor_homo)*100, number_decimals))
            assign(paste("hetero_", env_variable_table_5, "_", selected_env_level, sep=""), round((cases_hetero/total_hetero)*100, number_decimals))
            assign(paste("major_homo_", env_variable_table_5, "_", selected_env_level, sep=""), round((cases_major_homo/total_major_homo)*100, number_decimals))
        }

        #I have checked the figure of PTPN1 for obesity (rs2143511 - overdominant: pdf with FDR<0.05). I find it very strange to me. I calculated the number of individuals with a genotype in given PA and obesity level and then divide by the total number of individuals across all genotypes in that obesity and PA level. I have checked the percentage of the figure and it seems to be calculate in this way (see below), BUT it is very strange way to present these results. You can see how TT-CC (homo) frequency is higher in overweight under low PA compared to the hetero (TC). In high PA, the TT-CC is the less frequent genotype between overweight under high levels of PA, suggesting that TT-CC is protective only under high levels of physical activity. But it is very difficult to see in this way. It is better to calculate the percentage of obesity within a genotype and PA level, in that way you can see a decrease in the percentage of obesity across genotypes and PA levels.

            #percentages calculated in the initial figure

                #low PA, non-obesity
                #homo
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])
                #hetero
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])

                #low PA, obesity
                #homo
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])
                #hetero
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])

                #high PA, non-obesity
                #homo
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])
                #hetero
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])

                #high PA, obesity
                #homo
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])
                #hetero
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 %in% c("1/1", "2/2", "1/2") & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])

            #sample size in each category

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 0),])
                
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 0),])

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])
                
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])

            #percentage of individuals in each category

                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==0 & !is.na(myData_ptpn1$obesity)),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==0 & !is.na(myData_ptpn1$obesity)),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==0 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==0 & !is.na(myData_ptpn1$obesity)),])
                 
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/1" & myData_ptpn1$PA_factor==1 & !is.na(myData_ptpn1$obesity)),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "1/2" & myData_ptpn1$PA_factor==1 & !is.na(myData_ptpn1$obesity)),])
                #nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==1 & myData_ptpn1$obesity == 1),])/nrow(myData_ptpn1[which(myData_ptpn1$rs2143511 == "2/2" & myData_ptpn1$PA_factor==1 & !is.na(myData_ptpn1$obesity)),])
    }

    #if the number of selected rows is 2, and hence the spn-phenotype association is significant for additive and codominant model
    if(nrow(selected_rows) == 2){

        #extract the association results of additive model from selected rows
        p_value_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$p_value, 3)
        fdr_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$fdr, 3)
        r2_additive = round(selected_rows[which(selected_rows$heritage_model == "additive"),]$r2_percentage, 3)

        #extract the association results of codominant model from selected rows
        p_value_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$p_value, 3)
        fdr_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$fdr, 3)
        r2_codominant = round(selected_rows[which(selected_rows$heritage_model == "codominant"),]$r2_percentage, 3)
    } else { #if not, only one heritage model is significant

        #thus we should have only 1 selected row
        if(nrow(selected_rows) == 1){
            
            #save the name of the significant
            model_significant = selected_rows$heritage_model
                #the significant model is in selected row, because this is the association selected as significant
            
            #save the name of the non-significant model
            model_no_significant = ifelse(model_significant == "additive", "codominant", "additive")
                #if the significant model is additive, then the non-significant model is codominant. If not, then the significant model is codominant and hence the non-significant is additive.

            #extract the results of the significant model from selected rows (there you have the significant associations) and then assign these results to objects named with the name of the significant model
            assign(paste("p_value_", model_significant, sep=""), round(selected_rows$p_value, 3))
            assign(paste("fdr_", model_significant, sep=""), round(selected_rows$fdr, 3))
            assign(paste("r2_", model_significant, sep=""), round(selected_rows$r2_percentage, 3))

            #extract the results of the non-significant model from crude_interacts (there you have all associations)
            results_model_no_significant = crude_interacts[which(crude_interacts$pheno_snp_combination == selected_combination & crude_interacts$heritage_model == model_no_significant),]
                #we need the whole supple (non-significant associations), but with the combination pheno-snp selected and the non-significant model

            #then assign these results to objects named with the name of the non-significant model
            assign(paste("p_value_", model_no_significant, sep=""), round(results_model_no_significant$p_value,3))
            assign(paste("fdr_", model_no_significant, sep=""), round(results_model_no_significant$fdr,3))
            assign(paste("r2_", model_no_significant, sep=""), round(results_model_no_significant$r2_percentage,3))
        }
    }

    #extract the final name of the selected phenotype for the table
    pheno_table = pheno_names[which(pheno_names$var_name == selected_pheno),]$final_name

    #bind results into one row
    results = cbind.data.frame(
        selected_snp,
        pheno_table,
        minor_homo_PA_factor_0,
        hetero_PA_factor_0,
        major_homo_PA_factor_0,
        minor_homo_PA_factor_1,
        hetero_PA_factor_1,
        major_homo_PA_factor_1,
        p_value_additive,
        fdr_additive,
        r2_additive,
        p_value_codominant,
        fdr_codominant,
        r2_codominant)

    #change columns names to match names table 2
    names(results)[1] <- "SNP"
    names(results)[2] <- "Phenotype"
    names(results)[3] <- "11 PI"
    names(results)[4] <- "12 PI"
    names(results)[5] <- "22 PI"
    names(results)[6] <- "11 PA"
    names(results)[7] <- "12 PA"
    names(results)[8] <- "22 PA"
    names(results)[9] <- "P add"    
    names(results)[10] <- "FDR add"
    names(results)[11] <- "R2 add"
    names(results)[12] <- "P cod"    
    names(results)[13] <- "FDR cod"
    names(results)[14] <- "R2 cod"

    #add the results as a row in the final data.frame
    table_5 = rbind.data.frame(table_5, results)
}

#remove the first row with all NAs
table_5 = table_5[-which(rowSums(is.na(table_5)) == ncol(table_5)),]


## change some name entries
#add R squared and percentage to the column names with R2
colnames(table_5)[which(colnames(table_5) == "R2 add")] <- "R\\textsuperscript{2} add (\\%)"
colnames(table_5)[which(colnames(table_5) == "R2 cod")] <- "R\\textsuperscript{2} cod (\\%)"


##for the phenotypes with percentage (%) or squared (^2), make some changes to be acceptable for latex.

#We have to add slash (\\) to avoid problems in latex (two because one is an expression for R).

#pheno with slash
pheno_slash = which(grepl("%", table_5$Phenotype)) #rows with percentage as phenotype name
#change names of these phenotypes modifying "%" by "\\%"
table_5[pheno_slash,]$Phenotype <- gsub("%", "\\%", table_5[pheno_slash,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression

#pheno with ^2
pheno_squared = which(grepl("\\^2", table_5$Phenotype)) #rows with squared as phenotype name
#change names of these phenotypes modifying "^2" by "\\textsuperscript{2}"
table_5[pheno_squared,]$Phenotype <- gsub("^2", "\\textsuperscript{2}", table_5[pheno_squared,]$Phenotype, fixed=TRUE)
    #fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression




################################################
##### CONVERT TABLES TO LATEX AND COMPILE ######
################################################

#load the required package
require(xtable)

#path to save tex table
path_tex_table = "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/tables/"

#name tex table
name_tex_table = "tables_latex_v4.tex"

#name doc table
name_doc_table = "tables_latex_v4.odt"

#remove the previous file with tables in latex
system(paste("cd ", path_tex_table, "; rm ", name_tex_table, "; rm ", name_doc_table, sep=""))

#convert table 1 to a latex table
print.xtable(xtable(table_1_1KGP, caption="Table 1", label=NULL, align="llccccc", digits=2, display=c("s", "s", "s", "s", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_1_1KGP)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date())
        #arguments xtable
            #caption: caption of the table
            #label: label of the table
            #align: position of the content within the cell. The number is ncols + 1 because it also consider the row names, even if you remove them in print.xtable. The first position is for NO column, the second one is for the first column, the third is for the second column, etc.
            #digits: number of decimals. You can set negative decimals here to have scientific notation, but this can be also done in the next argument. You can set a different number of decimals per column (ncol+1)
            #display: way to show content. s for string, f for usual numbers xxxx.xxxx and E/e for scientific notation (in upper and lower case, respectively). The first position is for NO column, the second one is for the first column, the third is for the second column, etc.
        #arguments print.xtable
            #type: type of table, latex or html
            #file: file to save the table
            #append: logical to indicate if the table should be appended in the file or that file should be overwritten.
            #floating: logical indicating whether this a floating table of latex
            #table.placement: position of the floating table. Only valid with floating=TRUE. The default is [ht], indicating "here" and "top".
            #caption.placement: position of the caption. 
            #caption.width: width of the column.
            #latex.environments: environment in which the table is embedded. For example begin center.
            #hline.after: place to include bold lines in the table. The default is c(-1,0,nrow(table)) which indicated lines at both sides of the first row and in the bottom of the last row. This DOES NOT work in padoc to odf because of the odf file with the style.
            #NA.string: string to include in the cases with NA
            #include.rownames: logical indication whether to include row.names or not
            #comment: logical indicating whether a comment is included in the table or not.
            #timestamp: timestamp included in case "comment" is TRUE. The default is date().    

#convert table 2 to a latex table
print.xtable(xtable(table_2, caption="Table 2", label=NULL, align="llcccccc", digits=c(0,0,0,2,0,2,0,2), display=c("s", "s", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_2)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #digits: We select the number of digits because we want each column to have a different number of digits. 
        #Numeric vector of length equal to one (in which case it will be replicated as necessary) or to the number of columns of the resulting table plus 1 for the column with the row.names. These row names will be then removed setting FALSE the argument "include.rownames" from print.xtable. 
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...
    #Rest of the argument in the line of the table 1

#convert table 3 to a latex table
print.xtable(xtable(table_3, caption="Table 3", label=NULL, align="llcccccc", digits=2, display=c("s", "s", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_3)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...
    #Rest of the argument in the line of the table 1

#convert table 4 to a latex table
print.xtable(xtable(table_4, caption="Table 4", label=NULL, align="cccccccccccc", digits=3, display=c("s", "s", "s", "f", "f", "f", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_4)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #argument in the line of the table 1
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...

#convert table 5 to a latex table
print.xtable(xtable(table_5, caption="Table 5", label=NULL, align="ccccccccccccccc", digits=3, display=c("s", "s", "s", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_5)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #argument in the line of the table 1
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...

#compile to odt with pandoc
system(paste("cd ", path_tex_table, "; pandoc -s ", name_tex_table, " -o ", name_doc_table, sep=""))

#From the ".odt" file, you have to copy each table to excel, using advanced options for pasting ("text"). In that way, you can get the +- and other symbols correctly pasted. Once you have the first version in excel with the column names correctly displayed (squared, etc...), you can just paste the content (phenotype values, frequencies...).




#################################################
################### FIGURES #####################
#################################################
#Figures are removed in this script version.

#FOR FUTURE FIGURES, IMPORTANT: The SNP "rs6067472" is palindromic, so it has T/A in HELENA and A/T in ncbi (see "alleles" object). If you plot this snp with the current plot_assoc functions, the result will be TT TT TT, because the function changes T to A in the genotype data to set A as the major, but then when changing the minor, that A should be T, we have all A, so the result is all T. In these cases, you should first copy the SNP in a new vector, change AA to XX, then TT to AA, and then XX to AA. In the commit "5a18508cb30745ee26117b0aa5aa87985edcf268" of "analyses_fdr_bh_ptpn1_v2.R" you have code for doing that. 

#The figure of obesity are not really understandable. I calculated the percentage of individuals with each genotype within obesity and non-obesity. This is not intuitive. It is better to calculate the percentage of obesity within each genotype. See the script of table 4 (for factor phenotypes) for further details.




####################################################
######## SUPLE LEPTIN - OBESITY CORRELATION ########
####################################################

if(FALSE){

#open the pdf
pdf(paste("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/figures/figures_suple_pvals/figure_S", sequence_pheno[length(sequence_pheno)]+1, ".pdf", sep=""), height = 6, width = 12)
par(mfrow=c(1,2),  mar=c(6.5, 4, 2, 2) +0.1)

##plot body fat percentage vs. leptin
#make the plot
plot(myData_ptpn1$CRF_Body_fat_PC, myData_ptpn1$Leptin_ng_ml, type="p", xlab="Body fat %", ylab="Leptin (ng/ml)", cex.lab=1.5)

#make the correlation
tests_pc = cor.test(myData_ptpn1$CRF_Body_fat_PC, myData_ptpn1$Leptin_ng_ml, test="spearman")

#extract and plot the results of the correlation
tests_pc_p = bquote(italic(p.value) == .(format(tests_pc$p.value, digits = 3)))
text(x=17, y=160, labels = tests_pc_p, cex=1.3)
tests_pc_t = bquote(italic(t) == .(format(tests_pc$statistic, digits = 3)))
text(x=17, y=145, labels = tests_pc_t, cex=1.3)
tests_pc_rho = bquote(italic(rho) == .(format(tests_pc$estimate, digits = 3)))
text(x=17, y=130, labels = tests_pc_rho, cex=1.3)


##plot FMI vs. leptin
#make the plot
plot(myData_ptpn1$FMI, myData_ptpn1$Leptin_ng_ml, type="p", xlab="Fat mass index (kg/m^2)", ylab="Leptin (ng/ml)", cex.lab=1.5)

#make the correlation
tests_fmi = cor.test(myData_ptpn1$FMI, myData_ptpn1$Leptin_ng_ml, test="spearman")

#extract and plot the results of the correlation
tests_fmi_p = bquote(italic(p.value) == .(format(tests_fmi$p.value, digits = 3)))
text(x=17, y=160, labels = tests_fmi_p, cex=1.3)
tests_fmi_t = bquote(italic(t) == .(format(tests_fmi$statistic, digits = 3)))
text(x=17, y=145, labels = tests_fmi_t, cex=1.3)
tests_fmi_rho = bquote(italic(rho) == .(format(tests_fmi$estimate, digits = 3)))
text(x=17, y=130, labels = tests_fmi_rho, cex=1.3)

#add the title plot
mtext(paste("Online supplementary figure S", sequence_pheno[length(sequence_pheno)]+1, sep=""), side=1, font=2, cex=2, adj=0.015, padj=1.5, outer=TRUE, line=-3)

#close the pdf
dev.off()
}