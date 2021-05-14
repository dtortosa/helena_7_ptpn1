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
    c("Center (%)", "Age (%)", "Weight (kg)", "Height (cm)", "Triceps skinfold (mm)", "Subescapular skinfold (mm)", "% individuals", "BMI (kg/m^2)", "Waist circum. (cm)", "Waist/Height ratio", "Hip circum. (cm)", "Waist/Hip ratio", "Body fat (%)", "FMI (kg/m^2)", "Total cholesterol (mg/dL)","LDL-C (mg/dL)","HDL-C (mg/dL)","Total cholesterol/HDL-C","LDL-C/HDL-C","Triglycerides (mg/dL)","Triglycerides/HDL-C","ApoA1 (mg/dL)","ApoB (mg/dL)","ApoB/ApoA1","ApoB/LDL-C","Insulin (micro lU/mL)","HOMA","QUICKI","Leptin (ng/ml)","SBP (mm Hg)","DBP (mm Hg)"))
colnames(pheno_names) <- c("var_name", "final_name")

#### binomial phenotypes included in this set of results
binomial_pheno = c("obesity", "CVi_BP")



#### the same for SSB variable
ssb_names = cbind.data.frame(
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"),
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"))
colnames(ssb_names) <- c("var_name", "final_name")



##### level names of factor variables
factor_variables_levels = list(c("Non-overweight", "Overweight"), c("Low physic. act. (<60 min/day)", "High physic. act. (≥60 min/day)"))
names(factor_variables_levels) <- c("obesity", "PA_factor")



##### load allele and genes names of SNPs according to HELENA and ncbi
## gene names
gene_names = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/chromosome_snps.csv", sep=",", header=T)
#select snps from the studied gene
gene_names = gene_names[which(gene_names$selected_snp %in% labels(myData_ptpn1)),]
nrow(gene_names) == length(labels(myData_ptpn1))




####################
##### TABLE 1 ######
####################

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



### check minor alleles in HELENA respect to 1000 Genomes Project

## Minor allele frequencies obtained according to 1000 KGP obtained from the NCBI
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


## create a table with these minor allele frequencies
maf_1kgp = rbind.data.frame(NA, 0.3648, 0.0755, 0.4543, 0.2744, 0.4861)
row.names(maf_1kgp) <- c("PTPN1", "rs6067472", "rs10485614", "rs2143511", "rs6020608", "rs968701")
colnames(maf_1kgp) <- c("MAF 1KGP - Europe")


ESTO YA NO HACE FALTA!!! PORQUE LO REVISAMOS EN EL OTRO SCRIPT

## merge with table 1 
table_1_1KGP = merge(table_1, maf_1kgp, by="row.names")


## change the allele that it is different (rs6067472 - T should be the minor, not the major)
table_1_1KGP[which(table_1_1KGP$Row.names == "rs6067472"), which(colnames(table_1_1KGP) == "Major allele")] <- "A"
table_1_1KGP[which(table_1_1KGP$Row.names == "rs6067472"), which(colnames(table_1_1KGP) == "Minor allele")] <- "T"


## reorder the table again following the order in the chromosome 
table_1_1KGP = table_1_1KGP[match(c("PTPN1", as.character(ptpn1_snps$snp)), table_1_1KGP$Row.names),]
table_1_1KGP$Row.names == c("PTPN1", as.character(ptpn1_snps$snp))


## remove the column name for row names
colnames(table_1_1KGP)[which(colnames(table_1_1KGP) == "Row.names")] <- ""


## convert to a latex table
require(xtable)
print.xtable(xtable(table_1_1KGP, align="lcccccc"), include.rownames=FALSE, NA.string="", floating = FALSE)





####################
##### TABLE 2 ######
####################
response_pheno = c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap","obesity","CRF_BMI","CRF_waist","waist_height","CRF_hip","waist_hip","CRF_Body_fat_PC","FMI", "TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","HOMA","QUICKI","Leptin_ng_ml","SBP","DBP")
length(response_pheno) == nrow(pheno_names)
summary(response_pheno == pheno_names$var_name)

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
table_2 = data.frame(cbind(NA, NA, NA, NA, NA, NA, NA))
names(table_2)[1] <- "Phenotype"
names(table_2)[2] <- paste("All non-overweight (n=", n_all_healthy_weight, ")", sep="")
names(table_2)[3] <- paste("All overweight (n=", n_all_overweight, ")", sep="")
names(table_2)[4] <- paste("Male non-overweight (n=", n_male_healthy_weight, ")", sep="")
names(table_2)[5] <- paste("Male overweight(n=", n_male_overweight, ")", sep="")
names(table_2)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
names(table_2)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

#for each phenotype
for (i in 1:length(response_pheno)){

    #select the [i] phenotype
    pheno_selected = response_pheno[i]

    #extract the complete name
    pheno_table = pheno_names$final_name[which(pheno_names$var_name == pheno_selected)]

    #if the phenotype is CVi_BP or obesity (binomial):
    if(pheno_selected %in% c("CVi_BP", "obesity")){

        #extract summary of the variable (sums) across both sexes between overweight and normal weight
        all_sumarize_healthy_weight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))))
        all_sumarize_overweight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))))        

        #calculate the percentage with 1 in all separating by overweight
        all_mean_all_sumarize_healthy_weight = round((all_sumarize_healthy_weight[2]/(all_sumarize_healthy_weight[1]+all_sumarize_healthy_weight[2]))*100,2)
        all_sd_all_sumarize_healthy_weight = NA
        all_mean_all_sumarize_overweight = round((all_sumarize_overweight[2]/(all_sumarize_overweight[1]+all_sumarize_overweight[2]))*100,2)
        all_sd_all_sumarize_overweight = NA

        #extract summary of the variable (sums) across both sexes between overweight and normal weight
        male_sumarize_healthy_weight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))))
        male_sumarize_overweight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))))

        #calculate the percentage with 1 in male separating by overweight
        male_mean_male_sumarize_healthy_weight = round((male_sumarize_healthy_weight[2]/(male_sumarize_healthy_weight[1]+male_sumarize_healthy_weight[2]))*100,2)
        male_sd_male_sumarize_healthy_weight = NA
        male_mean_male_sumarize_overweight = round((male_sumarize_overweight[2]/(male_sumarize_overweight[1]+male_sumarize_overweight[2]))*100,2)
        male_sd_male_sumarize_overweight = NA

        #extract summary of the variable (sums) across both sexes between overweight and normal weight
        female_sumarize_healthy_weight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 0),]$", pheno_selected, sep=""))))))
        female_sumarize_overweight = summary(factor(na.omit(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 1),]$", pheno_selected, sep=""))))))

        #calculate the percentage with 1 in female separating by overweight
        female_mean_female_sumarize_healthy_weight = round((female_sumarize_healthy_weight[2]/(female_sumarize_healthy_weight[1]+female_sumarize_healthy_weight[2]))*100,2)
        female_sd_female_sumarize_healthy_weight = NA
        female_mean_female_sumarize_overweight = round((female_sumarize_overweight[2]/(female_sumarize_overweight[1]+female_sumarize_overweight[2]))*100,2)
        female_sd_female_sumarize_overweight = NA

        #bind results in one row
        results = cbind.data.frame(
            pheno_table,
            paste(all_mean_all_sumarize_healthy_weight, sep=""),
            paste(all_mean_all_sumarize_overweight, sep=""),
            paste(male_mean_male_sumarize_healthy_weight, sep=""),
            paste(male_mean_male_sumarize_overweight, sep=""),
            paste(female_mean_female_sumarize_healthy_weight, sep=""),
            paste(female_mean_female_sumarize_overweight, sep=""))

        #change columns names to match names table 2
        names(results)[1] <- "Phenotype"
        names(results)[2] <- paste("All non-overweight (n=", n_all_healthy_weight, ")", sep="")
        names(results)[3] <- paste("All overweight (n=", n_all_overweight, ")", sep="")
        names(results)[4] <- paste("Male non-overweight (n=", n_male_healthy_weight, ")", sep="")
        names(results)[5] <- paste("Male overweight(n=", n_male_overweight, ")", sep="")
        names(results)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
        names(results)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

        #add to the table
        table_2 = rbind(table_2, results)  
    } else {#if the phenotype is continuous

        #if phenotype selected is a grouping variable like percentage of individuals in two age classes
        if(pheno_selected %in% c("CRF_age", "center")){

            #extract the grouping variable and check if it is a factor
            check_factor = eval(parse(text=paste("myData_ptpn1$", pheno_selected, sep="")))

            #extrac group levels for the group variable (facotr or continuous)
            if(is.factor(check_factor)){

                #set the groups
                group_levels = eval(parse(text=paste("levels(myData_ptpn1$", pheno_selected, ")", sep="")))
                
                #add condition 
                condition_levels = paste("=='", group_levels, "'", sep="")
            } else {

                #set the threshold to divided the continuos variable
                threshold = eval(parse(text=paste("floor(median(myData_ptpn1$", pheno_selected, "))", sep="")))

                #set the groups
                condition_levels = as.character(c(paste(">=", threshold, sep=""), paste("<", threshold, sep="")))
            }

            #open empty data.frame
            percent_calc = data.frame(selected_level=NA, all_sumarize_healthy_weight=NA, all_sumarize_overweight=NA, male_sumarize_healthy_weight=NA, male_sumarize_overweight=NA, female_sumarize_healthy_weight=NA, female_sumarize_overweight=NA)

            #for each level
            for(l in 1:length(condition_levels)){

                #select the [l] level
                selected_level = condition_levels[l]

                #for these subsets we don't use "na.omit" because we include the varialbe of interest in the conditions
                #calculate number of individuals of both sexes for the [l] level
                all_sumarize_healthy_weight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 0 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))
                all_sumarize_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male', 'female') & myData_ptpn1$obesity == 1 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))

                #calculate number of males for the [l] level
                male_sumarize_healthy_weight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 0 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))
                male_sumarize_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('male') & myData_ptpn1$obesity == 1 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))                

                #calculate number of females for the [l] level
                female_sumarize_healthy_weight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 0 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))
                female_sumarize_overweight = nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$CRF_sex %in% c('female') & myData_ptpn1$obesity == 1 & myData_ptpn1$", pheno_selected, selected_level, "),]", sep=""))))  

                #save
                percent_calc = rbind.data.frame(percent_calc, cbind.data.frame(selected_level, all_sumarize_healthy_weight, all_sumarize_overweight, male_sumarize_healthy_weight, male_sumarize_overweight, female_sumarize_healthy_weight, female_sumarize_overweight))
            }
       
            #remove the first row with NAs
            percent_calc = percent_calc[-1,]

            #remove " ' " and the equal ("==") from the selected levels (was use to set the conditions with "which")
            percent_calc$selected_level = gsub("'", "", percent_calc$selected_level, fixed=TRUE)
            percent_calc$selected_level = gsub("==", "", percent_calc$selected_level, fixed=TRUE)#fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression

            #set final names of the levels
            percent_calc$selected_level = paste(pheno_table, ": ", percent_calc$selected_level, sep="")

            #copy to save final percentages
            final_percent = percent_calc

            #calculate percentage of each level
            for(c in 2:ncol(final_percent)){

                #calculat total individuals for the [c] category (sex and overweight)
                total_indv = sum(final_percent[,c])

                #calculate and save percentages
                final_percent[,c] = round((final_percent[,c]/total_indv)*100, 2)
            }       

            #change columns names to match names table 2
            names(final_percent)[1] <- "Phenotype"
            names(final_percent)[2] <- paste("All non-overweight (n=", n_all_healthy_weight, ")", sep="")
            names(final_percent)[3] <- paste("All overweight (n=", n_all_overweight, ")", sep="")
            names(final_percent)[4] <- paste("Male non-overweight (n=", n_male_healthy_weight, ")", sep="")
            names(final_percent)[5] <- paste("Male overweight(n=", n_male_overweight, ")", sep="")
            names(final_percent)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
            names(final_percent)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

            #bind to the table
            table_2 = rbind(table_2, final_percent)
        }else{

            #set the number decimals
            #if the phenotype is Age, weight, height, triceps and subscapular fold and BMI
            if(pheno_selected %in% c("CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap", "CRF_BMI", "CRF_waist", "CRF_hip", "FMI")){

                #1 decimal
                number_decimals = 1
            } else {#if not

                #and the phenotype is LDL, HDL, TG, Insulin, Leptine, SBP, DBP
                if(pheno_selected %in% c("LDL", "HDL", "TG", "Insulin", "Leptin_ng_ml", "SBP", "DBP")){

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
            names(results)[5] <- paste("Male overweight(n=", n_male_overweight, ")", sep="")
            names(results)[6] <- paste("Female non-overweight (n=", n_female_healthy_weight, ")", sep="")
            names(results)[7] <- paste("Female overweight (n=", n_female_overweight, ")", sep="")

            #bind to the table
            table_2 = rbind(table_2, results)
        }            
    }
}

#remove the first row
table_2 = table_2[-1,]

#remove the rows of obesity phenotype ("% individuals"; we now have columns separated by overweight) and two centers for which we don't have data: Modena and Birmingham
table_2 = table_2[-which(table_2$Phenotype %in% c("Center (%): Birmingham* in UK", "Center (%): Modena (Italy)", "% individuals")),]

#see the table
table_2

#for the phenotupes with percentahe (%) in the name, we add slash (\\) to avoid problemas i latex (two because 1 is a en expression for R).
#pheno with slash
pheno_slash = which(grepl("%", table_2$Phenotype))
#change names of these phenotypes modifying "%" by "\\%"
table_2[pheno_slash,]$Phenotype <- gsub("%", "\\%", table_2[pheno_slash,]$Phenotype, fixed=TRUE)#fixed: logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as is.  Overrides all conflicting arguments. fixed=TRUE prevents R from using regular expressions, which allow more flexible pattern matching but take time to compute. Without fixed=TRUE, gsub recognise \\ as a regular expression
table_2$Phenotype

#convert to a latex table
require(xtable)
print.xtable(xtable(table_2, align="lccccccc"), include.rownames=FALSE, NA.string="", floating = FALSE, sanitize.text.function=function(x) {x})

#Copia la tabla en tables_latex.tex y corre este comando
system("cd /media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/tables; pandoc -s tables_latex_v3.tex -o tables_latex_v3.odt")
    
#check for all phenotypes there is at least one data
#select phenotype columns    
only_pheno = myData_ptpn1[,which(!colnames(myData_ptpn1)%in%labels(myData_ptpn1))]
#which rows (individuals) have NA for all columns (phenotypes)
length(which(rowSums(is.na(only_pheno)) == ncol(only_pheno))) == 0 #it should be zero
    



####################
##### TABLE 4 ######
####################

#haz un filtro para comparar que aditivo/codominign pilla la mayoría de asociaciones significativas

#hay que cambiar el nombre del major/minor del SNP chungo en la tabla de snp helena names


#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2903808/
#table 2
#11 12  22
# genotyping success rate; Major allele: 1; Minor allele: 2; HW: P-value for Hardy-Weinberg equilibrium. Data are n (frequency).

#fenotipo, snp, media and sd de cada genotipo, FDR the moedl add y co, R2 de ambos modelos. 


#################################################
################### FIGURES #####################
#################################################
#Figures are removed in this script version.




####################################################
######## SUPLE LEPTIN - OBESITY CORRELATION ########
####################################################

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