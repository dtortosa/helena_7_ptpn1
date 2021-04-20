#HERE I CHANGE THE CODE FOR TABLE 2 in LPL code TO INCLUDE SOME SUGGESTIONS FOR REVIEWERS OF CNTF

#################################
##### LOAD ENVIRONMENT ##########
#################################
load("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/rdata/analysis.RData")
require(SNPassoc)
require(genetics)

###set wd
setwd("/Users/diegosalazar/My Drive/science/open_projects/helena_7")

#For this paper we will only use adiposity markers as CVD risk factors only show significant associations with FDR<0.1 for leptin and SBP and three snps (FDR>0.07 in all cases). No haplotype association nor interaction with physical activity.

##### create a data.frame withe the names of phenotypes in dataset and the real name por the figure
pheno_names = cbind.data.frame(
    c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap","obesity","CRF_BMI","CRF_waist","waist_height","CRF_hip","waist_hip","CRF_Body_fat_PC","FMI", "TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","HOMA","QUICKI","Leptin_ng_ml","SBP","DBP"),
    c("Center (%)", "Age (%)", "Weight (kg)", "Height (cm)", "Triceps skinfold (mm)", "Subescapular skinfold (mm)", "% individuals", "BMI (kg/m^2)", "Waist circum. (cm)", "Waist/Height ratio", "Hip circum. (cm)", "Waist/Hip ratio", "Body fat (%)", "FMI (kg/m^2)", "Total cholesterol (mg/dL)","LDL-C (mg/dL)","HDL-C (mg/dL)","Total cholesterol/HDL-C","LDL-C/HDL-C","Triglycerides (mg/dL)","Triglycerides/HDL-C","ApoA1 (mg/dL)","ApoB (mg/dL)","ApoB/ApoA1","ApoB/LDL-C","Insulin (micro lU/mL)","HOMA","QUICKI","Leptin (ng/ml)","SBP (mm Hg)","DBP (mm Hg)"))
colnames(pheno_names) <- c("var_name", "final_name")

#binomial phenotypes included in this set of results
binomial_pheno = c("obesity", "CVi_BP")

##### the same for SSB variable
ssb_names = cbind.data.frame(
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"),
    c("mean_Alspac_v42","CVi_softdrink_cont_2000"))
colnames(ssb_names) <- c("var_name", "final_name")

##### level names of factor variables
factor_variables_levels = list(c("Non-overweight", "Overweight"),
    c("Low physic. act. (<60 min/day)", "High physic. act. (≥60 min/day)"))
names(factor_variables_levels) <- c("obesity", "PA_factor")

##### load allele and genesnames of SNPs according to helena and ncbi
gene_names = read.table("/Users/diegosalazar/My Drive/science/open_projects/helena_7/data/snps/chromosome_snps.csv", sep=",", header=T)
#select snps from the studied gene
gene_names = gene_names[which(gene_names$selected_snp %in% labels(myData_ptpn1)),]

#load allele names
alleles = read.table("/Users/diegosalazar/My Drive/science/open_projects/helena_7/data/snps/alleles_ptpn1.csv", sep=",", header=T)
nrow(gene_names) == length(labels(myData_ptpn1))
nrow(alleles) == length(labels(myData_ptpn1))
    #IMPORTANT NOTE: I have matched the allele names of HELENA and ncbi (at 10/09/2019). Now, they all matched, but I have noted that UCP alleles have changed in ncbi. For example, alleles that I changed to match ncbi, now are in ncbi exactly as original HELENA. If you check this for these SNPs, and you see changes no panic! The important thing is you are using the complementary chain and the first allele is always the major. Indeed, in many cases both options (i.e., both chains) are included as synonimous in ncbi (i.e., HGVS).

####################
##### TABLE 1 ######
####################

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

## add alleles to table 1 
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

#extract the allele names frmo the table1 (combined major and minor) along with snp names
alleles_from_table_1 = cbind.data.frame(row.names(table_1), paste(table_1[,which(colnames(table_1) == "Major allele")], "/", table_1[,which(colnames(table_1) == "Minor allele")], sep=""))
colnames(alleles_from_table_1) <- c("snp", "alleles_from_ncbi")

#merge these alleles names with the original allele names from alleles df
df_check_alleles = merge(alleles_from_table_1, alleles)

#check that the colums of alleles names are similar
identical(df_check_alleles$alleles_from_ncbi, df_check_alleles$ncbi)

#reorder the table following the order in the chromosome 
table_1 = table_1[match(ptpn1_snps$snp, row.names(table_1)),]
row.names(table_1) == ptpn1_snps$snp

#loop for add the gene name to the table
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

#convert to a latex table
require(xtable)
print.xtable(xtable(table_1, align="lcccc"), include.rownames=TRUE, NA.string="", floating = FALSE)


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
system("cd /Users/diegosalazar/My\\ Drive/science/open_projects/helena_7/results/tables; pandoc -s tables_latex_v2.tex -o tables_latex_v2.doc")
    
#check for all phenotypes there is at least one data
#select phenotype columns    
only_pheno = myData_ptpn1[,which(!colnames(myData_ptpn1)%in%labels(myData_ptpn1))]
#which rows (individuals) have NA for all columns (phenotypes)
length(which(rowSums(is.na(only_pheno)) == ncol(only_pheno))) == 0 #it should be zero
    

####################
###### FIGURES #####
####################

##################
#### Figure 1 ####
##################
#HAPLOVIEW. IN PTPN1 this figure was included in the supplementary as first figure.


######################################
### Function to plot single assocs ###
######################################
plot_assoc = function(pheno_selected, selected_model, selected_snp, display_legend, print_model=FALSE, FDR_value = NULL, title_plot=FALSE, title_plot_below=FALSE, output=NULL){

    #SE function to plot erro bars in bar plots
    se <- function(x) sd(x)/sqrt(length(x)) 

    #for discrete phenotypes (obesity) 
    if(pheno_selected %in% c("CVi_BP", "obesity")){
        family_error = "binomial"
        response = pheno_selected
    } else{ #for continous phenotypes
        family_error = "gaussian"
        response = paste("log(", pheno_selected, ")", sep="") #apply log transformation      
    }

    #select the row of results with the interest snp
    selected_row = geno_pheno_results[which(geno_pheno_results$selected_pheno==pheno_selected &
        geno_pheno_results$selected_model==selected_model &
        geno_pheno_results$snp_to_test==selected_snp),]

    #extract p.val asociation
    pval_association = round(selected_row$pvals,5)
    FDR_association = round(selected_row$fdr,4)

    #extract allele names
    allele_names = as.character(alleles[which(alleles$snp == selected_snp),]$helena)
    allele_names = gsub(" ", "", allele_names) #drop the space between characters
    list_allels = strsplit(allele_names, split="/")[[1]]

    #extract genotype names directly from the database using SNPassoc
    if(selected_model == "additive"){
        geno_levels = levels(eval(parse(text=paste("codominant(myData_ptpn1$",selected_snp, ")", sep=""))))

    } else { #if the heritage model is not additive

        #extract genotype names
        geno_levels = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$",selected_snp, ")", sep=""))))  
    }

    #always A is the first, because number are afabetically ordered in genotype levels in Helena
    if("A" %in% list_allels){
        geno_levels = gsub("1", "A", geno_levels)
        geno_levels = gsub("2", list_allels[!"A" == list_allels], geno_levels)                                  
    } else { #if not, C first
        if("C" %in% list_allels){ 
            geno_levels = gsub("1", "C", geno_levels)
            geno_levels = gsub("2", list_allels[!"C" == list_allels], geno_levels)                                                  
        } else { #if not G first (T always will be the last)
            if("G" %in% list_allels){
                geno_levels = gsub("1", "G", geno_levels)
                geno_levels = gsub("2", list_allels[!"G" == list_allels], geno_levels)               
            }
        }            
    }
    
    #drop the slash ("/")
    geno_levels = gsub("/", "", geno_levels)
      
    #change helena names by ncbi names (see "https://mail.google.com/mail/u/0/#search/ruizj%40ugr.es/15ee2eeabdc53653" for further details)
    #load ncbi allele names
    allele_names_ncbi = as.character(alleles[which(alleles$snp == selected_snp),]$ncbi)
    allele_names_ncbi_major = strsplit(allele_names_ncbi, split="/")[[1]][1]
    allele_names_ncbi_major = gsub(" ", "", allele_names_ncbi_major)        
    allele_names_ncbi_minor = strsplit(allele_names_ncbi, split="/")[[1]][2]
    allele_names_ncbi_minor = gsub(" ", "", allele_names_ncbi_minor)
    #load helena allele names
    allele_names_helena = as.character(alleles[which(alleles$snp == selected_snp),]$helena)
    allele_names_helena_major = strsplit(allele_names_helena, split="/")[[1]][1]
    allele_names_helena_major = gsub(" ", "", allele_names_helena_major)        
    allele_names_helena_minor = strsplit(allele_names_helena, split="/")[[1]][2]
    allele_names_helena_minor = gsub(" ", "", allele_names_helena_minor)

    #if the allele name in helena for both minor and major, we change it to the complementary base (A-T, C-G). If one allele is different, the other will be also, because I only change by hand helena alleles when is clear that ncbi uses the complementary. If the major allele is A in helena and T in ncbi, and the minor is C in helena, then for sure I stablished the minor of ncbi as G (This is done BY HAND).
    if(!allele_names_ncbi_major == allele_names_helena_major){
        if(allele_names_ncbi_major == "T"){
            geno_levels = gsub("A", "T", geno_levels)
        }    
        if(allele_names_ncbi_major == "C"){
            geno_levels = gsub("G", "C", geno_levels)                
        }
        if(allele_names_ncbi_major == "G"){
            geno_levels = gsub("C", "G", geno_levels)                
        }   
        if(allele_names_ncbi_major == "A"){
            geno_levels = gsub("T", "A", geno_levels)                
        }                 
    }
    if(!allele_names_ncbi_minor == allele_names_helena_minor){
        if(allele_names_ncbi_minor == "T"){
            geno_levels = gsub("A", "T", geno_levels)
        }
        if(allele_names_ncbi_minor == "C"){
            geno_levels = gsub("G", "C", geno_levels)                
        }
        if(allele_names_ncbi_minor == "G"){
            geno_levels = gsub("C", "G", geno_levels)                
        }    
        if(allele_names_ncbi_minor == "A"){
            geno_levels = gsub("T", "A", geno_levels)                
        }             
    }
    
    #if the phenotype is not a factor, thus the family error is not binomial
    if(!family_error=="binomial"){
        #mean values for all genotype/PA_category combination. 
        mean_values = aggregate(as.formula(paste(pheno_selected, "~", selected_model, "(", selected_snp, ")", sep="")), data = myData_ptpn1, FUN=mean) 

        #mean values for all genotype/PA_category combination. 
        se_values = aggregate(as.formula(paste(pheno_selected, "~", selected_model, "(", selected_snp, ")", sep="")), data = myData_ptpn1, FUN=se) 

        #juntamos mean y se                   
        data_mean<-cbind(mean_values, se_values[,2])  
        
        #change columns names
        colnames(data_mean)[1] <- "selected_snp" 
        colnames(data_mean)[2] <- "mean"
        colnames(data_mean)[3] <- "se" 
        data_mean
                           
        #Creamos una matriz a #partir de los datos.
        xtabs.var<-xtabs(mean ~factor(selected_snp), data=data_mean)
        xtabs.var

        #set the plot limits depending on whether the variable has negative values or not
        if(!pheno_selected == "risk_score"){#if the variable is not the risk score and hence does not have negative values
            #max ylim of bar plot
            max_y = max(data_mean$mean)+ abs(((max(data_mean$mean)*5)/100))

            #max ylim of bar plot
            min_y = min(data_mean$mean)- abs(((min(data_mean$mean)*18)/100))
        } else {#if not
            #max ylim of bar plot
            max_y = max(data_mean$mean)+ abs(((max(data_mean$mean)*400)/100))

            #max ylim of bar plot
            min_y = min(data_mean$mean)- abs(((min(data_mean$mean)*85)/100))
        }
        
        #make the barplot
        if(display_legend==TRUE){ #if there is legend we need different colors for each bar genotype
            xs<-barplot(xtabs.var, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, col=gray.colors(length(unique(data_mean$selected_snp))), lwd=1:2, ylim=c(min_y, max_y)) 
        } else {
            xs<-barplot(xtabs.var, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, lwd=1:2, ylim=c(min_y, max_y)) 
        }
        
        #add line at zero if required        
        if(pheno_selected == "risk_score"){#if the variable is the risk score and hence has negative values add line at zero       
            abline(h=0, lty=2)
        } 

        #add x axis
        if(length(geno_levels) == 2){ #if the model include 2 genotypes
            axis(side=1, at=c(-10, xs[1,1], xs[2,1], 40), labels=c("", bquote(.(geno_levels[1])), bquote(.(geno_levels[2])), "") , cex.lab=0.9, cex.axis=1.2, lwd=1.5, lwd.ticks=0, font=1, mgp=c(3, 1, 0), padj=1) 
        } else { #if the model include 3 genotypes
            axis(side=1, at=c(-10, xs[1,1], xs[2,1], xs[3,1], 40), labels=c("", bquote(.(geno_levels[1])), bquote(.(geno_levels[2])), bquote(.(geno_levels[3])), "") , cex.lab=0.9, cex.axis=1.2, lwd=1, lwd.ticks=0, font=1, mgp=c(3, 1, 0), padj=1)         
        }

        #add x title
        if(print_model==TRUE){
            mtext(text=paste( toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), " ", selected_snp, " (", selected_model, ")", sep=""), side = 1, line=3, outer=FALSE, cex=0.8, font=1)  
        } else {
            mtext(text=paste( toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), " ", selected_snp, sep=""), side = 1, line=3, outer=FALSE, cex=0.8, font=1)         
        }

        #add y axis
        if(min_y < 10){ #if y axis has low values (lower than 10) number used to calculate the tick positions (by argument in seq function) will be smaller
            axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,2)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
        } else {
            axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,1)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
        }

        #add y axis title
        y_title = pheno_names[which(pheno_names$var_name == pheno_selected),]$final_name
        mtext(text=y_title, side = 2, line=2.8, outer=FALSE, cex=0.9, font=1, adj=0.7)

        #add p.values
        if(length(FDR_value)==0){
            mtext(text=paste("P=", pval_association, "; FDR=",  FDR_association, sep=""), side = 3, line=0.45, outer=FALSE, cex=0.8, font=1, adj=0.48)
        } else {
            mtext(text=paste("P=", pval_association, "; FDR=",  FDR_value, sep=""), side = 3, line=0.45, outer=FALSE, cex=0.8, font=1, adj=0.48)
        }

        #calculate n observations for each category
        if(selected_model == "additive"){

            #empty vector
            n_observations = NULL

            #extract snp levels
            levels_snp = levels(eval(parse(text=paste("na.omit(codominant(myData_ptpn1$", selected_snp, "))", sep=""))))

            #loop for extracting nº observations
            for(l in 1:length(levels_snp)){

                #selected snp
                selected_level = levels_snp[l]

                #n obs
                n_observations = append(n_observations, length(eval(parse(text=paste("na.omit(myData_ptpn1[which(codominant(myData_ptpn1$", selected_snp, ") == '", selected_level, "'),]$", pheno_selected, ")", sep="")))))            
            }        
        } else {

            #empty vector
            n_observations = NULL 

            #extract snp levels               
            levels_snp = levels(eval(parse(text=paste("na.omit(", selected_model, "(myData_ptpn1$", selected_snp, "))", sep=""))))
            
            #loop for extracting nº observations
            for(l in 1:length(levels_snp)){
                
                #selected snp                
                selected_level = levels_snp[l]
                
                #n obs                
                n_observations = append(n_observations, length(eval(parse(text=paste("na.omit(myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ") == '", selected_level, "'),]$", pheno_selected, ")", sep="")))))            
            }  
        }

        #plot the number of observations
        for(j in 1:nrow(xs)){ #for each bar 

            #select the bar
            selected_xs = xs[j,]
            
            #plot the correponding n
            mtext(paste("(n = ", n_observations[j], ")", sep=""), side = 1, at = selected_xs, line=0.4, outer=FALSE, cex=0.8)        
        }

        #add a legend
        if(display_legend==TRUE){
            legend("topright", legend=geno_levels, fill =gray.colors(nrow(xs)), bty="n", title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gen), selected_snp), cex=0.85)
        }

        #plot error bars
        for(l in 1:nrow(data_mean)){

            #for each genotype
            subset_data_error_bars = data_mean[l,]
            
            #add the arrow
            arrows(xs[l,1],subset_data_error_bars$mean-subset_data_error_bars$se, xs[l,1], subset_data_error_bars$mean+ subset_data_error_bars$se,angle=90,code=3,length=0.02, lwd=1, xpd=TRUE)  
        }
    } else{#if not, and the phenotype is a factor
        #extract genotype names
        if(selected_model == "additive"){
            geno_levels_number = sort(unique(eval(parse(text=paste("as.factor(na.omit(additive(myData_ptpn1$", selected_snp, ")))", sep="")))))
        } else { #if the heritage model is not additive

            #extract genotype names
            geno_levels_number = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$", selected_snp, ")", sep=""))))  
        }

        #extract levels of the phenotype
        pheno_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", pheno_selected, sep="")))))

        #extract label names of the pheno factor
        unique_pheno_labels = factor_variables_levels[which(names(factor_variables_levels) == pheno_selected)]

        #calculate the height of the bars
        height = NULL
        response_labels = NULL
        geno_labels = NULL
        for(l in 1:length(geno_levels_number)){ #for each genotype

            #select the [l] genotype
            selected_level_l = geno_levels_number[l]

            #for each phenotype level
            for(p in 1:length(pheno_levels)){

                #select the [p] level
                selected_level_p = pheno_levels[p]

                #extract number of individuals with the [l] genotype and the [p] phenotype
                geno_subset = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")== '", selected_level_l, "' & myData_ptpn1$", pheno_selected, "==", selected_level_p, "),]", sep="")))

                #extract number of individuals with the [p] phenotype y without NA for the [l] genotype
                full_geno = eval(parse(text=paste("myData_ptpn1[which(!is.na(", selected_model, "(myData_ptpn1$", selected_snp, ")) & myData_ptpn1$", pheno_selected, "=='", selected_level_p, "'),]", sep="")))

                #calculate the percentage of individuals with the [l] genotype in the group of the [p] phenotype
                percent = round((nrow(geno_subset)/nrow(full_geno))*100, 2)
                
                #save that value as height of the bar
                height = append(height, percent) 

                #save the label of the response variable
                response_labels = append(response_labels, unique_pheno_labels[[1]][p])

                #save the genotype
                geno_labels = append(geno_labels, selected_level_l)               
            }
        }
                
        #convert geno_labels to factor to select exactly the order of the levels and thus the order of the columns
        geno_labels = factor(geno_labels, levels=unique(geno_labels))

        #save all of them into data
        data = cbind.data.frame(response_labels, geno_labels, height)
        
        #create a matrix from the data
        xtabs.var<-xtabs(height ~ geno_labels + response_labels, data=data)
        xtabs.var

        #max ylim of bar plot
        max_y = max(data$height)+ ((max(data$height)*30)/100)

        #max ylim of bar plot
        min_y = min(data$height)- ((min(data$height)*13)/100) 

        #make the barplot
        xs<-barplot(xtabs.var, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, width=3, col=gray.colors(length(unique(geno_levels))), lwd=1:2, ylim=c(min_y, max_y))
       
        #add ylab
        y_lab = as.vector(pheno_names[which(pheno_names$var_name == pheno_selected),]$final_name)
        mtext(text=expression(""), side = 3, line=5.3, outer=FALSE, cex=1)
        mtext(text=paste(y_lab), side = 2, line=2.29, outer=FALSE, cex=0.9, font=1, adj=0.7)

        #add p.values
        if(length(FDR_value)==0){
            mtext(text=paste("P=", pval_association, "; FDR=",  FDR_association, sep=""), side = 3, line=0.45, outer=FALSE, cex=0.8, font=1, adj=0.48)
        } else {
            mtext(text=paste("P=", pval_association, "; FDR=",  FDR_value, sep=""), side = 3, line=0.45, outer=FALSE, cex=0.8, font=1, adj=0.48)
        }

        #add y axis
        if(min_y < 10){ #if y axis has low values (lower than 10) number used to calculate the tick positions (by argument in seq function) will be smaller
            axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,2)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
        } else {
            axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,1)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
        }

        #add x axis
        if(selected_model %in% c("additive", "codominant")){
            axis(side=1, at=c(-10, xs[2,1], xs[2,2], 40), labels=c("", "", "", ""), cex.lab=4, cex.axis=1, lwd=1, lwd.ticks=0, font=1, mgp=c(3, 1, 0)) 
            mtext(side = 1, at=c(-10, xs[2,1], xs[2,2], 40), text = c("", unique_pheno_labels[[1]][1], unique_pheno_labels[[1]][2], ""), line = 1.7, cex=0.8)
        } else {
            axis(side=1, at=c(-10, (xs[1,1]+xs[2,1])/2, (xs[1,2]+xs[2,2])/2, 40), labels=c("", "", "", ""), cex.lab=2, cex.axis=1, lwd=1, lwd.ticks=0, font=1, mgp=c(3, 1, 0))
            mtext(side = 1, at=c(-10, (xs[1,1]+xs[2,1])/2, (xs[1,2]+xs[2,2])/2, 40), text = c("", unique_pheno_labels[[1]][1], unique_pheno_labels[[1]][2], ""), line = 1.8, cex=0.8)                                     
        }    

        #add a legend
        legend("topright", legend=geno_levels, fill =gray.colors(nrow(xs)), selected_snp, bty="n")
   
        #add model and gene name
        if(print_model==TRUE){
            mtext(text=paste( toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), " ", selected_snp, " (", selected_model, ")", sep=""), side = 1, line=3, outer=FALSE, cex=0.8, font=1)  
        } else {
            mtext(text=paste( toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), " ", selected_snp, sep=""), side = 1, line=3, outer=FALSE, cex=0.8, font=1)         
        } 

        #calculate n observations for each category of the factor
        # pheno_levels and geno_levels_number were calculated before
        #for each phenotype level. 
        n_observations = NULL
        for(p in 1:length(pheno_levels)){

            #select the [p] level
            selected_level_p = pheno_levels[p]                

            #for each genotype
            for(i in 1:length(geno_levels_number)){

                #selected level
                seletec_level = geno_levels_number[i]
        
                #no_obs with first level of PA factor
                n_observations = append(n_observations, nrow(eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")=='", seletec_level, "' & myData_ptpn1$", pheno_selected, "=='", selected_level_p, "'),]", sep="")))))
            }                           
        }

        #plot the number of observations
        for(i in 1:length(as.vector(xs))){ #for each bar 

            #select the bar
            selected_xs = as.vector(xs)[i]
        
            #plot the correponding n
            mtext(paste("(", n_observations[i], ")", sep=""), side = 1, at = selected_xs, line=0.4, outer=FALSE, cex=0.7)        
        }
    }

    #if we want to include a title for the page
    if(title_plot_below == TRUE){
        
        #set number figure
        number_figure = strsplit(output, split="_")[[1]][2]

        #insert figure name
        mtext(paste("Online supplementary figure ", number_figure, sep=""), side=1, font=2, cex=1.5, adj=0.015, padj=10.75, outer=TRUE, line=-41)
    }

    #if we want to include a title for the page more below
    if(title_plot == TRUE){
        
        #set number figure
        number_figure = strsplit(output, split="_")[[1]][2]

        #insert figure name
        mtext(paste("Online supplementary figure ", number_figure, sep=""), side=1, font=2, cex=1.5, adj=0.015, padj=1.5, outer=TRUE, line=-3.3)
    }
}


##################
#### Figure 2 ####
##################

#significance threshold for barplots in the main text
significance_threshold_main_text = "fdr<0.1"

#rows to select according to the significance threshold
rows_to_select = eval(parse(text=paste("which(geno_pheno_results$", significance_threshold_main_text, ")", sep="")))

#extract those rows from the crude associations SNP - phenotype
assocs_significant = geno_pheno_results[rows_to_select,]

#check that all associaciones included matches the condition
length(eval(parse(text=paste("which(assocs_significant$", significance_threshold_main_text, ")", sep="")))) == nrow(assocs_significant)
length(eval(parse(text=paste("which(assocs_significant$", significance_threshold_main_text, ")", sep="")))) == nrow(eval(parse(text=paste("geno_pheno_results[which(geno_pheno_results$", significance_threshold_main_text, "),]", sep=""))))
summary(eval(parse(text=paste("geno_pheno_results[-rows_to_select,]$", significance_threshold_main_text, sep=""))))#all false, no case matches the condition after removing those included in assocs_significant

#select only those phenotypes that I decided to include in the manuscript
assocs_significant = assocs_significant[which(assocs_significant$selected_pheno %in% pheno_names$var_name),]
#check
summary(unique(assocs_significant$selected_pheno) %in% pheno_names[-which(pheno_names$var_name %in% c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap")),]$var_name)#all TRUE

#reorder now rows to get bind those rows with the same phenotype
assocs_significant = assocs_significant[order(assocs_significant$selected_pheno, assocs_significant$snp_to_test, assocs_significant$selected_model),]

#If we have binomial phenotypes: reorder assocs_significant to put together interactions of each phenotype. This will be the order used for the interactions. We order interactions in basis on phenotyp to have those cases with obesity or other binomial phenotype bind, because interactions with obesity have two panels each one.
#if we have any of the binomial phenotypes within significant phenotypes
if(TRUE %in% (binomial_pheno %in% assocs_significant$selected_pheno)){

    #extract the row number with the binomial phneotypes and concatenate them after that the rest of row numbers
    new_row_order = c(which(assocs_significant$selected_pheno %in% binomial_pheno), which(!assocs_significant$selected_pheno %in% binomial_pheno))
    print(length(new_row_order) == nrow(assocs_significant))

    #reorder putting binomial phenotypes first
    assocs_significant = assocs_significant[new_row_order,]
}

#number of pages
ceiling(nrow(assocs_significant)/12)

#plot
pdf("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/figures/figures_significant_main_text.pdf", height = 10, width = 7)
par(mfrow=c(6,2), mar=c(5, 4, 2, 2) +0.1)
for(i in 1:nrow(assocs_significant)){

    #select the [i] row
    selected_row = assocs_significant[i,]

    #run the function to plot
    plot_assoc(pheno_selected=selected_row$selected_pheno, selected_model=selected_row$selected_model, selected_snp=selected_row$snp_to_test, FDR_value = NULL, display_legend=FALSE, print_model=TRUE, title_plot=FALSE, title_plot_below=FALSE, output=NULL)
}
dev.off()


##############################
#### Barplots suplementary ###
##############################

#if the significance threshold for barplots in the main text is not FDR<0.1 and hence some of the barplots accepted in the exploratory screening are not showed in the main text
if(!significance_threshold_main_text == "fdr<0.1"){

    ####plot associations significant according to FDR<0.1
    assocs_fdr_0.1 = geno_pheno_results[which(geno_pheno_results$fdr<0.1),]
    #check that all association with FDR<0.1 are included and the rest have FDR>0.1
    length(which(assocs_fdr_0.1$fdr<0.1)) == nrow(assocs_fdr_0.1)
    summary(geno_pheno_results[-which(geno_pheno_results$fdr<0.1),]$fdr>0.1)
    
    #select only those phenotypes that I decided to include in the manuscript
    assocs_fdr_0.1 = assocs_fdr_0.1[which(assocs_fdr_0.1$selected_pheno %in% pheno_names$var_name),]
    
    #reorder
    assocs_fdr_0.1 = assocs_fdr_0.1[order(assocs_fdr_0.1$selected_pheno, assocs_fdr_0.1$snp_to_test, assocs_fdr_0.1$selected_model),]

    #now reorder assocs_significant to put together interactions of each phenotype. This will be the order used for the interactions. We order interactions in basis on phenotyp to have those cases with obesity or other binomial phenotype bind, because interactions with obesity have two panels each one.
    #if we have any of the binomial phenotypes within significant phenotypes
    if(TRUE %in% (binomial_pheno %in% assocs_fdr_0.1$selected_pheno)){
        #extract the row number with the binomial phneotypes and concatenate them after that the rest of row numbers
        new_row_order_fdr_0.1 = c(which(assocs_fdr_0.1$selected_pheno %in% binomial_pheno), which(!assocs_fdr_0.1$selected_pheno %in% binomial_pheno))
        print(length(new_row_order_fdr_0.1) == nrow(assocs_fdr_0.1))

        #reorder putting binomial phenotypes first
        assocs_fdr_0.1 = assocs_fdr_0.1[new_row_order_fdr_0.1,]  
    }

    #number of pages
    number_pages_single_assoc_fdr_0.1 = ceiling(nrow(assocs_fdr_0.1)/12)
    
    #figures in which we have to add the title figure
    figures_title = seq(from=11, to=nrow(assocs_fdr_0.1), by=12)
    
    #last plot where to add the title of the page
    #if the last plot is not an even number (1,3,5,...)
    if(nrow(assocs_fdr_0.1) %% 2 != 0){

        #we add the title of the last page in that plot
        last_plot_to_add_title = nrow(assocs_fdr_0.1)#The first column of plots os always 1,3,5,... so we want there the title of the figure. We check if     the last plot is even or uneve with the module. The module gives the rest of a fraction, for example, 5 mod 2 gives 1, because the fraction 5/2 has     1 as rest. In constrast, 9 mod 3 gives 0, because the fraction 9/3 has 0 as rest. Any number divided by 2 without rest is an even number, so if nrow(assocs_fdr_0.1) %% 2 is not equal to zero, the number of associations (and hence plots) is not an even number and the last plot will be in the first column of plots
    } else {

        #if not, and the last plot is even, and hence is in the second column, we select the previous plot to add the title of the last page
        last_plot_to_add_title = nrow(assocs_fdr_0.1)-1 
    }

    #add last figure to add title to figures titles if it is not included in figures_title
    if(!last_plot_to_add_title %in% figures_title){
        figures_title = c(figures_title, last_plot_to_add_title)
    }

    #plot
    pdf("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/figures/figures_fdr_0.1.pdf", height = 10, width = 7)
    par(mfrow=c(6,2), mar=c(6, 4, 1.5, 2) +0.1)
    for(i in 1:nrow(assocs_fdr_0.1)){
    
        #select the [i] row
        selected_row = assocs_fdr_0.1[i,]
    
        #select the name of the figure
        if(i %in% figures_title){#if the number of the plot is any of the figures where we include the title
    
            #extract in what page is  
            name_suple = which(figures_title == i)
    
            #establish the name
            output_name = paste("figure_S", name_suple, sep="")
    
            #if this is the last page and is not complete (we don't have 12 plots in the last page and thus it is not complete). If the rest of the division between the number of plots and 12 is not zero (nº plots mod 12), then the last page will be incomplete
            if(name_suple == length(figures_title) & ((last_plot_to_add_title+1) %% 12) != 0){

                #we put the title at a upper position (because we don't have the page complete)
                title_plot = FALSE
                title_plot_below = TRUE
            } else {#if not, i.e. it is not the last page or it is the last page but is complete

                #we put the title at the bottom
                title_plot = TRUE
                title_plot_below = FALSE
            }    
        } else {

            #in case the figure is not one of the selected, we do not add title
            title_plot = FALSE
            title_plot_below = FALSE
            output_name = NULL
        }    
         
        #run the function to plot
        plot_assoc(pheno_selected=selected_row$selected_pheno, selected_model=selected_row$selected_model, selected_snp=selected_row$snp_to_test, FDR_value = NULL, display_legend=FALSE, print_model=TRUE, title_plot=title_plot, title_plot_below=title_plot_below, output=output_name)
    }
    dev.off()
} else{ #supplementary is only for p.value figures
    number_pages_single_assoc_fdr_0.1=0
}

#####################################
### Function to plot interactions ###
#####################################
plot_interact = function(env_variable, env_var_factor, pheno_selected, selected_model, selected_snp, display_legend, legen_pos){

    #SE function to plot erro bars in bar plots
    se <- function(x) sd(x)/sqrt(length(x)) 

    #for discrete phenotypes (obesity) 
    if(pheno_selected %in% c("CVi_BP", "obesity")){
        family_error = "binomial"
        response = pheno_selected
    } else{ #for continous phenotypes
        family_error = "gaussian"
        response = paste("log(", pheno_selected, ")", sep="") #apply log transformation      
    }

    #select the row of results with the interest snp
    selected_row = interact_results_fdr_0.05[which(interact_results_fdr_0.05$pheno_selected==pheno_selected &
        interact_results_fdr_0.05$selected_model==selected_model &
        interact_results_fdr_0.05$snps==selected_snp),]

    #extract allele names
    allele_names = as.character(alleles[which(alleles$snp == selected_snp),]$helena)
    allele_names = gsub(" ", "", allele_names) #drop the space between characters
    list_allels = strsplit(allele_names, split="/")[[1]]

    #extract genotype names
    if(selected_model == "additive"){
        geno_levels = levels(eval(parse(text=paste("codominant(myData_ptpn1$",selected_snp, ")", sep=""))))

    } else { #if the heritage model is not additive

        #extract genotype names
        geno_levels = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$",selected_snp, ")", sep=""))))  
    }

    #always A is the first, because number are afabetically ordered in genotype levels in Helena
    if("A" %in% list_allels){
        geno_levels = gsub("1", "A", geno_levels)
        geno_levels = gsub("2", list_allels[!"A" == list_allels], geno_levels)                                  
    } else { #if not, C first
        if("C" %in% list_allels){ 
            geno_levels = gsub("1", "C", geno_levels)
            geno_levels = gsub("2", list_allels[!"C" == list_allels], geno_levels)                                                  
        } else { #if not G first (T always will be the last)
            if("G" %in% list_allels){
                geno_levels = gsub("1", "G", geno_levels)
                geno_levels = gsub("2", list_allels[!"G" == list_allels], geno_levels)               
            }
        }            
    }
    
    #drop the slash ("/")
    geno_levels = gsub("/", "", geno_levels)
      
    #change helena names by ncbi names (see "https://mail.google.com/mail/u/0/#search/ruizj%40ugr.es/15ee2eeabdc53653" for further details)
    #load ncbi allele names
    allele_names_ncbi = as.character(alleles[which(alleles$snp == selected_snp),]$ncbi)
    allele_names_ncbi_major = strsplit(allele_names_ncbi, split="/")[[1]][1]
    allele_names_ncbi_major = gsub(" ", "", allele_names_ncbi_major)        
    allele_names_ncbi_minor = strsplit(allele_names_ncbi, split="/")[[1]][2]
    allele_names_ncbi_minor = gsub(" ", "", allele_names_ncbi_minor)
    #load helena allele names
    allele_names_helena = as.character(alleles[which(alleles$snp == selected_snp),]$helena)
    allele_names_helena_major = strsplit(allele_names_helena, split="/")[[1]][1]
    allele_names_helena_major = gsub(" ", "", allele_names_helena_major)        
    allele_names_helena_minor = strsplit(allele_names_helena, split="/")[[1]][2]
    allele_names_helena_minor = gsub(" ", "", allele_names_helena_minor)

    #if the allele name in helena for both minor and major is different from ncbi, we change it to the complementary base (A-T, C-G). If one allele is different, the other will be also, because I only change by hand helena alleles when is clear that ncbi uses the complementary. If the major allele is A in helena and T in ncbi, and the minor is C in helena, then for sure I stablished the minor of ncbi as G (This is done BY HAND).
    if(!allele_names_ncbi_major == allele_names_helena_major){
        if(allele_names_ncbi_major == "T"){
            geno_levels = gsub("A", "T", geno_levels)
        }    
        if(allele_names_ncbi_major == "C"){
            geno_levels = gsub("G", "C", geno_levels)                
        }
        if(allele_names_ncbi_major == "G"){
            geno_levels = gsub("C", "G", geno_levels)                
        }   
        if(allele_names_ncbi_major == "A"){
            geno_levels = gsub("T", "A", geno_levels)                
        }                 
    }
    if(!allele_names_ncbi_minor == allele_names_helena_minor){
        if(allele_names_ncbi_minor == "T"){
            geno_levels = gsub("A", "T", geno_levels)
        }
        if(allele_names_ncbi_minor == "C"){
            geno_levels = gsub("G", "C", geno_levels)                
        }
        if(allele_names_ncbi_minor == "G"){
            geno_levels = gsub("C", "G", geno_levels)                
        }    
        if(allele_names_ncbi_minor == "A"){
            geno_levels = gsub("T", "A", geno_levels)                
        }             
    }

    #if the environmental variable is discrete
    if(env_var_factor==TRUE){
        #if the phenotype is not a binomial factor
        if(!family_error == "binomial"){

            ###plot a barplot of phenotype ~ snp+env_var ###
            #mean values for all genotype/env_category combination. 
            mean_values = aggregate(as.formula(paste(pheno_selected, "~", selected_model, "(", selected_row$snps, ")+", env_variable, sep="")), data = myData_ptpn1, FUN=mean) 

            #mean values for all genotype/PA_category combination. 
            se_values = aggregate(as.formula(paste(pheno_selected, "~", selected_model, "(", selected_row$snps, ")+", env_variable, sep="")), data = myData_ptpn1, FUN=se) 

            #juntamos mean y se                   
            data_mean<-cbind(mean_values, se_values[,3])  
            
            #change columns names
            colnames(data_mean)[1] <- "snp" 
            colnames(data_mean)[2] <- env_variable 
            colnames(data_mean)[3] <- "mean"
            colnames(data_mean)[4] <- "se" 
            data_mean
                               
            #Creamos una matriz a #partir de los datos.
            xtabs.var<-xtabs(mean ~factor(snp) + factor(eval(parse(text=env_variable))), data=data_mean)
            xtabs.var

            #set the plot limits depending on whether the variable has negative values or not
            if(!pheno_selected == "risk_score"){#if the variable is not the risk score and hence does not have negative values
                #max ylim of bar plot
                max_y = max(data_mean$mean)+ abs(((max(data_mean$mean)*9)/100))

                #max ylim of bar plot
                min_y = min(data_mean$mean)- abs(((min(data_mean$mean)*16)/100))
            } else {#if not
                #max ylim of bar plot
                max_y = max(data_mean$mean)+ abs(((max(data_mean$mean)*400)/100))

                #max ylim of bar plot
                min_y = min(data_mean$mean)- abs(((min(data_mean$mean)*85)/100))
            }
            
            #make the barplot
            xs<-barplot(xtabs.var, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, width=3, col=gray.colors(length(unique(data_mean$snp))), lwd=1:2, ylim=c(min_y, max_y))

            #extract levels of the env factor
            unique_env_var_labels = factor_variables_levels[which(names(factor_variables_levels) == env_variable)]

            #add x axis
            if(selected_model %in% c("codominant", "additive")){
                axis(side=1, at=c(-10, xs[2,1], xs[2,2], 40), labels=c("", unique_env_var_labels[[1]][1], unique_env_var_labels[[1]][2] , ""), cex.lab=0.9, cex.axis=0.85, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=0.8) 
            } else {
                axis(side=1, at=c(-10, (xs[1,1]+xs[2,1])/2, (xs[1,2]+xs[2,2])/2, 40), labels=c("", unique_env_var_labels[[1]][1], unique_env_var_labels[[1]][2], ""), cex.lab=0.9, cex.axis=0.85, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=0.8)                         
            }    

            #add line at zero if required        
            if(pheno_selected == "risk_score"){#if the variable is the risk score and hence has negative values add line at zero       
                abline(h=0, lty=2)
            }

            #add y axis
            if(min_y < 10){ #if y axis has low values (lower than 10) number used to calculate the tick positions (by argument in seq function) will be rounded to more decimals
                axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,2)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
            } else {
                axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,1)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
            }

            #add p.values
            mtext(text=paste("P=", round(selected_row$p.val, 6), "; FDR=",  round(selected_row$FDR, 6), sep=""), side = 3, line=0.45, outer=FALSE, cex=0.8, font=1, adj=0.48)

            #add model 
            mtext(text=paste("Under ", selected_model, " model", sep=""), side = 1, line=3.4, outer=FALSE, cex=0.9, font=1) 
            
            #y axis titles
            mtext(text=paste(pheno_names[which(pheno_names$var_name == pheno_selected),]$final_name), side = 2, line=2.8, outer=FALSE, cex=0.8)

            #add a legend
            if(display_legend==TRUE){
                if(legen_pos == "right"){
                    x_coord = xs[nrow(xs), ncol(xs)]
                    x_coord = x_coord - ((x_coord*6)/100)
                    y_coord = max_y + ((max_y*8)/100)
                    legend(x_coord, y_coord, legend=geno_levels, fill =gray.colors(nrow(xs)), bty="n", title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gen),selected_snp), cex=0.65, xpd=TRUE, xjust=-0.05)
                } else {
                    x_coord = xs[nrow(xs), ncol(xs)]
                    x_coord = x_coord - ((x_coord*95)/100)
                    y_coord = max_y + ((max_y*32)/100)
                    legend(x_coord, y_coord, legend=geno_levels, fill =gray.colors(nrow(xs)), bty="n", title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gen),selected_snp), cex=0.65, xpd=TRUE, xjust=-0.05)
                }
            }

            #extract genotype names in numbers (1/2) to calculate the number of observation
            if(selected_model == "additive"){
                geno_levels_number = sort(unique(eval(parse(text=paste("as.factor(na.omit(additive(myData_ptpn1$", selected_snp, ")))", sep="")))))
            } else { #if the heritage model is not additive

                #extract genotype names
                geno_levels_number = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$", selected_snp, ")", sep=""))))  
            }

            #extract levels of the env variable as coded in the database
            env_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", env_variable, sep="")))))

            #calculate n observations for each category of the factor
            # pheno_levels and geno_levels_number were calculated before
            #for each phenotype level. 
            n_observations = NULL
            for(p in 1:length(env_levels)){

                #select the [p] level
                selected_level_p = env_levels[p]                

                #for each genotype
                for(i in 1:length(geno_levels_number)){

                    #selected level
                    seletec_level = geno_levels_number[i]
            
                    #no_obs with first level of PA factor
                    n_observations = append(n_observations, nrow(eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")=='", seletec_level, "' & myData_ptpn1$", env_variable, "=='", selected_level_p, "' & !is.na(myData_ptpn1$", pheno_selected, ")),]", sep="")))))
                }                           
            }

            #plot the number of observations
            for(i in 1:length(as.vector(xs))){ #for each bar 

                #select the bar
                selected_xs = as.vector(xs)[i]
            
                #plot the correponding n
                mtext(paste("(", n_observations[i], ")", sep=""), side = 1, at = selected_xs, line=0.4, outer=FALSE, cex=0.7)        
            }

            #plot error bars
            for(l in 1:length(unique(eval(parse(text=paste("data_mean$", env_variable, sep="")))))){

                #for each PA level
                level_PA = unique(eval(parse(text=paste("data_mean$", env_variable, sep=""))))[l]
                subset_data = eval(parse(text=paste("data_mean[data_mean$", env_variable, "==", level_PA, ",]", sep="")))
                
                #for each genotype
                for (f in 1:nrow(subset_data)){
                    arrows(xs[f,l],subset_data$mean[f]-subset_data$se[f], xs[f,l], subset_data$mean[f]+ subset_data$se[f],angle=90,code=3,length=0.02, lwd=1)
                }                               
            } 
        } else {#else, and hence the response is a factor of two levels

            #extract genotype names
            if(selected_model == "additive"){
                geno_levels_number = sort(unique(eval(parse(text=paste("as.factor(na.omit(additive(myData_ptpn1$", selected_snp, ")))", sep="")))))
            } else { #if the heritage model is not additive

                #extract genotype names
                geno_levels_number = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$", selected_snp, ")", sep=""))))  
            }

            #extract levels of the phenotype
            pheno_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", pheno_selected, sep="")))))

            #extract levels of the env variable
            env_variable_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", env_variable, sep="")))))

            #extract label names of the pheno factor
            unique_pheno_labels = factor_variables_levels[which(names(factor_variables_levels) == pheno_selected)]

            #extract label names of the pheno factor
            unique_env_var_labels = factor_variables_levels[which(names(factor_variables_levels) == env_variable)]

            #calculate the height of the bars
            height = NULL
            env_labels = NULL            
            response_labels = NULL
            geno_labels = NULL
            for(l in 1:length(geno_levels_number)){ #for each genotype

                #select the [l] genotype
                selected_level_l = geno_levels_number[l]

                #for each level of the environmental variable
                for(e in 1:length(env_variable_levels)){

                    #select the [e] level
                    selected_level_e = env_variable_levels[e]

                    #for each phenotype level
                    for(p in 1:length(pheno_levels)){

                        #select the [p] level
                        selected_level_p = pheno_levels[p]

                        #extract number of individuals with the [l] genotype and the [p] phenotype
                        geno_subset = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")== '", selected_level_l, "' & myData_ptpn1$", env_variable, "==", selected_level_e, " & myData_ptpn1$", pheno_selected, "==", selected_level_p, "),]", sep="")))

                        #extract number of individuals with the [p] phenotype y without NA for the [l] genotype
                        full_geno = eval(parse(text=paste("myData_ptpn1[which(!is.na(", selected_model, "(myData_ptpn1$", selected_snp, ")) & myData_ptpn1$", env_variable, "==", selected_level_e, " & myData_ptpn1$", pheno_selected, "==", selected_level_p, "),]", sep="")))

                        #calculate the percentage of individuals with the [l] genotype in the group of the [p] phenotype
                        percent = round((nrow(geno_subset)/nrow(full_geno))*100, 2)
                        
                        #save that value as height of the bar
                        height = append(height, percent) 

                        #save the label of the response variable
                        response_labels = append(response_labels, unique_pheno_labels[[1]][p])

                        #save the label of the env variable
                        env_labels = append(env_labels, unique_env_var_labels[[1]][e])

                        #save the genotype
                        geno_labels = append(geno_labels, selected_level_l)               
                    }
                }    
            }

            #convert geno_labels to factor to select exactly the order of the levels and thus the order of the columns
            geno_labels = factor(geno_labels, levels=unique(geno_labels))

            #save all of them into data
            data = cbind.data.frame(response_labels, geno_labels, env_labels, height)
            data

            #for each level of the environmental variables
            for(e in 1:length(env_variable_levels)){

                #extract the number of the [e] level for the environmental variable
                selected_level_env_variable_levels = env_variable_levels[e]               
                
                #select the [e] level
                selected_level_env = unique_env_var_labels[[1]][selected_level_env_variable_levels]

                #select the bar height for the [e] level of the environmental variable
                subset_data_by_level_env = data[which(data$env_labels == selected_level_env),]

                #create a matrix from the data
                xtabs.var_env<-xtabs(height ~ geno_labels + response_labels, data=subset_data_by_level_env)
                xtabs.var_env

                #max ylim of bar plot
                max_y = max(data$height)+ ((max(data$height)*10)/100)

                #max ylim of bar plot
                min_y = min(data$height)- ((min(data$height)*13)/100)#we use the maximum and minimum value of y from the whole data with high and low physical activity to calculate the same y axis (same breaks) for both plots

                #make the barplot
                xs<-barplot(xtabs.var_env, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, width=3, col=gray.colors(length(unique(geno_levels))), lwd=1:2, ylim=c(min_y, max_y))

                #add x axis
                if(e == 1){
                    axis(side=1, at=c(2.495, ((xs[1,1]+xs[nrow(xs), ncol(xs)])/2), 40), labels=c("", selected_level_env, ""), cex.lab=0.9, cex.axis=0.85, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=1.8, xpd=TRUE) 
                } else {
                    axis(side=1, at=c(-10, ((xs[1,1]+xs[nrow(xs), ncol(xs)])/2), 18.7), labels=c("", selected_level_env, ""), cex.lab=0.9, cex.axis=0.85, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=1.8, xpd=TRUE)                          
                }    

                #add labels of the binomial phenotype 
                if(selected_model %in% c("additive", "codominant")){
                    mtext(side = 1, at=c(-10, xs[2,1], xs[2,2], 40), text = c("", unique_pheno_labels[[1]][1], unique_pheno_labels[[1]][2], ""), line = 1.1, cex=0.6)
                } else {
                    mtext(side = 1, at=c(-10, (xs[1,1]+xs[2,1])/2, (xs[1,2]+xs[2,2])/2, 40), text = c("", unique_pheno_labels[[1]][1], unique_pheno_labels[[1]][2], ""), line = 1.1, cex=0.6)                                     
                }    

                #add line at zero if required        
                if(pheno_selected == "risk_score"){#if the variable is the risk score and hence has negative values add line at zero       
                    abline(h=0, lty=2)
                }

                #add y axis, p.values and model only in the first plot
                if(e == 1){
                    #add y axis
                    if(min_y < 10){ #if y axis has low values (lower than 10) number used to calculate the tick positions (by argument in seq function) will be rounded to more decimals
                        axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,2)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
                    } else {
                        axis(side=2, at=seq(from=round(min_y,2), to=round(max_y,2), by=round((max_y-min_y)/10,1)), cex.lab=1.6, cex.axis=1, font=1, mgp=c(3, 1, 0)) #we set "at" with the sequence from y min to y max, by the difference between these values divided by 10.
                    }

                    #add the p.value and FDR
                    mtext(text=paste("P=", round(selected_row$p.val, 6), "; FDR=",  round(selected_row$FDR, 6), sep=""), side = 3, line=0, outer=FALSE, cex=0.8, font=1, adj=2)

                    #add model 
                    mtext(text=paste("Under ", selected_model, " model", sep=""), side = 1, line=3.4, outer=FALSE, cex=0.9, font=1, adj=2.4) 
                    
                    #y axis titles
                    mtext(text=paste(pheno_names[which(pheno_names$var_name == pheno_selected),]$final_name), side = 2, line=2.8, outer=FALSE, cex=0.8)
                }

                #add a legend
                if(display_legend==TRUE & e == 2){
                    if(legen_pos == "right"){
                        x_coord = xs[nrow(xs), ncol(xs)]
                        x_coord = x_coord - ((x_coord*5)/100)
                        y_coord = max_y + ((max_y*13)/100)
                        legend(x_coord, y_coord, legend=geno_levels, fill =gray.colors(nrow(xs)), bty="n", title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gen),selected_snp), cex=0.65, xpd=TRUE, xjust=-0.05)
                    } else {
                        x_coord = xs[nrow(xs), ncol(xs)]
                        x_coord = x_coord - ((x_coord*95)/100)
                        y_coord = max_y + ((max_y*32)/100)
                        legend(x_coord, y_coord, legend=geno_levels, fill =gray.colors(nrow(xs)), bty="n", title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gen),selected_snp), cex=0.65, xpd=TRUE, xjust=-0.05)
                    }
                }

                #extract genotype names in numbers (1/2) to calculate the number of observation
                if(selected_model == "additive"){
                    geno_levels_number = sort(unique(eval(parse(text=paste("as.factor(na.omit(additive(myData_ptpn1$", selected_snp, ")))", sep="")))))
                } else { #if the heritage model is not additive

                    #extract genotype names
                    geno_levels_number = levels(eval(parse(text=paste(selected_model, "(myData_ptpn1$", selected_snp, ")", sep=""))))  
                }

                #extract levels of the binomial phenotype as coded in the database
                pheno_levels = sort(unique(eval(parse(text=paste("myData_ptpn1$", pheno_selected, sep="")))))

                #calculate n observations for each category of the factor
                # pheno_levels and geno_levels_number were calculated before
                #for each phenotype level. 
                n_observations = NULL
                for(p in 1:length(pheno_levels)){

                    #select the [p] level
                    selected_level_p = pheno_levels[p]                

                    #for each genotype
                    for(i in 1:length(geno_levels_number)){

                        #selected level
                        seletec_level = geno_levels_number[i]
                
                        #no_obs with first level of pheno factor for the [e] level of the environmental variable
                        n_observations = append(n_observations, nrow(eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", env_variable, "==", selected_level_env_variable_levels, " & ", selected_model, "(myData_ptpn1$", selected_snp, ")=='", seletec_level, "' & myData_ptpn1$", pheno_selected, "==", selected_level_p, "),]", sep="")))))
                    }                           
                }

                #plot the number of observations
                for(i in 1:length(as.vector(xs))){ #for each bar 

                    #select the bar
                    selected_xs = as.vector(xs)[i]
                
                    #plot the correponding n
                    mtext(paste("(", n_observations[i], ")", sep=""), side = 1, at = selected_xs, line=0.1, outer=FALSE, cex=0.7)        
                }
            }               
        }     
    } else { #FROM HERE CODE NOT REVISED
        if(!pheno_selected=="CVi_BP"){

            #set name of phenotype
            ylab = pheno_names[which(pheno_names$var_name==pheno_selected),]$final_name
            xlab=ssb_names[which(ssb_names$var_name==env_variable),]$final_name

            #plot response versus env_variable
            plot(eval(parse(text=response))~eval(parse(text=env_variable)), data=myData_ptpn1, ylab=ylab, xlab="", cex.lab=1.15)
            
            #extract SNP levels
            snp_levels = unique(eval(parse(text=paste("na.omit(", selected_model, "(myData_ptpn1$", selected_row$snps, "))", sep=""))))

            #for each SNP level
            for(l in 1:length(snp_levels)){

                #select the [i] level
                selected_level = snp_levels[l]

                #subset in basis on that level
                subset_snp = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_row$snps, ") =='", selected_level, "'),]", sep="")))

                #plot points of that subset with the [l] color of the current palette
                points(eval(parse(text=response))~eval(parse(text=env_variable)), data=subset_snp, col=palette()[l])

                #plot regression line
                if(nrow(subset_snp)>=6){#print lines only if there are a minimum of points (6 per predictor) 
                    abline(glm(eval(parse(text=response))~eval(parse(text=env_variable)), data=subset_snp, family=family_error), col=palette()[l])
                }    
            }

            #add the legend
            if(display_legend==TRUE){
                #legend to get the position of topright in this plot
                first_legend = legend("topright", legend=geno_levels, fill=palette()[1:length(snp_levels)], title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), selected_row$snps, sep=" "), plot=FALSE)
                #set exact position of the legend
                min_x = min(eval(parse(text=paste("na.omit(myData_ptpn1$", env_variable, ")", sep=""))))
                max_x = max(eval(parse(text=paste("na.omit(myData_ptpn1$", env_variable, ")", sep=""))))
                min_y = min(eval(parse(text=paste(transformation, "(na.omit(myData_ptpn1$", pheno_selected, "))", sep=""))))
                max_y = max(eval(parse(text=paste(transformation, "(na.omit(myData_ptpn1$", pheno_selected, "))", sep=""))))            
                x_pos_legend = first_legend$rect$left-(abs(max_x-min_x)*(-6)/100)
                y_pos_legend = first_legend$rect$top-(abs(max_y-min_y)*1/100)
                #plot the legend
                legend(x_pos_legend, y_pos_legend, legend=geno_levels, fill=palette()[1:length(snp_levels)], title=paste(toupper(gene_names[which(gene_names$selected_snp == selected_snp),]$gene), selected_row$snps, sep=" "), plot=TRUE, cex=0.8, bty="n")
            }

            #add results as text
            mtext(paste("P=", round(selected_row$p.val,4), ", ", "FDR=", round(selected_row$FDR,4), sep=""), side=3, line=-1.5, cex=0.6)

            #add xlab
            mtext(text=xlab, side = 1, line=2, outer=FALSE, cex=0.75, font=1, adj=0.449) 

            #add model
            mtext(text=paste("Under ", selected_model, " model", sep=""), side = 1, line=3.2, outer=FALSE, cex=0.75, font=1, adj=0.449)             
        } else {
            
            ###plot a barplot of phenotype ~ snp+PA ###
            #mean values for all genotype/PA_category combination. 
            mean_values = aggregate(as.formula(paste(env_variable, "~", response, "+", selected_model, "(", selected_row$snps, ")", sep="")), data = myData_ptpn1, FUN=mean) 

            #mean values for all genotype/PA_category combination. 
            se_values = aggregate(as.formula(paste(env_variable, "~", response, "+", selected_model, "(", selected_row$snps, ")", sep="")), data = myData_ptpn1, FUN=se)

            #juntamos mean y se                   
            data_mean<-cbind(mean_values, se_values[,3])  
    
            #change columns names
            colnames(data_mean)[1] <- response 
            colnames(data_mean)[2] <- "snp" 
            colnames(data_mean)[3] <- "mean"
            colnames(data_mean)[4] <- "se" 
            data_mean
                               
            #Creamos una matriz a #partir de los datos.
            xtabs.var<-xtabs(mean ~factor(snp) + factor(eval(parse(text=response))), data=data_mean)
            xtabs.var

            #max ylim of bar plot
            max_y = max(data_mean$mean)+ ((max(data_mean$mean)*40)/100)

            #make the barplot
            xs<-barplot(xtabs.var, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, width=3, col=gray.colors(length(unique(data_mean$snp))), lwd=1:2, ylim=c(0, max_y))

            #add axis
            if(selected_model %in% c("additive", "codominant")){
                axis(side=1, at=c(-10, xs[2,1], xs[2,2], 40), labels=c("", "Non-overweight", "Overweight", ""), cex.lab=0.9, cex.axis=1, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=0.8) 
            } else {
                axis(side=1, at=c(-10, (xs[1,1]+xs[2,1])/2, (xs[1,2]+xs[2,2])/2, 40), labels=c("", "Non-overweight", "Overweight", ""), cex.lab=0.9, cex.axis=1, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0), padj=0.8)                         
            }    
            axis(side=2, cex.lab=1.6, cex.axis=1, font=2, mgp=c(3, 1, 0))

            #add main title
            mtext(text=paste("P=", round(selected_row$p.val,4), ", ", "FDR=", round(selected_row$FDR,4), sep=""), side = 3, line=1, outer=FALSE, cex=0.8, font=4)
    
            #y axis titles
            #set name of ssb variable
            xlab=ssb_names[which(ssb_names$var_name==env_variable),]$final_name            
            mtext(text=xlab, side = 2, line=2.8, outer=FALSE, cex=0.8)

            #add model
            mtext(text=paste("Under ", selected_model, " model", sep=""), side = 1, line=3.4, outer=FALSE, cex=0.9, font=1, adj=0.449) 

            #add the legend
            #legend to get the position of topright in this plot
            first_legend = legend("topleft", legend=geno_levels, fill =gray.colors(nrow(xs)), plot=FALSE)
            #set exact position of the legend
            min_y = min(eval(parse(text=paste("na.omit(myData_ptpn1$", env_variable, ")", sep=""))))
            max_y = max(eval(parse(text=paste("na.omit(myData_ptpn1$", env_variable, ")", sep=""))))
           
            x_pos_legend = first_legend$rect$left
            y_pos_legend = first_legend$rect$top-(abs(max_y-min_y)*1/100)
            #plot the legend
            legend(x_pos_legend, y_pos_legend, legend=geno_levels, fill =gray.colors(nrow(xs)), plot=TRUE)

            #create subsets of physical activity
            no_PA_subset = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", response, "==0 & !is.na(myData_ptpn1$", selected_row$snps, ")),]", sep="")))
            PA_subset = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$", response, "==1 & !is.na(myData_ptpn1$", selected_row$snps, ")),]", sep="")))                    

            #calculate n observations for each category of the factor
            if(!selected_model == "additive"){ #if the model is not additive we use genotype object of SNPassoc

                #extract levels of the snp
                levels_snp = levels(eval(parse(text=paste("na.omit(", selected_model, "(myData_ptpn1$", selected_row$snp, "))", sep=""))))
                
                #loop for extracting n_obs of each genotype
                n_obs_no_PA = NULL
                n_obs_PA = NULL
                for(i in 1:length(levels_snp)){

                    #selected level
                    seletec_level = levels_snp[i]
                
                    #no_obs with first level of PA factor
                    n_obs_no_PA = append(n_obs_no_PA, length(eval(parse(text=paste("na.omit(no_PA_subset[", selected_model, "(no_PA_subset$", selected_row$snp, ") == '", seletec_level, "',]$", pheno_selected,")", sep="")))))
                    
                    #no_obs with second level of PA factor
                    n_obs_PA = append(n_obs_PA, length(eval(parse(text=paste("na.omit(PA_subset[", selected_model, "(PA_subset$", selected_row$snp, ") == '", seletec_level, "',]$", pheno_selected,")", sep="")))))                        
                }
            } else { #if model is additive we have to transform the variable to a factor
                
                #extract levels of the snp                    
                levels_snp = levels(as.factor(eval(parse(text=paste("na.omit(", selected_model, "(myData_ptpn1$", selected_row$snp, "))", sep="")))))

                #loop for extracting n_obs of each genotype                        
                n_obs_no_PA = NULL
                n_obs_PA = NULL
                for(i in 1:length(levels_snp)){

                    #selected level                            
                    seletec_level = levels_snp[i]
                
                    #no_obs with first level of PA factor                        
                    n_obs_no_PA = append(n_obs_no_PA, length(eval(parse(text=paste("na.omit(no_PA_subset[", selected_model, "(no_PA_subset$", selected_row$snp, ") == '", seletec_level, "',]$", pheno_selected,")", sep="")))))

                    #no_obs with second level of PA factor                            
                    n_obs_PA = append(n_obs_PA, length(eval(parse(text=paste("na.omit(PA_subset[", selected_model, "(PA_subset$", selected_row$snp, ") == '", seletec_level, "',]$", pheno_selected,")", sep="")))))                        
                }                 
            }

            #bind all observations. 
            n_observations = c(n_obs_no_PA, n_obs_PA)

            #plot the number of observations
            for(i in 1:length(as.vector(xs))){ #for each bar 

                #select the bar
                selected_xs = as.vector(xs)[i]
        
                #plot the correponding n
                mtext(paste("(", n_observations[i], ")", sep=""), side = 1, at = selected_xs, line=0.4, outer=FALSE, cex=0.7)        
            }

            #plot error bars
            for(l in 1:length(unique(eval(parse(text=paste("data_mean$", response, sep="")))))){

                #for each PA level
                level_PA = unique(eval(parse(text=paste("data_mean$", response, sep=""))))[l]
                subset_data = eval(parse(text=paste("data_mean[data_mean$", response, "==", level_PA, ",]")))
                
                #for each genotype
                for (f in 1:nrow(subset_data)){
                    arrows(xs[f,l],subset_data$mean[f]-subset_data$se[f], xs[f,l], subset_data$mean[f]+ subset_data$se[f],angle=90,code=3,length=0.02, lwd=1)
                }                               
            }                  
        }      
    }
} #IMPORTANTE: CODE FOR CONTINUOUS ENVIRONMENTAL VARIABLE IS NOT REVISED


###################
#### Figures 5 ####
###################

####plot associations significant according to FDR<0.1
final_interactions

#select only those phenotypes that I decided to include in the manuscript
final_interactions = final_interactions[which(final_interactions$pheno_selected %in% pheno_names$var_name),]
#check
summary(unique(final_interactions$pheno_selected) %in% pheno_names[-which(pheno_names$var_name %in% c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap")),]$var_name)#all TRUE

#reorder now rows to get bind those rows with the same phenotype
final_interactions = final_interactions[order(final_interactions$pheno_selected, final_interactions$snps, final_interactions$selected_model),]

#if we have binomial phenotypes: reorder final_interactions to put together interactions of each phenotype. We order in basis on phenotyp to have those cases with obesity or other binomial phenotype bind.
#if we have any of the binomial phenotypes within significant phenotypes
if(TRUE %in% (binomial_pheno %in% final_interactions$pheno_selected)){

    #extract the row number with the binomial phneotypes and concatenate them after that the rest of row numbers
    new_row_order = c(which(final_interactions$pheno_selected %in% binomial_pheno), which(!final_interactions$pheno_selected %in% binomial_pheno))
    length(new_row_order) == nrow(final_interactions)

    #reorder putting binomial phenotypes first
    final_interactions = final_interactions[new_row_order,]
}

#number of panels
n_panels = (length(which(final_interactions$pheno_selected %in% binomial_pheno))*2) + length(which(!final_interactions$pheno_selected %in% binomial_pheno))#calculate the number of panels, two panels for each significant association of binomial phenotypes) and then one panel for each of the other associations

#number of pages
number_pages_interact_PA = ceiling(n_panels/12)

#plot
cairo_pdf("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/figures/figures_interact_PA_factor.pdf", height = 10, width = 7, onefile=TRUE)#onefile=TRUE to get all pages in the same file
par(mfrow=c(6,2), mar=c(5, 4, 2, 1.5) +0.1)
for(i in 1:nrow(final_interactions)){

    #select the [i] row
    selected_row = final_interactions[i,]

    #plot
    plot_interact(env_variable="PA_factor", env_var_factor=TRUE, pheno_selected=selected_row$pheno_selected, selected_model=selected_row$selected_model, selected_snp=selected_row$snps, display_legend=TRUE, legen_pos="right")
}
dev.off()



############################################################
########### FIGURES P.VALUES GEN-PHENOTYPE #################
############################################################
plot_pvals = function(pheno, output){

    #extract pvals for the phenotype
    pvals = geno_pheno_results[which(geno_pheno_results$selected_pheno == pheno),]

    #open a pdf
    pdf(paste("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals/", output, ".pdf", sep=""), width=7, height=9)

    #set graphic parameters
    par(mfcol=c(5,1), mar=c(4,2,0.5,2)+0.1, oma = c(1.5,0.5,1.1,0))

    #models to plot
    heritage_models = c("codominant", "dominant", "recessive", "overdominant", "additive")

    #loop for making a figure for each model (each column)
    for(i in 1:length(heritage_models)){

        #select the [i] model
        selected_model = heritage_models[i]

        #extract the pvals for the [i] model
        selected_pvals = pvals[which(pvals$selected_model==selected_model),]

        #change mar in the bottom of the figure (last plot)
        if(i == length(heritage_models)){
            par(mar=c(9,2,0.5,2)+0.1)
        }

        #selected model
        model_name = selected_model

        #create a dataset with -log10 p vals and false x index to make a scatter plot instead of a box plot
        x_index = 0:(nrow(selected_pvals)-1) #0 is the first SNPs, 17 is the last (we have 18)
        log_pvals = -log(selected_pvals$pvals, 10) #calcualte -log in 10 degree of p.values
        data = cbind.data.frame(x_index, log_pvals) #bind  both 

        #extract names of SNPs
        xlabels = selected_pvals$snp_to_test

        #calculate nominal and bofferoni threshold
        nominal_limit = -log(0.05,10)
        bonferroni_limit = -log((0.05/length(na.omit(data$log_pvals))), 10) #only take into account snps without NA, thus in non-dominant models ill be 15 only
        
        #calcualte highest p.value
        max_pval = max(na.omit(data$log_pvals))

        #make plot
        plot(log_pvals~x_index, data, type="b", axes=FALSE, main=model_name, ylab="", xlab="", ylim=c(0, round(max_pval+1.1,0)))

        #add x axis
        if(!i == 5){ #only add SNPs names for the last plot 
            axis(1, at=seq(0,(length(xlabels)-1),by=1), labels=rep("", length(xlabels)), las = 2)
        } else {
            axis(1, at=seq(0,(length(xlabels)-1),by=1), labels=xlabels, las = 2, cex.axis=1)
        }

        #add Y axis
        axis(2, at=seq(0, round(max_pval+1.1,0), 1), line=-1.7)

        #add Y title if it is the plot 3
        if(i == 3){
            mtext(expression(-log[10]~"(p value)"), side=2, padj=-0.5)
        }

        #add lines of threshold only if they are close to the highest p.value
        if((nominal_limit - max_pval) < 0.7){
            axis(side=1, at=c(0, 34.005), labels=FALSE, cex.lab=1, cex.axis=1, font=1, lty=2, pos=nominal_limit, lwd.ticks=0, outer=FALSE)
        }
        if((bonferroni_limit - max_pval) < 0.7){           
            axis(side=1, at=c(0, 34.005), labels=FALSE, cex.lab=1, cex.axis=1, font=1, lty=3, pos=bonferroni_limit, lwd.ticks=0, outer=FALSE)
        }  

        if(i == 2){             
            legend(x=length(labels(myData_ptpn1))-2, y=1.2*round(max_pval+1.1,0), legend=c("Nominal p value", "Bonferroni correction"), lty=c(2,3), bty="n", xpd=TRUE)#X position is the position of the antepenultimate SNP, i.e., length of SNP names less 2
        }
    }

    #add lines of genes
    #axis(1,at=c(-0.2,7.2),col="blue",line=5.7,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
    #axis(1,at=c(7.5, 11.3),col="blue",line=5.7,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
    #axis(1,at=c(11.7, 17.3),col="blue",line=5.7,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)

    #add names of genes
    #mtext("ptpn11", line= -11.5, adj=0.2, cex=1.5, font=4)
    #mtext("ptpn12", line= -11.5, adj=0.55, cex=1.5, font=4)
    #mtext("ptpn13", line= -11.5, adj=0.87, cex=1.5, font=4)

    #set number figure
    number_figure = strsplit(output, split="_")[[1]][2]

    #insert figure name
    mtext(paste("Online supplementary figure ", number_figure, sep=""), side=1, font=2, cex=1.5, adj=0.015, padj=1.5, outer=TRUE, line=-3)

    dev.off()
}

#list of phenotypes
pheno_figures = c("obesity","CRF_BMI","CRF_waist","waist_height","CRF_hip","waist_hip","CRF_Body_fat_PC","FMI","TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","HOMA","QUICKI","Leptin_ng_ml","SBP","DBP")
length(pheno_figures) == nrow(pheno_names[-which(pheno_names$var_name %in% c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap")),])

#set if you have additional suplementary figures. For example, the haplotype was included in supplemenetary in PTPN1 paper because we had too many associations and corresponding figures in main text
additional_supple_figures = 1

#sequence of figure numbers: number list from the number of previous supplementaries "MORE 1" to the total number of phenotypes "MORE the number of previous supplementaries"
sequence_pheno = (number_pages_single_assoc_fdr_0.1+additional_supple_figures+1):(length(pheno_figures)+number_pages_single_assoc_fdr_0.1+additional_supple_figures)

#for each phenotype
for(k in 1:length(pheno_figures)){

    #select the [k] phenotype
    selected_phenotype = pheno_figures[k]
    
    #paste with figure_S
    output_name = paste("figure_S", sequence_pheno[k], sep="")   

    #print the number of the phenotype and the figure number
    print(paste(selected_phenotype, "_", output_name, sep=""))

    #plot
    plot_pvals(pheno=selected_phenotype, output=output_name)
}


#######################################
######## BIND ALL PDS OF SUPLE ########
#######################################

#extract path of all species with phylo ensamble. In with and without debería haber el mismo número de especies, por eso cogemos una de las carpetas 
path_suples = list.files("/Users/diegosalazar/My Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals", pattern="figure_S")
#sort supple names by the number (we don't want 10 as the first)
library(gtools)
path_suples = mixedsort(sort(path_suples))

#check that all figures has been created
n_vars = nrow(pheno_names[-which(pheno_names$var_name %in% c("center", "CRF_age", "CRF_weight", "CRF_height", "CRF_trici", "CRF_subscap")),])
length(path_suples) == (n_vars)

#extract the number of the first and last page of the pval supplementaries
if(length(path_suples) == (n_vars)){#if we have all the phenotypes in the folder

    #extract the number of supple figures of pvals using the number of pages with supple barplots and the number of phenotypes studied "n_vars"
    first_pval_suple_number = number_pages_single_assoc_fdr_0.1+additional_supple_figures+1
    second_pval_suple_number = number_pages_single_assoc_fdr_0.1+additional_supple_figures+n_vars
} else { #if not
    
    #print an error and stop
    stop("ALL PHENOTYPES ANALYZED ARE NOT INCLUDED IN THE FOLDER")
}
#bind all pdf generated into one single pdf file
system(paste("rm /Users/diegosalazar/My\\ Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals/supplementary_material_full_geno_pheno_pvals.pdf;

    cd /Users/diegosalazar/My\\ Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals; 

    pdftk $(for n in {", first_pval_suple_number, "..", second_pval_suple_number, "}; do echo figure_S$n.pdf; done)  cat output /Users/diegosalazar/My\\ Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals/supplementary_material_full_geno_pheno_pvals.pdf", sep="")) #first delete the prevous full file in full directory; next change the drectory to the that with separate plots; finally run pdftk to bind all pdf in that directory and save the result in full directory
    #lo de echo es para que se ornden del 1 al 19 las figuras.

#calculate the number of pages in the pdf created
library(pdftools)
info_pdf <- pdf_info("/Users/diegosalazar/My\ Drive/science/open_projects/helena_7/results/figures/figures_suple_pvals/supplementary_material_full_geno_pheno_pvals.pdf")
info_pdf$pages

#check if this number match the number of studied phenotypes
if(info_pdf$pages == n_vars){

    #print message that everything is ok
    print("All phenotypes are included in the pdf")
} else {

    #print err message
    stop("All phenotypes are NOT included in the pdf")
}