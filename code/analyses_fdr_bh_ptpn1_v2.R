#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



#######################################################################################################################
######### SCRIPT FOR MODELING GENE-PHENOTYPE AND GENE*ENV INTERACTIONS FOR PTPN1 AND OBESITY-RELATED TRAITS ###########
#######################################################################################################################

#In this script, we model the association between PTPN1 polymorphisms and obesity related traits in the HELENA cohort. We also consider the association between haplotypes of these SNPs and traits, and the interaction between these SNPs and physical activity. 




##############################################################
######### DIFFERENCES RESPECT TO PREVIOUS VERSIONS ###########
##############################################################

#Version 1
    #We add here the calculation of R2 for each gene-phenotype association and each gene*physical activity interactions to create supplementary data 1 and 2, respectively. 




#######################################
######### REQUIRED PACKAGES ###########
#######################################

require(foreign) #for importing the HELENA database from SPSS format
require(SNPassoc) #for using the heritage model functions and create setupSNP objects
    #Not all the packages are listed!



#######################################
######### DATA PREPARATION ############
#######################################

######clean the working space#########
#We remove all the elements except wideScreen, which is a command to expand the space of the terminal in R
rm(list = ls()[-which(ls()=="wideScreen")])

########### Data for analysing ucp genotype data from helena
setwd("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7")

#### import data from SPSS database
require(foreign)
helena <- read.spss("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_sweet/data/snps/Helena_CSS_no_UK_Modena_13_12_20-2.sav", to.data.frame=TRUE)
str(helena)
summary(helena)

#database used in the first two helenas
if(FALSE){
    helena = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_ucp_obesity/data/snps/helena_db.txt", sep="\t", header=T)
    str(helena)
    summary(helena)
}

### faster way to extract the data
if(FALSE){
    n<-3551 #number of rows without the header (e.g. number of individuals) 
    dat<-scan("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/helena_db.txt", list("character"),skip=1) 
    variables<-scan("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/helena_db.txt", list("character"), n=1) 
    ncols<-length(dat[[1]])/n
    temp<-matrix(dat[[1]],nrow=n,ncol=ncols,byrow=TRUE) 
    HapMap<-data.frame(temp, stringsAsFactors = FALSE) 
    dimnames(HapMap)[[2]]<-variables[[1]]
}

### load SNPassoc
require(SNPassoc)

##############################################################
######## DATA MANIPULATION AND DESCRIPTIVE ANALYSIS ##########
##############################################################

###########################################
######## Create PA covariate ##############
###########################################
#This variable is used in GEN X ENVIRONMENT INTERACTION
PA_factor = NULL
for(i in 1:length(helena$MVPA_mean)){ 
    
    #for each value of MVPA_mean
    MVPA_selected = helena$MVPA_mean[i] 

    #if [i] value of MVPA_mean is NA
    if(is.na(MVPA_selected)){

        #PA_factor is NA
        PA_factor = append(PA_factor, NA)

    #if [i] value of MVPA_mean is not NA        
    } else {

        #if [i] value of MVPA_mean lower than 60
        if(MVPA_selected < 60){
            
            #PA_factor is 0 
            PA_factor = append(PA_factor, 0)        

        #if [i] value of MVPA_mean is not lower than 60 (i.e. equal or higher than 60)
        } else {

            #PA_factor is 1
            PA_factor = append(PA_factor, 1)

        }
    }
}
helena$PA_factor = factor(PA_factor)
rm(PA_factor)

## test
unique(helena[helena$MVPA_mean >= 60, ]$PA_factor) #only active individuals (1) or NA
unique(helena[helena$MVPA_mean < 60, ]$PA_factor) #only non active individual (0) or NA

##############################################################
######## DATA MANIPULATION AND DESCRIPTIVE ANALYSIS ##########
##############################################################

###########################################
######## Create some phenotypes ###########
###########################################

##### recode BMI cat variable with the purpose of converting in a obesity variable
## recode BMI_cat variable to two levels: i) obesity-overweight with BMI equal or higher than 25; ii) normal weight with BMI lower than 25
new_obesity = NULL #empty vector
for(i in 1:length(helena$BMI_cat)){

    #select the [i] BMI value
    BMI = helena$BMI_cat[i]

    #if BMI indicate obese or overweight
    if(BMI == "Obese [>30 in adults]" | BMI == "Overweight [25-30 in adults]"){
        new_obesity = append(new_obesity, 1)
    } else { #if not
        new_obesity = append(new_obesity, 0)
    }
}
helena$obesity = factor(new_obesity)#save it as a factor
rm(new_obesity) #remove the object

## test 
unique(helena[which(helena$obesity == 1), ]$BMI_cat) #only obese and overweight
unique(helena[which(helena$obesity == 0), ]$BMI_cat) #only optimal BMI and BMI too low (cases with too low BMI were very few)

##### create FMI (fat mass index)
#fat mass in kg divided by (height in meter)^2
helena$FMI = helena$fat_mass_kg/((helena$CRF_height/100)^2)
#check: obtain height from FMI reverting the calculus (sqrt(fat_mass/FMI)*100)
summary(round(helena$CRF_height - (sqrt(helena$fat_mass_kg/helena$FMI))*100,4) == 0)

####create waist/height ratio
helena$waist_height = helena$CRF_waist/helena$CRF_height
#check reversing the formula
summary(round((helena$waist_height*helena$CRF_height) - helena$CRF_waist,2) == 0)

####create waist/hip ratio
helena$waist_hip = helena$CRF_waist/helena$CRF_hip
#check reversing the formula
summary(round((helena$waist_hip*helena$CRF_hip) - helena$CRF_waist, 2) == 0)

####info about calculation of QUICKI here: "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_6/code/check_quicki"

#### create BP variables (mean of the two values of each variable)
#We two measurements of each blood pressure variable. in both cases, first and second measurements of blood pressures are similar
plot(helena$CRF_BPsys1,helena$CRF_BPsys2)
plot(helena$CRF_BPdias1,helena$CRF_BPdias2)
dev.off()

#mean of the two values of each variable
helena$SBP = (helena$CRF_BPsys1+helena$CRF_BPsys2)/2
helena$DBP = (helena$CRF_BPdias1+helena$CRF_BPdias2)/2
#check
summary(helena$SBP == rowMeans(helena[c('CRF_BPsys1', 'CRF_BPsys2')], na.rm=TRUE))
summary(helena$DBP == rowMeans(helena[c('CRF_BPdias1', 'CRF_BPdias2')], na.rm=TRUE))

#### TC_HDL variable
helena$TC_HDL = helena$TC/helena$HDL
#check reversing the formula
summary(round((helena$TC_HDL*helena$HDL) - helena$TC, 2) == 0)

#### sum of two skinfolds: tricipital skin fold (in triceps) and  sub-scapular skin fold (behind the scapula)
helena$sum_two_skinfolds = helena$CRF_trici + helena$CRF_subscap #I have selected these skinfolds because Jonatan used them for calculating body fat percentage in "http://www.spanishexernet.com/pdf/Adolescentes%20nuevo/RuizJR_2010APAM_FTO%20Fatness%20&%20PA-HELENA.pdf"
#check reversing the formula
summary(round((helena$sum_two_skinfolds-helena$CRF_subscap) - helena$CRF_trici, 2) == 0)

#VLDL
summary((helena$TG/5) == helena$VLDL)#In helena, VLDL has been estimated from TG, as TG/5. Therefore, we don't have a direct measure of VLDL, and we are repeating analysis with exactly the same variable. Decision: We remove VLDL from the analysis.

### create heart risk score
#variable used to calculate the score
score_variables = c(
    "TC_HDL",
    "TG",
    "HOMA",
    "SBP", 
    "sum_two_skinfolds")

#create a data.frame with standarized variables (value-mean)/sd
standarized_variables = as.matrix(rep(NA,nrow(helena))) #as many rows as rows has helena

#for each score variable
for(i in 1:length(score_variables)){

    #selected var name
    selected_var_name = score_variables[i]

    #select variable extracting data with parse y eval
    selected_var = eval(parse(text=paste("helena$", selected_var_name, sep="")))
    
    #calculate mean of the [i] variable
    mean_var = mean(na.omit(selected_var))

    #calculate sd of the [i] variable    
    sd_var = sd(na.omit(selected_var))

    #standarize each value of the [i] variable
    standarized_value = NULL
    for(k in 1:length(selected_var)){

        #select [k] value
        selected_value = selected_var[k]

        #substract the mean of the [i] variable from the [k] value y divide by sd of the [i] variable
        standarized_value = append(standarized_value, ((selected_value - mean_var)/sd_var))#you calculate the difference between the interest value and the mean. If the value si similar to the mean, the result would close to zero (small), while the value is very different, the difference would be bigger. This difference would be then corrected by the standard deviation, a small sd would reduce the total (as sd is the denominator is smaller). This makes becuase a small SD indicate little variation around the mean. In contrast, high sd would increase the result, incresing the difference, which makes sense as there is more variation around the mean.
    }

    #check that the calculation of the variable is correct reverting the formula for calculating the standarized_value and the obtain the selected_var
    print(summary(((standarized_value*sd_var)+mean_var) - selected_var))

    #save it
    standarized_variables = cbind.data.frame(standarized_variables, standarized_value)
}
standarized_variables = standarized_variables[,-1]
colnames(standarized_variables) <- paste(score_variables, "_standarized", sep="")
str(standarized_variables)

#calculate the score for each sample (individual)
risk_score = NULL
for(i in 1:nrow(standarized_variables)){
    #select the sample and extract standarized values of four variables as a vector
    selected_row = as.numeric(as.vector(standarized_variables[i,]))
    
    #if there is not any NA
    if(!TRUE %in% is.na(selected_row)){

        #calculate the mean of the standarized variables
        risk_score = append(risk_score, mean(selected_row))
    } else { #if not

        #add a NA
        risk_score = append(risk_score, NA)
    }
}

#negative values of risk score, Jonatan response: El risk score puede tener datos negativos! of course! un valor negativo en el score significa que la persona está por debajo de la media. Es un valor REAL y no se puede modificar. Piensa en cómo se calcula: valor individual menos el valor medio de la población dividido todo por la SD. Por lo tanto un valor negativo quiere decir que el individuo está por debajo de la media de su población.

#Therefore, we calculate a risk score without negative values for applying the log transformation (you cannot calculate the log of negative values). This transformed variable will be used for analyses, but not for graphical representation
risk_score_for_log = risk_score + abs(floor(min(na.omit(risk_score))))
#check that risk_score_for_log less the minum value of ris_score is equal to risk score
summary((risk_score - (risk_score_for_log-abs(floor(min(na.omit(risk_score)))))) < 0.000005)

#see correlation between score variables
plot(risk_score~risk_score_for_log)
dev.off()
length(risk_score) == length(risk_score_for_log)

#add the score to helena
helena$risk_score <- risk_score
helena$risk_score_for_log <- risk_score_for_log

#see correlation of risk score with score variables
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/risk_score/risk_score_cors.pdf")
par(mfrow=c(2,2))
plot(risk_score~TC_HDL, helena)
plot(risk_score~TG, helena)
plot(risk_score~HOMA, helena)
plot(risk_score~SBP, helena)
plot(risk_score~sum_two_skinfolds, helena)
dev.off()

#remove vector of the operations
rm(risk_score)
rm(risk_score_for_log)
rm(standarized_variables)

#### convert into factors those covariables that are factors. Response that are 0-1 will not convert because they have to be numeric for SNPassoc. 
helena$CRF_sex <- as.factor(helena$CRF_sex)
helena$center <- as.factor(helena$center)

##################################
######## Subset helena ###########
##################################

## extract the columns of SNPs
col_first_snp = which(colnames(helena) == "Ch_1_6_rs1537516")
col_second_snp = which(colnames(helena) == "Ch_15_X_seqSNP9")
col_snps = col_first_snp:col_second_snp

## create vector with SNPs varialbes without chromosome position
subset_helena = helena[, col_snps] #snps de helena
require(stringr)
ancient_labels = NULL
new_labels= NULL
for (i in 1:length(colnames(subset_helena))){ #loop for create the new label variable
    subset=colnames(subset_helena)[i]
    ancient_labels = append(ancient_labels, subset)
    new_labels = append(new_labels, str_split_fixed(subset, "_", 4)[4])
}
snps_labels = data.frame(cbind(ancient_labels, new_labels))

#loop for checking that the changing of the snps occurred well
check_change_names = NULL
for (i in 1:nrow(snps_labels)){

    #selected row
    selected_row = snps_labels[i,]

    #extract the chromosome of the ancient name of the [i] row
    chr_selected_ancient_name = paste(strsplit(as.character(selected_row$ancient_labels), "_", 3)[[1]][1], "_", strsplit(as.character(selected_row$ancient_labels), "_", 3)[[1]][2], "_", strsplit(as.character(selected_row$ancient_labels), "_", 3)[[1]][3], sep="")

    #the combination of the chromosome name with the rs number produce the ancient label? Save the result of the test in the vector
    check_change_names = append(check_change_names, paste(chr_selected_ancient_name, "_", selected_row$new_labels, sep="") == selected_row$ancient_labels)
}
summary(check_change_names) #All TRUE, perfect.

#change names of SNPs, now only rs. 
require(data.table)
setnames(helena, old = as.vector(snps_labels$ancient_labels), new = as.vector(snps_labels$new_labels))

#indicate snps of ptpn1 genes
ptpn1_snps = data.frame(cbind(c(rep("ptpn1", 7)), c(c("rs6067472", "rs10485614", "rs2143511", "rs6020608", "rs968701", "rs6512654", "rs3787345"))))
colnames(ptpn1_snps) <- c("gen", "snp")
nrow(ptpn1_snps)

#check that are all correct with snp-gene data from the whoel panel (extracted from ncbi in "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_sweet/code/analyses_fdr_bh.R")
chromosome_snps_helena = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/chromosome_snps.csv", sep=",", header=T)

#snps not included in chromosome_snps_helena or in ptpn1_snps (we have to check both possibilities)
different_snps = c(
    setdiff(as.vector(chromosome_snps_helena[which(chromosome_snps_helena$gene==toupper("ptpn1")),]$selected_snp), as.vector(ptpn1_snps$snp)),
    setdiff(as.vector(ptpn1_snps$snp), as.vector(chromosome_snps_helena[which(chromosome_snps_helena$gene==toupper("ptpn1")),]$selected_snp)))#gene name in upper case becuse gene names are in that way in "chromosome_snps_helena"

#from these snps, what are included in the helena
different_snps[which(different_snps %in% snps_labels$new_labels)] #these snps are errors from myself, you have to included in ptpn1_snps. ADD the snps not incldued by hand.

#The rest are not genotyped in helena so we are not interested in them
snps_not_genotyped = different_snps[which(!different_snps %in% snps_labels$new_labels)] #these snps are not genptyped

#drop these snps
ptpn1_snps = ptpn1_snps[which(!ptpn1_snps$snp %in% snps_not_genotyped),]
rownames(ptpn1_snps) <- seq(length=nrow(ptpn1_snps)) #reset numeration of rows
ptpn1_snps

#check that all snps of ptpn1_snps are included in helena
summary(ptpn1_snps$snp %in% snps_labels$new_labels)

#check that all snps of ptpn1 of helena are included in ptpn1_snps
summary(chromosome_snps_helena[which(chromosome_snps_helena$gene==toupper("ptpn1")),]$selected_snp %in% ptpn1_snps$snp)#gene name in upper case becuse gene names are in that way in "chromosome_snps_helena"

#what columns of helena have the selected SNPs
helena_7_cols = which(colnames(helena) %in% ptpn1_snps$snp) 
length(helena_7_cols) == nrow(ptpn1_snps)

#create a data frame with position of genes
ptpn1_position = data.frame(cbind(as.vector(ptpn1_snps$gen), as.vector(ptpn1_snps$snp), c(rep("Chr20", nrow(ptpn1_snps))), c("50515979", "50516743", "50528132", "50573995", "50578711")))
colnames(ptpn1_position)[1] <- "gene"
colnames(ptpn1_position)[2] <- "snp"
colnames(ptpn1_position)[3] <- "chr"
colnames(ptpn1_position)[4] <- "pos" #this data frame can be used for create setup objects with position
ptpn1_position = ptpn1_position[order(ptpn1_position$pos),] #order snps by position in chr
#check
ptpn1_position$pos == ptpn1_position[order(ptpn1_position$pos),]$pos

#save it
write.table(ptpn1_position, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/tables/ptpn1_snps_ordered.csv", sep=",", col.names = TRUE, row.names = FALSE)

## select phenotypes of interest
pheno_selected = c(
    "CRF_sex",
    "CRF_age",
    "center",
    "CRF_weight",
    "CRF_height",
    "CRF_trici",
    "CRF_subscap",
    "MVPA_mean",
    "PA_factor",
    "energy_intake",    
    "mean_Alspac_v42",
    "CVi_softdrink_cont_2000",
    "obesity",
    "CRF_waist",
    "CRF_hip",    
    "waist_height",
    "waist_hip",
    "CRF_BMI",
    "CRF_Body_fat_PC",
    "FMI",                
    "CVi_BP",
    "SBP",
    "DBP",
    "TG",
    "TC",
    "LDL",
    "HDL",
    "LDL_HDL",
    "Apo_A1",
    "Apo_B",
    "ApoB_ApoA1",
    "apoB_LDL",
    "TG_HDL",
    "Insulin",
    "Leptin_ng_ml",
    "HOMA",
    "QUICKI", 
    "TC_HDL", 
    "risk_score",
    "risk_score_for_log")
cols_pheno = NULL
for (i in 1:length(pheno_selected)){
    pheno = pheno_selected[i]
    cols_pheno = append(cols_pheno, which(colnames(helena)==pheno))
}

## select final phenotypes along snps (all or ucp subseted)
final_cols_total = c(col_snps, cols_pheno)
final_cols_target_snps = c(helena_7_cols, cols_pheno)

## Subset rows total helena: 
#remove rows that have NA for all SNPs or all phenotypes in our set of interest snps
#select snps columns
helena_snps_subset = helena[,col_snps]
#subset pheno columns
helena_pheno_subset = helena[,cols_pheno]
#which rows have all NAs for all SNPs or all phenotypes in our set of interest snps
rows_to_remove_total = which(rowSums(is.na(helena_snps_subset)) == ncol(helena_snps_subset) | rowSums(is.na(helena_pheno_subset)) == ncol(helena_pheno_subset)) #if you ask for NAs and sum the result by row, you will obtain the number of cases with NA for each row. If the sum per row is equal to the number of columns, a give row will have NAs for all phenotypes or SNPs

#if we have rows without SNPs or phenotypes
if(length(rows_to_remove_total) != 0){

    #remove these rows
    helena_total = helena[-rows_to_remove_total,]
}

## Subset columns total helena
helena_total = helena_total[,final_cols_total]

#interest snps
#subset snps columns
helena_7_snps_subset = helena[,which(colnames(helena) %in% ptpn1_snps$snp)]
#subset pheno columns
helena_7_pheno_subset = helena[,which(colnames(helena) %in% pheno_selected)]
#which rows have all NAs for all SNPs or all phenotypes in our set of interest snps
rows_to_remove_ptpn1 = which(rowSums(is.na(helena_7_snps_subset)) == ncol(helena_7_snps_subset) | rowSums(is.na(helena_7_pheno_subset)) == ncol(helena_7_pheno_subset)) #if you ask for NAs and sum the result by row, you will obtain the number of cases with NA for each row. If the sum per row is equal to the number of columns, a give row will have NAs for all phenotypes or SNPs

#if we have rows without SNPs or phenotypes
if(length(rows_to_remove_ptpn1) != 0){

    #remove these rows
    helena_7 = helena[-rows_to_remove_ptpn1,]
}

## Subset columns helena subset of target snps
helena_7 = helena_7[,final_cols_target_snps]

#check that we have the 1057 individuals and the columns with snps and phenotypes
nrow(helena_total) == 1057
nrow(helena_7) == 1057
ncol(helena_total) == length(col_snps) + length(cols_pheno)
ncol(helena_7) == length(ptpn1_snps$snp) + length(pheno_selected)

## select SNPs columns in the subseted data.frames
# all snps
col_first_snp_final = which(colnames(helena_total) == "rs1537516")
col_second_snp_final = which(colnames(helena_total) == "seqSNP9")
col_snps_final = col_first_snp_final:col_second_snp_final

#only target snps
final_helena_7_cols = which(colnames(helena_7) %in% ptpn1_snps$snp)
length(final_helena_7_cols) == nrow(ptpn1_position)

## create a setupSNP object 
myData<-setupSNP(data=helena_total, colSNPs=col_snps_final, sep="")
summary(myData)
labels(myData)

## create a setupSNP object with only ucp SNPs
myData_ptpn1<-setupSNP(data=helena_7, colSNPs=final_helena_7_cols, sep="")
summary(myData_ptpn1)
labels(myData_ptpn1)

#save summary
resume_snps = summary(myData_ptpn1)
resume_snps

#search for snps with MAF lower than 0.9
snps_maf_low_0.9 = row.names(resume_snps[which(!resume_snps$major.allele.freq > 90),])#snps names

#check that with MAF>90 are not included in snps_maf_low_0.9
!row.names(resume_snps[which(resume_snps$major.allele.freq > 90),]) %in% snps_maf_low_0.9#ALL TRUE

#see snps with MAF higher 0.9
snps_maf_high_0.9 = row.names(resume_snps[which(resume_snps$major.allele.freq > 90),])
snps_maf_high_0.9
#check that these snps haf MAF higher 90
row.names(resume_snps[which(resume_snps$major.allele.freq > 90),]) %in% snps_maf_high_0.9


##quick check about how the function snp create the SNP variables

#load the required package
require(stringi)

#check that a SNP is different in the original dataframe of ptpn1 and in the SNPassoc object
head(myData_ptpn1$rs10485614)
head(helena_7$rs10485614)

#save the SNP from the original dataframe in a vector
new_geno_check = helena_7$rs10485614

#for each element of that vector
for(i in 1:length(new_geno_check)){

    #select the [i] element of the vector
    stri_sub(str=new_geno_check[i], from=2, to=1) <- "/"
        #add a slash in the middle of the entry, that is, in the middle of the genotype.
        #https://stackoverflow.com/questions/13863599/insert-a-character-at-a-specific-location-in-a-string
}

#check that the new SNP variable is exactly the same than the one created by SNPassoc
summary(new_geno_check == myData_ptpn1$rs10485614)

#detach the used package
detach("package:stringi")




#################################
######## Missing data ###########
#################################
plotMissing(myData_ptpn1) 
dev.off() #We may have a look at the information we have for each SNPs using plotMissing function. This function requires the data to be an object of class setupSNP. The top plot in Figure 2 shows the missing information for SNPs data set.




########################################
############ Genotyping rate ###########
########################################
mean_genotyping_rate = mean(100 - resume_snps$missing)

##selceted snps with genotyping rate equal or lower than 90%
snps_geno_rate_low = row.names(resume_snps[which(resume_snps$missing >= 10),])

#drop these snps
if(!length(snps_geno_rate_low) == 0){
    myData_ptpn1 = myData_ptpn1[,-which(colnames(myData_ptpn1) %in% snps_geno_rate_low)]

    #check
    !snps_geno_rate_low %in% colnames(myData_ptpn1) #all TRUE
}

#check that these snps were removed
length(which(labels(myData_ptpn1) %in% snps_geno_rate_low)) == 0



### Hardy-Weinberg equilibrium (HWE) 
## across groups
res <- tableHWE(myData_ptpn1)
res #The column indicated by flag shows those SNPs that are statistically significant at level 0.05. This significance level may be changed using the argument sig in the print function (e.g. print(myData_ptpn1, sig=0.0001)). The number of decimals may also be changed using the digits parameter.

## change significance level to one less conservative.
print(res, sig=0.01) #Only 0 with signals of deviation. It is a limit rasonable ("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/pdf/ukmss-33586.pdf"), see for example ("http://www.sciencedirect.com/science/article/pii/S0006322307006087?via%3Dihub#sec5")

HWE_pvals = as.vector(summary(myData_ptpn1)$HWE) #Extracted from here becuase tableHWE gives a SNPs with p.value = 0 as NA (rs1050450). 

#snps with deviation
snps_no_HW = row.names(summary(myData_ptpn1))[which(HWE_pvals < 0.01)]

#drop these snps for both data.sets
if(!length(snps_no_HW) == 0){
    myData_ptpn1 = myData_ptpn1[,-which(colnames(myData_ptpn1) %in% snps_no_HW)]

    #check
    !snps_no_HW %in% colnames(myData_ptpn1) #all false
}

#check that these snps were removed
length(which(labels(myData_ptpn1) %in% snps_no_HW)) == 0

#check that no snp remaining SNP have HWE < 0.01
length(which(summary(myData_ptpn1)$HWE < 0.01)) == 0

##A stratified analysis by
#sex
tableHWE(myData_ptpn1,strata=myData_ptpn1$CRF_sex)
#obesity
tableHWE(myData_ptpn1,strata=myData_ptpn1$obesity)
#CVi_BP
tableHWE(myData_ptpn1,strata=myData_ptpn1$CVi_BP)



### Final SNPs
final_snps = snps_labels[which(snps_labels$new_labels %in% colnames(myData_ptpn1)),]$new_labels
length(final_snps)

#snps not overlapped between snps_no_HW and snps_geno_rate_low
total_snps_removed = c(setdiff(snps_no_HW, snps_geno_rate_low), setdiff(snps_geno_rate_low, snps_no_HW))# you have to repeat setdiff two times with different order to include not overlapped elements of both vectors

#total snps removed less total snps gives the final number of snps
length(final_helena_7_cols) - length(total_snps_removed) == length(final_snps) 

## check that there are not SNPs with deviation form HWE or genotyping rate equal or lower than 90%
summary(!summary(myData_ptpn1)$missing >= 10)#all TRUE
summary(!summary(myData_ptpn1)$HWE < 0.01)#ALL TRUE

## plot missing
plotMissing(myData_ptpn1) 
dev.off()



###################################################################
######## COMPARE ALLELES HELENA VS 1000 GENOMES PROJECT ###########
###################################################################

#check minor alleles in HELENA respect to 1000 Genomes Project for the Europeans populations. We should have similar allele frequencies. 


## load the alleles names
alleles = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/alleles_ptpn1_v1.csv", sep=",", header=T)
nrow(alleles) == length(labels(myData_ptpn1))
    #First allele is major, and second is minor.
    #IMPORTANT NOTE: I have matched the allele names of HELENA and ncbi (at 10/09/2019). Now, they all matched, but I have noted that UCP alleles have changed in ncbi. For example, alleles that I changed to match ncbi, now are in ncbi exactly as original HELENA. If you check this for these SNPs, and you see changes no panic! The important thing is you are using the complementary chain and the first allele is always the major. Indeed, in many cases both options (i.e., both chains) are included as synonimous in ncbi (i.e., HGVS).
    #UPDATE APRIL 2021: One of the alleles seems to be wrong. From now on, you must always check that the minor allele in HELENA is the allele with the lowest frequency in 1000 Genomes Project for Europeans.


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


## change allele names IF NEEDED
#In our table, rs6067472 has T as the major and A as the minor allele. However, this is the opposite according to ncbi. T=0.3648. We have to exchange alleles A/T instead of T/A.

#I usually compare the alleles of each SNP with the alleles in HELENA and in ncbi. We have a word file for HELENA ("Appendix list of SNPs genotyped by GoldenGate[1].doc") that shows the major and minor allele of each SNP in HELENA. The first allele is the major, while the second allele is the minor. If for example, we have G and A as alleles in HELENA, while the same SNP in ncbi has C and T, I would understand that we are reporting the opposite strand in HELENA, so I would switch to the other strand to match the alleles of ncbi.

#In the case of rs6067472, the alleles are T and A, so this is a palindromic SNP. We have T and A as alleles in HELENA, while in ncbi we also have T and A, this leaded me to an error. I assumed that we are reporting the same strand in HELENA, but if we check the allele frequencies we can see they do not match. T is the major allele and A is the minor allele in HELENA, but ncbi shows T as the minor and A as the major according to the 1000 Genomes Project for European populations. Therefore, we have to also switch the alleles of this SNP to match the strand of ncbi.

#For the future, you should ALWAYS compare, not only the allele names, but also the allele frequencies between HELENA and ncbi. If we are reporting the same strand, the same allele should be the major, if we are reporting the opposite strand, the complementary allele of the major in HELENA should be the major in ncbi. If the SNP is palindromic (A/T - C/G), the same major allele would indicate that we are reporting the same strand, if not, we are reporting the opposite strand.

#Therefore, in this case, HELENA alleles should be T/A, the 1 (A) is less frequent, while 2 (T) is more frequent. In the case of ncbi, it should be the opposite: A/T. A is the major allele, while T is the minor.
summary(myData_ptpn1$rs6067472)

#add first the new combination of alleles in the ncbi columns
alleles$ncbi = factor(alleles$ncbi, levels=c(levels(alleles$ncbi), "A/T"))

#now change the allele order for the problematic SNP.
alleles[which(alleles$snp == "rs6067472"),]$ncbi <- "A/T"

#remove unused levels
alleles$ncbi = droplevels(alleles$ncbi)

#save the table with alleles names
write.table(alleles, "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/alleles_ptpn1_v2.csv", sep=",", col.names=TRUE, row.names=FALSE)



### Linkage Disequilibrium
#D prima cercano a 1 es LD, si ademas el pvalue es significativo entonces LD es signifciatvo: lo observado es diferente de lo esperado si hubiera random associations between alleles. See for more information: 
    #"http://pbgworks.org/sites/pbgworks.org/files/measuresoflinkagedisequilibrium-111119214123-phpapp01_0.pdf"
    #https://www.researchgate.net/post/I_have_Linkage_Disequilibrium_LD_data_for_two_SNPs-r2_is_about_014_D_is_around_08_Could_these_SNPs_be_said_to_be_in_strong_LD
require(genetics)

#Sselect snps
ptpn1_snp = myData_ptpn1[,which(colnames(myData_ptpn1) %in% labels(myData_ptpn1))]
#check that both objects have the same column order to used then labels(myData_positive_new) as columns names for the final dataset with genotype data
colnames(ptpn1_snp) == labels(myData_ptpn1)

#transform genotyping data of each snp to "genetics" format and then save it
ptpn1_snp_geno = data.frame(rep(NA, nrow(myData_ptpn1)))
for (i in 1:ncol(ptpn1_snp)){
    snp_selected = ptpn1_snp[,i]
    ptpn1_snp_geno = cbind(ptpn1_snp_geno, genotype(snp_selected))
}
ptpn1_snp_geno[,1] <- NULL
names(ptpn1_snp_geno) <- labels(myData_ptpn1)
str(ptpn1_snp_geno)

#LD between all snps
linkage_target_snps = LD(ptpn1_snp_geno)

#summary
summary(linkage_target_snps)
plot(linkage_target_snps, which="D'")
LDplot(linkage_target_snps, which="D'")
dev.off()

#plot LD
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/linkage_dis/plot_LD.pdf", width=20, height=20)
LDtable(linkage_target_snps, digits=2, colorize="D'", colorcut=seq(0,1, 0.1), colors=colorRampPalette(c("yellow", "red"))(length(seq(0,1, 0.1))), cex=0.8)
dev.off()




#########################################################################
################### Previous functions written ##########################
#########################################################################

#to make normality plots and test of normality under several transformations
normality_test_log = function(phenotype, data){

    variable = data[,which(colnames(data) == phenotype)]
    variable_name = colnames(data)[which(colnames(data) == phenotype)]

    #plots
    par(mfcol=c(2,2))
    plot(density(na.omit(variable)), main=paste("raw", variable_name))
    qqnorm((variable))
    qqline((variable), col = "red")        
    plot(density(log(na.omit(variable+1))), main=paste("Log", variable_name))
    qqnorm(log(variable+1))
    qqline(log(variable+1), col = "red")

    #normality tests
    require(nortest)

    #tests for normality with and without log transformation
    ad_log = ad.test(log(variable+1))
    ad = ad.test((variable))    
    cvm_log = cvm.test(log(variable+1))
    cvm = cvm.test((variable))    
    ks_log = lillie.test(log(variable+1))
    ks = lillie.test((variable))    
    pear_log = pearson.test(log(variable+1))
    pear = pearson.test((variable))    
    shap_log = sf.test(log(variable+1))
    shap = sf.test((variable))    

    #bind all pvalues of the tests
    significations = rbind.data.frame(cbind(ad$p.value, ad_log$p.value), cbind(cvm$p.value, cvm_log$p.value), cbind(ks$p.value, ks_log$p.value), cbind(pear$p.value, pear_log$p.value), cbind(shap$p.value, shap_log$p.value))
    colnames(significations)[1] <- "raw_variable"
    colnames(significations)[2] <- "log_transformed"

    print("###############################################")
    print(paste("#### P.values of", variable_name, "######"))
    print("###############################################")
    print(significations)
}
normality_test_sqrt = function(phenotype, data){

    variable = data[,which(colnames(data) == phenotype)]
    variable_name = colnames(data)[which(colnames(data) == phenotype)]

    #plots
    par(mfcol=c(2,2))
    plot(density(na.omit(variable)), main=paste("raw", variable_name))
    qqnorm((variable))
    qqline((variable), col = "red")        
    plot(density(sqrt(na.omit(variable))), main=paste("sqrt", variable_name))
    qqnorm(sqrt(variable))
    qqline(sqrt(variable), col = "red")

    #normality tests
    require(nortest)

    #tests for normality with and without sqrt transformation
    ad_sqrt = ad.test(sqrt(variable))
    ad = ad.test((variable))    
    cvm_sqrt = cvm.test(sqrt(variable))
    cvm = cvm.test((variable))    
    ks_sqrt = lillie.test(sqrt(variable))
    ks = lillie.test((variable))    
    pear_sqrt = pearson.test(sqrt(variable))
    pear = pearson.test((variable))    
    shap_sqrt = sf.test(sqrt(variable))
    shap = sf.test((variable))    

    #bind all pvalues of the tests
    significations = rbind.data.frame(cbind(ad$p.value, ad_sqrt$p.value), cbind(cvm$p.value, cvm_sqrt$p.value), cbind(ks$p.value, ks_sqrt$p.value), cbind(pear$p.value, pear_sqrt$p.value), cbind(shap$p.value, shap_sqrt$p.value))
    colnames(significations)[1] <- "raw_variable"
    colnames(significations)[2] <- "sqrt_transformed"

    print("###############################################")
    print(paste("#### P.values of", variable_name, "######"))
    print("###############################################")
    print(significations)
}
normality_test_squared = function(phenotype, data){

    variable = data[,which(colnames(data) == phenotype)]
    variable_name = colnames(data)[which(colnames(data) == phenotype)]

    #plots
    par(mfcol=c(2,2))
    plot(density(na.omit(variable)), main=paste("raw", variable_name))
    qqnorm((variable))
    qqline((variable), col = "red")        
    plot(density((na.omit(variable))^2), main=paste("squared", variable_name))
    qqnorm((variable)^2)
    qqline((variable)^2, col = "red")

    #normality tests
    require(nortest)

    #tests for normality with and without squared transformation
    ad_squared = ad.test((variable)^2)
    ad = ad.test((variable))    
    cvm_squared = cvm.test((variable)^2)
    cvm = cvm.test((variable))    
    ks_squared = lillie.test((variable)^2)
    ks = lillie.test((variable))    
    pear_squared = pearson.test((variable)^2)
    pear = pearson.test((variable))    
    shap_squared = sf.test((variable)^2)
    shap = sf.test((variable))    

    #bind all pvalues of the tests
    significations = rbind.data.frame(cbind(ad$p.value, ad_squared$p.value), cbind(cvm$p.value, cvm_squared$p.value), cbind(ks$p.value, ks_squared$p.value), cbind(pear$p.value, pear_squared$p.value), cbind(shap$p.value, shap_squared$p.value))
    colnames(significations)[1] <- "raw_variable"
    colnames(significations)[2] <- "squared_transformed"

    print("###############################################")
    print(paste("#### P.values of", variable_name, "######"))
    print("###############################################")
    print(significations)
}




#########################################################################
######################### Normality #####################################
#########################################################################
phenotypes_normality = c(
    "CRF_waist",
    "CRF_hip",    
    "waist_height",
    "waist_hip",
    "CRF_BMI",
    "CRF_Body_fat_PC",
    "FMI",    
    "SBP",
    "DBP",
    "TG",
    "TC",
    "LDL",
    "HDL",
    "LDL_HDL",
    "Apo_A1",
    "Apo_B",
    "ApoB_ApoA1",
    "apoB_LDL",
    "TG_HDL",
    "Insulin",
    "Leptin_ng_ml",
    "HOMA",
    "QUICKI",
    "TC_HDL",    
    "risk_score_for_log")

#### plot the distribution and qqplot of response variables
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/normality_plots_log.pdf")
for (i in 1:length(phenotypes_normality)){
    pheno = phenotypes_normality[i]
    normality_test_log(pheno, data=myData_ptpn1)
}
dev.off()

#### plot the distribution and qqplot of response variables
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/normality_plots_sqrt.pdf")
for (i in 1:length(phenotypes_normality)){
    pheno = phenotypes_normality[i]
    normality_test_sqrt(pheno, data=myData_ptpn1)
}
dev.off()

#### plot the distribution and qqplot of response variables
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/normality_plots_squared.pdf")
for (i in 1:length(phenotypes_normality)){
    pheno = phenotypes_normality[i]
    normality_test_squared(pheno, data=myData_ptpn1)
}
dev.off()

#### plot the distribution and qqplot of residuals of all models (used with continuous response variables)
type_heritage = c("recessive", "dominant", "overdominant", "additive", "codominant")
require(car)

#without transformation
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/residual_models_without_trans_plots.pdf")
#for each phenotype
for (i in 1:length(phenotypes_normality)){
    
    #select the [i] phenotype
    pheno = phenotypes_normality[i]
    
    #for each snp
    for(j in 1:length(labels(myData_ptpn1))){

        #select the [j] snp
        snp = labels(myData_ptpn1)[j]

        #for each heritage model
        for(k in 1:length(type_heritage)){
            
            #select the [k] model
            heritage = type_heritage[k]

            #extract number of snp levels
            n_levels_snp = length(unique(eval(parse(text=paste(heritage, "(na.omit(myData_ptpn1$", snp, "))", sep="")))))

            #if we have more than 1 level
            if(n_levels_snp > 1){

                #fit a model. the script is made for oth obesity and CV variables, in the latter case we include CRF_BMI as predictor
                if (pheno %in% c("CRF_waist", "CRF_hip", "waist_height", "waist_hip", "CRF_BMI", "CRF_Body_fat_PC", "FMI")){  
                    model = glm(as.formula(paste(pheno, "~", heritage, "(", snp, ")", "+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                } else {
                    model = glm(as.formula(paste(pheno, "~", heritage, "(", snp, ")", "+CRF_BMI+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                }  

                #plot residuals
                par(mfcol=c(2,2))
                Res = residuals(model, type="pearson")
                Fit <- fitted(model)
                plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted")
                abline(h=0) 
                plot(density(Res), main=paste("res", pheno, snp, heritage))
                shapiro.test(Res) 
                qqnorm(Res)
                qqline(Res)
            }    
        }    
    }
}
dev.off()

#with log transformation
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/residual_models_log_plots.pdf")
#for each phenotype
for (i in 1:length(phenotypes_normality)){
      
    #select the [i] phenotype
    pheno = phenotypes_normality[i]
    
    #for each snp
    for(j in 1:length(labels(myData_ptpn1))){

        #select the [j] snp
        snp = labels(myData_ptpn1)[j]

        #for each heritage model
        for(k in 1:length(type_heritage)){
            
            #select the [k] model
            heritage = type_heritage[k]

            #extract number of snp levels
            n_levels_snp = length(unique(eval(parse(text=paste(heritage, "(na.omit(myData_ptpn1$", snp, "))", sep="")))))

            #if we have more than 1 level
            if(n_levels_snp > 1){

                #fit a model. the script is made for oth obesity and CV variables, in the latter case we include CRF_BMI as predictor
                if (pheno %in% c("CRF_waist", "CRF_hip", "waist_height", "waist_hip", "CRF_BMI", "CRF_Body_fat_PC", "FMI")){  
                    model = glm(as.formula(paste("log(", pheno, ")", "~", heritage, "(", snp, ")", "+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                } else {
                    model = glm(as.formula(paste("log(", pheno, ")", "~", heritage, "(", snp, ")", "+CRF_BMI+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                }  

                #plot residuals
                par(mfcol=c(2,2))
                Res = residuals(model, type="pearson")
                Fit <- fitted(model)
                plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted")
                abline(h=0) 
                plot(density(Res), main=paste("res", pheno, snp, heritage))
                shapiro.test(Res) 
                qqnorm(Res)
                qqline(Res)
            }    
        }    
    }
}
dev.off()

#with sqrt transformation
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/residual_models_sqrt_plots.pdf")
#for each phenotype
for (i in 1:length(phenotypes_normality)){
      
    #select the [i] phenotype
    pheno = phenotypes_normality[i]
    
    #for each snp
    for(j in 1:length(labels(myData_ptpn1))){

        #select the [j] snp
        snp = labels(myData_ptpn1)[j]

        #for each heritage model
        for(k in 1:length(type_heritage)){
            
            #select the [k] model
            heritage = type_heritage[k]

            #extract number of snp levels
            n_levels_snp = length(unique(eval(parse(text=paste(heritage, "(na.omit(myData_ptpn1$", snp, "))", sep="")))))

            #if we have more than 1 level
            if(n_levels_snp > 1){

                #fit a model. the script is made for oth obesity and CV variables, in the latter case we include CRF_BMI as predictor
                if (pheno %in% c("CRF_waist", "CRF_hip", "waist_height", "waist_hip", "CRF_BMI", "CRF_Body_fat_PC", "FMI")){  
                    model = glm(as.formula(paste("sqrt(", pheno, ")", "~", heritage, "(", snp, ")", "+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                } else {
                    model = glm(as.formula(paste("sqrt(", pheno, ")", "~", heritage, "(", snp, ")", "+CRF_BMI+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                }  

                #plot residuals
                par(mfcol=c(2,2))
                Res = residuals(model, type="pearson")
                Fit <- fitted(model)
                plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted")
                abline(h=0) 
                plot(density(Res), main=paste("res", pheno, snp, heritage))
                shapiro.test(Res) 
                qqnorm(Res)
                qqline(Res)
            }    
        }    
    }
}
dev.off()

#with squares transformation
pdf("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/normality_phenotypes/residual_models_squared_plots.pdf")
#for each phenotype
for (i in 1:length(phenotypes_normality)){
      
    #select the [i] phenotype
    pheno = phenotypes_normality[i]
    
    #for each snp
    for(j in 1:length(labels(myData_ptpn1))){

        #select the [j] snp
        snp = labels(myData_ptpn1)[j]

        #for each heritage model
        for(k in 1:length(type_heritage)){
            
            #select the [k] model
            heritage = type_heritage[k]

            #extract number of snp levels
            n_levels_snp = length(unique(eval(parse(text=paste(heritage, "(na.omit(myData_ptpn1$", snp, "))", sep="")))))

            #if we have more than 1 level
            if(n_levels_snp > 1){

                #fit a model. the script is made for oth obesity and CV variables, in the latter case we include CRF_BMI as predictor
                if (pheno %in% c("CRF_waist", "CRF_hip", "waist_height", "waist_hip", "CRF_BMI", "CRF_Body_fat_PC", "FMI")){  
                    model = glm(as.formula(paste("(", pheno, ")^2", "~", heritage, "(", snp, ")", "+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                } else {
                    model = glm(as.formula(paste("(", pheno, ")^2", "~", heritage, "(", snp, ")", "+CRF_BMI+CRF_sex+CRF_age+center", sep="")), data=myData_ptpn1, family="gaussian")
                }  

                #plot residuals
                par(mfcol=c(2,2))
                Res = residuals(model, type="pearson")
                Fit <- fitted(model)
                plot(Res ~ Fit, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted")
                abline(h=0) 
                plot(density(Res), main=paste("res", pheno, snp, heritage))
                shapiro.test(Res) 
                qqnorm(Res)
                qqline(Res)
            }    
        }    
    }
}
dev.off()




####################################################################
##### Gene - phenotype association with biochemical variables  #####
####################################################################

#NOTE: All SNPs are used with dominant model (including those with MAF < 0.1). The rest of models are using in SNPs with MAF > 0.1.

#Note: We run additive model, which is a type of codominant model in which the heterozygous have an intermediate value of the phenotype between the other two genotypes. A codominant model assumes a different phenotype for heterozygous, but it has not to be exactly intermediate. See "http://www.bio.net/mm/gen-link/1998-January/001415.html" for further information. 

#We will also use this information to obtain a supplementary dataset with as many rows as phenotype*model*SNP combinations, showing for each one the phenotype, polymorphism, heritage model, minimum sample size, p-value, false discovery rate and R2 for the gene-phenotype association.

#We will apply a loop to run the models of each association. From that model, we will obtain: phenotype, polymorphism, heritage model, minimum sample size (the group with less individuals), p-value, false discovery rate and R2. Respect to R2: We can calculate the R2 (adjusted? Zimmerman and Guisan function?) of the model including the SNP. Maybe, we can calculate R2 from the likelihood ratio test, which would give us the exact R2 of the SNP after accounting for the rest of confounding factors.

#required packages for this
Dsquared_mod = function(model, null_model, adjust = TRUE, n_levels_predictor_test=NULL, predictor_is_factor=NULL) {
    # version 1.1 (13 Aug 2013)
    # calculates the explained deviance of a likelihood ratio test between a complex model and a nested one without one predictor. I MADE THIS MODIFICATION. 
    # model: a model object of class "glm"
    # the null model which is equal than model but without 1 predictor. The predictor removed will be the predictor under study  
    # adjust: logical, whether or not to use the adjusted deviance taking into account the number of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000). More data used for fitting enhance d2, whilst more parameters lower it.
    # n_levels_predictor_test: set the number of the levels of the predictor that is removed in the null model to calculate the adjusted R2.
    d2 <- (null_model$deviance - model$deviance) / null_model$deviance #If the model with the predictor of interest (SNP) has the same deviance (desviacion) than the null model (exactly equal to the model but without the SNP), then model (and hence the predictor) does not reduce the deviance and consequently R2 = 0 as null_model$deviance - model$deviance = 0. As the deviance in the model decreases respect to the null model, the denominator is higher and hence R2 is higher. 
    if (adjust) { #adjust by the sample size and the number of parameters tested. PROBLEM HERE: We are setting here as the number of parameters, the total number of predictors in the complex model, when the null model share all of them except 1. Indeed the change in deviance is for the change of only 1 predictor, so putting as number of coefficients ALL could be too stringent. You are taking the change of deviance calculated in a likelihood ratio test between two models. The difference in degrees of freedom between models is the number of levels of the factor removed (SNP in this case) less 1. So that would be to set p as the number of levels of the SNP less 1, as you have one coefficient for each level respect to the reference level. If the predictor changed between models is continuous, the number of levels is 1.
        n <- length(model$fitted.values)
        #p <- length(model$coefficients)
        p <- ifelse(predictor_is_factor, n_levels_predictor_test - 1, n_levels_predictor_test) #if the predictor is a factor we set the number of coefficients as the number of levels less 1. If the predictor is not a factor, then we set as just the number of levels that for a continuous variable should be 1.
        #IMPORTANT: An alternative for obtaining "p" would be to extract the "df" of the likelihood ratio tests between the models: anova(model, null_model, test="Chi")$Df[2]. Note that the change in df between models is 1 for each continuous predictor, for factors is n_levels minus 1. This is similar to the number of coefficients. If you have a three-level factor you have 2 coefficients respect to the reference level (3-1=2). 
        d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2) #higher R2 gives smaller (1-d2) and hence less is subtracted from 1, leading to higher R2. Lower n (lower sample size) and higher p (more parameters) leads to smaller numerator and smaller denominator being the result higher and then the result of multiplying by (1-d2) would be bigger, subtracting more from 1 and giving less R2. 
    }
    return(d2)
} #D2 function of Nick. If you have problem with this you can use modEvA ("http://modeva.r-forge.r-project.org/") or MuMIn ("https://cran.r-project.org/web/packages/MuMIn/index.html"). 
require(MuMIn) #for calculating the R2 from a likelihood ratio test. This will be compared with Dsquared_mod

#pheno to model
pheno_to_model = c("obesity","CRF_waist","CRF_hip","waist_height","waist_hip","CRF_BMI","CRF_Body_fat_PC","FMI","TC","LDL","HDL","TC_HDL","LDL_HDL","TG","TG_HDL","Apo_A1","Apo_B","ApoB_ApoA1","apoB_LDL","Insulin","Leptin_ng_ml","HOMA","QUICKI","SBP","DBP","risk_score_for_log")

#snps to test
snp_to_test = labels(myData_ptpn1)
length(snp_to_test) == nrow(ptpn1_snps)

#models to test
models = c("codominant", "dominant", "recessive", "overdominant", "additive")

#run ALL associations geno-pheno
geno_pheno_results = data.frame(selected_pheno=NA, selected_model=NA, snp_to_test=NA, min_n=NA, pvals=NA, BF=NA, fdr=NA, d2_mine=NA, adjust_d2_mine=NA, d2_mumin=NA, adjust_d2_mumin=NA)
for(p in 1:length(pheno_to_model)){

    #select the [p] phenotype
    selected_pheno = pheno_to_model[p]

    #select the family and transformation
    if(!selected_pheno %in% c("obesity", "CVi_BP")){
        family="gaussian"
        transformation="log("
    } else{
        family="binomial"
        transformation="("
    }

    #select control variables: BMI is a control variable only for CVD variables
    if(!selected_pheno %in% c("CVi_BP", "SBP", "DBP", "TG", "TC", "LDL", "HDL", "LDL_HDL", "Apo_A1", "Apo_B", "ApoB_ApoA1", "apoB_LDL", "TG_HDL", "Insulin", "Leptin_ng_ml", "HOMA", "QUICKI",  "TC_HDL",  "risk_score_for_log")){
        control_variables = "CRF_sex+CRF_age+center"
    } else{
        control_variables = "CRF_BMI+CRF_sex+CRF_age+center"         
    }

    #for each model
    for(m in 1:length(models)){

        #select the [m] model
        selected_model = models[m]
        
        #open empty vectors to save the results of all SNPs in the [m] heritage model
        d2_mine=NULL
        adjust_d2_mine=NULL
        d2_mumin=NULL
        adjust_d2_mumin=NULL
        min_n=NULL
        pvals = NULL
        
        #for each model               
        for(k in 1:length(snp_to_test)){
    
            #select the [k] snp
            selected_snp = snp_to_test[k]
    
            #extract genotype data of this SNP
            geno_data = eval(parse(text=paste(selected_model, "(na.omit(myData_ptpn1$", selected_snp, "))", sep="")))

            #extract genotype levels
            levels_genotypes = unique(geno_data)

            #number of levels of the [i] snp
            n_levels_snp = length(levels_genotypes)

            #calculate the number of individuals with each genotype and with data for the [p] phenotype
            sample_size_per_geno = NULL
            #for each genotype
            for(l in 1:length(levels_genotypes)){

                #select the [l] genotype
                selected_level = levels_genotypes[l]

                #if the variable is continuous
                if(!family=="binomial"){

                    #select those rows with the [l] level of the genotype and without na for the select phenotype
                    subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")=='", selected_level, "' & !is.na(myData_ptpn1$", selected_pheno, ")),]", sep="")))

                    #save the number of rows
                    sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))                    
                } else {#if not and then the variable is discrete
                    
                    #extract the levels of the phenotype
                    levels_factor = sort(unique(eval(parse(text=paste("myData_ptpn1$", selected_pheno, sep="")))))

                    #for each level of the discrete phenotype
                    for(f in 1:length(levels_factor)){

                        #select the [f] level of the factor
                        selected_level_factor = levels_factor[f]

                        #select those rows with the [l] level of the genotype and without na for the select phenotype
                        subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", selected_snp, ")=='", selected_level, "' & myData_ptpn1$", selected_pheno, "=='", selected_level_factor, "'),]", sep="")))
                        
                        #save the number of rows
                        sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))
                    }
                }
            }

            #run the models if: snp has MAF highen 0.9 and the model is dominant (these snps cannot be fitted wit the other models); if there are more than 10 individuals per level of genoype; if there is only ONE genotype for this SNP in the cohort
            if(selected_model %in% c("recessive", "overdominant", "codominant", "additive") & !selected_snp %in% snps_maf_low_0.9 | TRUE %in% (sample_size_per_geno < 10) | n_levels_snp==1){ #FIJAMOs un mínimo de 10 individuos por genotipo y con dato del correspondiente fenotipo para correr el modelo. Así reducimos (aunque no eliminamos) el riesgo de calcular un p.value muy significativo por un nivel con 2-3 individuos que sale muy diferente. Ese p.value luego afectaría al calculo del FDR. Además añadimos el filtro del número de niveles, porque podemos tener solo un genotipo en la cohorte para un snp dado, entonces lógicamente para ese nivel tenemos más de 10 individuos, pero es que solo hay un nivel, en ese caso no hay nada que hacer. Al menos que haya dos niveles.
                d2_mine=append(d2_mine, NA)
                adjust_d2_mine=append(adjust_d2_mine, NA)
                d2_mumin=append(d2_mumin, NA)
                adjust_d2_mumin=append(adjust_d2_mumin, NA)
                min_n=append(min_n, NA)                
                pvals = append(pvals, NA)
            } else{
    
                #run models with and without the snp
                model1 = glm(paste(transformation, selected_pheno, ") ~", selected_model, "(", selected_snp, ")+", control_variables, sep=""), data=eval(parse(text=paste("myData_ptpn1[which(!is.na(myData_ptpn1$", selected_snp, ")),]", sep=""))), family=family)
                model2 = glm(paste(transformation, selected_pheno, ") ~", control_variables, sep=""), data=eval(parse(text=paste("myData_ptpn1[which(!is.na(myData_ptpn1$", selected_snp, ")),]", sep=""))), family=family)#data is filtered by NAs in the snp, because in the second model we don't have the snp, so rows with NAs in that snp could be included and then the two models would be fitted with different data

                #extract the pvals                
                results = anova(model1, model2, test="Chi")$"Pr(>Chi)"[2]


                ##calculate the R2. WE WANT the R2 of the SNP after adjusting by the controlling variables.
                #first using the modified function of Nick Zimmermann
                d2_nick = Dsquared_mod(model=model1, null_model = model2, adjust=FALSE)
                adjust_d2_nick = Dsquared_mod(model=model1, null_model = model2, adjust=TRUE, n_levels_predictor_test=n_levels_snp, predictor_is_factor=TRUE) #we calculate the adjusted R2 and for that we need to set TRUE for adjust and consider the number of genotypes because that number will be the number of coefficients (one estimate for each level) that differ between the null and the complex model. We indicate that the predictors tested is a factor (the SNP) to calculate the the difference in coefficient between the null and the complex model (each level has an estimate respect to the reference level).

                #second with the MuMIn package. This statistic is is one of the several proposed pseudo-R^2's for nonlinear regression models. It is based on an improvement from _null_ (intercept only) model to the fitted model, and calculated as R^2 = 1 - exp(-2/n * logL(x) - logL(0)) where logL(x) and logL(0) are the log-likelihoods of the fitted and the _null_ model respectively. ML estimates are used if models have been fitted by REstricted ML (by calling ‘logLik’ with argument ‘REML = FALSE’). Note that the _null_ model can include the random factors of the original model, in which case the statistic represents the ‘variance explained’ by fixed effects. For OLS models the value is consistent with classical R^2. In some cases (e.g. in logistic regression), the maximum R_LR^2 is less than one.  The modification proposed by Nagelkerke (1991) adjusts the R_LR^2 to achieve 1 at its maximum: Radj^2 = R^2 / max(R^2) where max(R^2) = 1 - exp(2 / n * logL(0)) . ‘null.fit’ tries to guess the _null_ model call, given the provided fitted model object. This would be usually a ‘glm’. The function will give an error for an unrecognised class.
                d2_mumin_raw = r.squaredLR(object=model1, null=model2, null.RE=FALSE)
                    #object: a fitted model object (a.k.a. the full model)
                    #null: a fitted null model. This is the model to compare with. In my case, the simpler and nested model without the SNP under study.
                    #null.RE: logical, should the null model contain random factors?  Only used if no _null_ model is given, otherwise omitted, with a warning.
                    #the rest of arguments are for fitting the null model, but we have already fit it. 
                    #the result is: ‘r.squaredLR’ returns a value of R_LR^2, and the attribute ‘"adj.r.squared"’ gives the Nagelkerke's modified statistic. Note that this is not the same as nor equivalent to the classical ‘adjusted R squared’. IT IS NOT THE SAME ADJUST THAN THE ZIMMERMAN ADJUSTMENT
                #extract d2 and adjust d2 from the mumin results
                d2_mumin_pkg = d2_mumin_raw[1]
                adjust_d2_mumin_pkg = attributes(d2_mumin_raw)$adj.r.squared


                ##save the results
                d2_mine=append(d2_mine, d2_nick)
                adjust_d2_mine=append(adjust_d2_mine, adjust_d2_nick)
                d2_mumin=append(d2_mumin, d2_mumin_pkg)
                adjust_d2_mumin=append(adjust_d2_mumin, adjust_d2_mumin_pkg)
                min_n=append(min_n, min(sample_size_per_geno))
                pvals = append(pvals, results)
            }    
        }
        
        #calculate FDR
        fdr = p.adjust(pvals, method="BH")#correction for BH, it is valid use this correction when markerse are correlated directly, if two markers are correlated, a higher significance in one, will entail higher significance in the other. In the case the correlation was negative, we should used other correction, I think remember that Benjamini & Yekutieli (2001), but check. For further details see: "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_ucp_cv/p_value_correction/pvalue_correction.pdf"

        #calculate BF
        bf_correction = 0.05/length(na.omit(pvals))
        BF=NULL
        for(i in 1:length(pvals)){
            selected_pval = pvals[i]
            if(!is.na(selected_pval)){
                if(selected_pval < bf_correction){
                    BF=append(BF, "YES")
                } else {
                    BF=append(BF, "NO")
                }
            } else {
                BF=append(BF, NA)
            }    
        }

        #bind all
        final_results = cbind.data.frame(rep(selected_pheno, length(snp_to_test)), rep(selected_model, length(snp_to_test)), snp_to_test, min_n, pvals, BF, fdr, d2_mine, adjust_d2_mine, d2_mumin, adjust_d2_mumin)
        colnames(final_results)[1] <- "selected_pheno"
        colnames(final_results)[2] <- "selected_model"
        colnames(final_results)[3] <- "snp_to_test"
        colnames(final_results)[4] <- "min_n"
        colnames(final_results)[5] <- "pvals"
        colnames(final_results)[6] <- "BF"                
        colnames(final_results)[7] <- "fdr"
        colnames(final_results)[8] <- "d2_mine"
        colnames(final_results)[9] <- "adjust_d2_mine"
        colnames(final_results)[10] <- "d2_mumin"
        colnames(final_results)[11] <- "adjust_d2_mumin"

        #select significant results according to FDR<0.1       
        geno_pheno_results = rbind.data.frame(geno_pheno_results, final_results)
    }
}

#remove the first row with NAs
geno_pheno_results = geno_pheno_results[-which(rowSums(is.na(geno_pheno_results)) == ncol(geno_pheno_results)),]

#check we have all the associations
nrow(geno_pheno_results) == length(pheno_to_model) * length(snp_to_test) * length(models)


##Comparisons R2
#R2 mine and with mumin: They are very similar.
summary(geno_pheno_results$d2_mine - geno_pheno_results$d2_mumin) 
plot(geno_pheno_results$d2_mine, geno_pheno_results$d2_mumin)
cor.test(geno_pheno_results$d2_mine, geno_pheno_results$d2_mumin, method="spearman") #rho=0.999
dev.off()

#R2 adjusted is not equal with both methods.
summary(geno_pheno_results$adjust_d2_mine - geno_pheno_results$adjust_d2_mumin) 
plot(geno_pheno_results$adjust_d2_mine, geno_pheno_results$adjust_d2_mumin)
cor.test(geno_pheno_results$adjust_d2_mine, geno_pheno_results$adjust_d2_mumin, method="spearman") #rho=0.999
    #The non-adjusted values are similar between methods, but not when we adjust. Note that the method for adjusting used in the MuMIn function is different from the typical adjust. Mumin uses the Nagelkerke adjustment (see below). I have made the adjust manually modifying the function of Zimmermann. I set as the number of coefficients the number of genotypes less 1, because we have 1 coefficient for each level of the factor respect to the reference level. This is congruent with the fact that we are comparing the decrease in deviance between a model with the SNP and a simpler model with all the confounding factors but without the SNP. So we are checking the change in deviance caused by the SNP, and thus we have to consider the changes in degree of freedom (i.e., coefficients) caused by that SNP.
dev.off()


#For OLS models (linear regression?) the value is consistent with classical R^2. In some cases (e.g. in logistic regression), the maximum R_LR^2 is less than one.  The modification proposed by Nagelkerke (1991) adjusts the R_LR^2 to achieve 1 at its maximum: Radj^2 = R^2 / max(R^2). So you are basically adjusting the R2 by the maximum R2 present in your data. I see the point because the R2 max is not 1 in logistic. I see the point, you could underestimate the R2 for logistic models, but I do not fully understand how the maximum R2 can be established for a given model comparison.

#there are differences between adjusted and non-adjusted in mumin, but this is not caused by the logistic
plot(geno_pheno_results$adjust_d2_mumin, geno_pheno_results$d2_mumin, col="red")
cor.test(geno_pheno_results$adjust_d2_mumin, geno_pheno_results$d2_mumin)
dev.off()

#you can see here how when considering only obesity (factor with 2 levels, that is, logistic), there is a perfect correlation between adjusted and un-adjusted in mumin. 
plot(geno_pheno_results[which(geno_pheno_results$selected_pheno == "obesity"),]$adjust_d2_mumin, geno_pheno_results[which(geno_pheno_results$selected_pheno == "obesity"),]$d2_mumin, col="red")
cor.test(geno_pheno_results[which(geno_pheno_results$selected_pheno == "obesity"),]$adjust_d2_mumin, geno_pheno_results[which(geno_pheno_results$selected_pheno == "obesity"),]$d2_mumin, col="red")
dev.off()

#Given that I am not 100% sure if this a correct way to adjust and R2 does NOT seem underestimated for obesity (logistic), I am going to use the traditional pseudo R2 of glms, which is similar with the formula of Zimmermamnn and the function of MuMIn. In addition, my d2 and the mumin d2 are not so different respect to adjusted d2 using the traditional approach
summary(geno_pheno_results$d2_mine-geno_pheno_results$adjust_d2_mine)
plot(geno_pheno_results$d2_mine, geno_pheno_results$adjust_d2_mine)
summary(geno_pheno_results$d2_mumin-geno_pheno_results$adjust_d2_mine)
plot(geno_pheno_results$d2_mumin, geno_pheno_results$adjust_d2_mine)
    #There are differences, for each R2 value not adjusted, there is a very similar R2 adjusted value but also another value that differ in 0.001. My guess is that can be caused because recessive, dominant and overdominant models only have 2 levels, so 2-1 gives 1 parameter, thus the adjusted R2 is bigger and exactly similar to non-adjusted R2 (1 - ((n - 1) / (n - p)) * (1 - d2)). You can check it by calculating R2 adjusted and un-adjusted with Dsquared_mod, but setting the number of level to 2. The result is the same. 
plot(geno_pheno_results[which(!geno_pheno_results$selected_model %in% c("codominant", "additive")),]$d2_mumin, geno_pheno_results[which(!geno_pheno_results$selected_model %in% c("codominant", "additive")),]$adjust_d2_mine, col="red")
points(geno_pheno_results[which(geno_pheno_results$selected_model %in% c("codominant", "additive")),]$d2_mumin, geno_pheno_results[which(geno_pheno_results$selected_model %in% c("codominant", "additive")),]$adjust_d2_mine, col="blue")
    #the red dots are from models with two levels and have the similar R2 and adjusted R2. However, blue dots comes from additive and codominant model, being the adjusted R2 smaller, because now you have more parameters (3 genotypes - 1 = 2 levels).
dev.off()

#Given that both R2 (mine and mumin) are VERY similar and that mumin is more reproducible because it is in a published package, we will use the R2 of Mumin. 

#Summary: We will use the classical R^2 calculated with the mumin package.


## use these results to create the supplementary data 1
#set the folder to save
folder_to_save_supple_data_1 = "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/supple_data"
system(paste("mkdir -p ", folder_to_save_supple_data_1, sep=""))
    #p: no error if existing, make parent directories as needed

#select the rows and columns we are interested
suppl_data_1 = geno_pheno_results[which(geno_pheno_results$selected_pheno != "risk_score_for_log"), which(colnames(geno_pheno_results) %in% c("selected_pheno", "selected_model", "snp_to_test", "min_n", "pvals", "fdr", "d2_mumin"))]
    #We remove all rows belonging to the CVD risk score because this variable was finally not used in the manuscript. 
        #IMPORTANT: If we include risk_score_for_log, then we should change to risk_score in the supple data. 
    #We select the columns that includes the variables selected for the supplementary dataset 1

#convert R2 from 0-1 to percentage
suppl_data_1$d2_mumin = (suppl_data_1$d2_mumin*100)/1
    #For example: If over 1, we have 0.5, over 100 we would have X. X being (100*0.5)/1=50 -> 50% 

#change columns names
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "selected_pheno")] <- "phenotype"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "selected_model")] <- "heritage_model"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "snp_to_test")] <- "snp"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "min_n")] <- "min_sample_size"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "pvals")] <- "p_value"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "fdr")] <- "fdr"
colnames(suppl_data_1)[which(colnames(suppl_data_1) == "d2_mumin")] <- "r2_percentage"

#set the names of the files for saving
suppl_data_1_file_name = "suplementary_data_1_v2.csv"
suppl_data_1_file_name_zip = "suplementary_data_1_v2.zip"

#save the table
write.table(suppl_data_1, paste(folder_to_save_supple_data_1, "/", suppl_data_1_file_name, sep=""), col.names=TRUE, row.names=FALSE, sep=",")

#compress the text file and remove it after compression
system(paste("cd ", folder_to_save_supple_data_1, "; rm ", suppl_data_1_file_name_zip, "; zip ", suppl_data_1_file_name_zip, " ", suppl_data_1_file_name, " ; rm ", suppl_data_1_file_name, sep=""))
    #we could save directly as ".gz" using gzfile() around the file path with write.table. I avoid this option because this could give problems with the reviewers and readers, because they have to use a third party software to open it in windows.


## compare this version of the supplementary with the first one.

#Initially, I made some changes in this version for recoding the a problematic SNP respect allele names between HELENA and ncbi. But after thinking, I discover a more minimalistic way to solve this problem without touching the dataset.

#read the file of the previous version
suppl_data_1_previous_version = read.table(paste(folder_to_save_supple_data_1, "/suplementary_data_1_v1.txt.gz", sep=""), sep="\t", header=TRUE)

#read the file of the last version
suppl_data_1_last_version = read.table(unz(paste(folder_to_save_supple_data_1, "/suplementary_data_1_v2.zip", sep=""), suppl_data_1_file_name), sep=",", header=TRUE)

#select the non-numeric columns
no_numeric_cols = which(colnames(suppl_data_1_last_version) %in% c("phenotype", "heritage_model", "snp"))

#check that the non-numeric and numeric columns are the same between versions
summary(suppl_data_1_previous_version[,no_numeric_cols] == suppl_data_1_last_version[,no_numeric_cols])
summary(suppl_data_1_previous_version[,-no_numeric_cols] - suppl_data_1_last_version[,-no_numeric_cols])


##extract snps for further analyses
#FDR<0.1
geno_pheno_results_fdr_0.1 = geno_pheno_results[which(geno_pheno_results$fdr<0.1),]
geno_pheno_results_fdr_0.1
#check
summary(geno_pheno_results_fdr_0.1$fdr<0.1)
summary(geno_pheno_results[-which(geno_pheno_results$fdr<0.1),]$fdr>0.1)

#Combinations of phenotypes*snps that are significant for crude associations
significan_snp_pheno_crude = interaction(geno_pheno_results_fdr_0.1$selected_pheno, geno_pheno_results_fdr_0.1$snp_to_test)
    #this will be used for select haplotype and analyses
#check
summary(significan_snp_pheno_crude == paste(geno_pheno_results_fdr_0.1$selected_pheno, ".", geno_pheno_results_fdr_0.1$snp_to_test, sep=""))




###############################################
########### BLOQUES HAPLOTIPOS ################
###############################################

##### prepare data from haploview
##sex variable
sex = factor(myData_ptpn1$CRF_sex, levels = c(1,2,levels(myData_ptpn1$CRF_sex)))
sex[which(sex=="male")] <- 1
sex[which(sex=="female")] <- 2
sex = droplevels(sex)
summary(sex)
unique(myData_ptpn1[which(sex==1),]$CRF_sex)#all male
unique(myData_ptpn1[which(sex==2),]$CRF_sex)#all female

##create a matrix with samples as rows and markers as columns
#Los SNPs (columnas) tienen que estar ORDENADOS según su posición en el cromosoma.
# extract snp data form myData_ptpn1
geno_data = myData_ptpn1[,which(colnames(myData_ptpn1) %in% labels(myData_ptpn1))]
#extract snps orderded by position in the chromosome
snps_ordered = ptpn1_position[order(ptpn1_position$pos),]$snp
ptpn1_position$pos == sort(ptpn1_position$pos) #snps of ptpn1_position are ordered in basis on position
#reorder geno data in basis chrosome position
geno_data = geno_data[as.character(snps_ordered)]
str(geno_data)
#check order
colnames(geno_data) == snps_ordered #all TRUE
#loop for changing genotype data to code 1,2,3 
final_geno_data = data.frame(rep(NA, nrow(myData_ptpn1)))
for(i in 1:ncol(geno_data)){ #for each row
    
    #select the row
    seletec_col = geno_data[,i]

    #add three levels to the column (marker)
    seletec_col = factor(seletec_col, levels = c(1,2,3,levels(seletec_col)))
    
    #change the genotypes by numbers
    seletec_col[which(seletec_col=="1/1")] <- 1
    seletec_col[which(seletec_col=="1/2")] <- 3
    seletec_col[which(seletec_col=="2/2")] <- 2 #codification required by makeHaploviewInputFile

    #drop ancient levels
    seletec_col = droplevels(seletec_col)

    #check
    unique(eval(parse(text=paste("myData_ptpn1[which(seletec_col =='1'),]$", colnames(geno_data)[i], sep=""))))#all 1/1
    unique(eval(parse(text=paste("myData_ptpn1[which(seletec_col =='3'),]$", colnames(geno_data)[i], sep=""))))#all 1/2
    unique(eval(parse(text=paste("myData_ptpn1[which(seletec_col =='2'),]$", colnames(geno_data)[i], sep=""))))#all 2/2  

    #save new marker variable 
    final_geno_data = cbind.data.frame(final_geno_data, seletec_col)
}
#drop NA column
final_geno_data = final_geno_data[,-1]
#set column names
colnames(final_geno_data) <- colnames(geno_data) 
str(final_geno_data)

#separamos los analiss entre cromosomas no he visto ningún paper en el que analicen haplotipos de diferentes cromosomas juntos. Los bloques no cambian de todas maneras. Además, los algoritmos buscan listas de SNPs contiguos sin señales de recombinación, por tanto no tiene sentido que pongamos juntos el último SNP de UCP1 y el primero de UCP2 

#extract the chromosome of the snps implicated in significant associations
significant_chr = unique(chromosome_snps_helena[which(chromosome_snps_helena$selected_snp %in% unique(geno_pheno_results_fdr_0.1$snp_to_test)),]$chr_snps)#to what genes belongs snps associated with phenotypes according to FDR<0.1

#required package to prepare dta for haploview
require(HapEstXXR)

#for each gene we want to create datafile ready for haploview
for(i in 1:length(significant_chr)){

    #select the [i] gene
    selected_chr = paste("Chr", significant_chr[i], sep="")

    #selected the markers of the [i] gene      
    selected_markers = ptpn1_position[ptpn1_position$chr==selected_chr,]$snp

    #extract from geno_data the columns of the markers for the [i] gene
    final_geno_data_markers = final_geno_data[,which(colnames(final_geno_data) %in% selected_markers)]
    summary(colnames(final_geno_data_markers) == selected_markers)
    str(final_geno_data_markers)
    
    #check that the position of the columns is correct
    if(!FALSE %in% (colnames(final_geno_data_markers) == ptpn1_position[ptpn1_position$chr==selected_chr,]$snp)){
    
        #be sure you don't have a previous object called markers_pos
        rm(markers_pos)

        #if the position is correct, you can extract snp positions directly
        markers_pos = ptpn1_position[ptpn1_position$chr==selected_chr,]$pos
    } else {

        #be sure you don't have a previous object called markers_pos
        rm(markers_pos)

        #if not print error and stop
        stop("ERROR!!: colnames of final_geno_data_markers and snps of the dataframe with positions are not in the same order")        
    }

    #create haploviewdata
    makeHaploviewInputFile(famid = rep(NA, nrow(myData_ptpn1)), patid = row.names(myData_ptpn1), aff = rep(0, nrow(myData_ptpn1)), fid = rep(NA, nrow(myData_ptpn1)), mid = rep(NA, nrow(myData_ptpn1)), sex = sex, geno = final_geno_data_markers, marker.name=selected_markers, marker.position=markers_pos, haploview.pedfile = paste("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/haplotype/haplo_", selected_chr, ".pedfile", sep=""), haploview.infofile = paste("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/haplotype/haplo_", selected_chr, ".infofile", sep=""))
}

##open haploview
system("cd /media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/software/haploview; java -jar Haploview.jar", wait = FALSE) ## https://www.broadinstitute.org/haploview/haploview. If you want the plot with rs number, you have to include the ".infofile" in haploview. 

#COMPARAR REGLA DE LOS CUATRO GAMETOS Y CONFIDENCE INTERVAL Y TENER EN CUENTA SNPS CON LOW MAF. 

#RESULTS: 1 block according to the four gamete rule and the confidence interval (gabriel et al). 

#Los bloques de haplotipos se definen mediante el algoritmo de la regla de los cuatro gametos, el cual busca sets de SNPs contiguos y ordenados en los que no hay evidencias de recombinación, ya que considera que la recombinación se puede dar entre haplitpos, pero no dentro de haplotipos. See "Distribution of recombination crossovers and the origin of haplotype blocks: The interplay of population history, recombination, and mutation". Con los confidence intervals de gabriel et al sale igual. Lo he comprobado porque Si hay MAF muyy bajo no puedes usar la regla de los 4 gametos y tendría que suar el por defecto, porque te va amarcar cosas como blques por ser pocos frecuentes. En nuestro caso no hay problema. 

#Los valores de LD obtenidos de haploview son los mismos obtenidos con LD de genetics (solo diferencias de 1 decimal en algunos casos (e.g., 0.93 en vez de 0.94)). Los valores de MAF, genotyping sucess son identicos. Los p.values de HW son casi identicos. El orden de los SNPs es el correcto, posicion física. Por tanto, se están suando los datos correctos. 

#EN LA PESTAÑA HAPLOTYPES SE VEN LAS FRECUENCIAS DE LOS AHPLOTIPOS Y LA RECOMBINACIÓN ENTRE BLOQUES.

#los colores indicate support for linkage. Se basan en LOD y D' values: 1) LOD score (logarithm (base 10) of odds), developed by Newton Morton,[9] is a statistical test often used for linkage analysis in human, animal, and plant populations. The LOD score compares the likelihood of obtaining the test data if the two loci are indeed linked, to the likelihood of observing the same data purely by chance. Positive LOD scores favour the presence of linkage, whereas negative LOD scores indicate that linkage is less likely; 2) D': Suppose that among the gametes that are formed in a sexually reproducing population, allele A occurs with frequency p_{A}} p_{A} at one locus (i.e. p_{A}} p_{A} is the proportion of gametes with A at that locus), while at a different locus allele B occurs with frequency p_{B}} p_{B}. Similarly, let p_{AB}} p_{AB} be the frequency with which both A and B occur together in the same gamete (i.e. p_{AB}} p_{AB} is the frequency of the AB haplotype). The association between the alleles A and B can be regarded as completely random—which is known in statistics as independence—when the occurrence of one does not affect the occurrence of the other, in which case the probability that both A and B occur together is given by the product p_{A}p_{B}} p_{A}p_{B} of the probabilities. There is said to be a linkage disequilibrium between the two alleles whenever {p_{AB}} p_{AB} differs from p_{A}p_{B}} p_{A}p_{B} for any reaso: D=Pab - Pa*Pb. Then D'=D/D_max

    #white: LOD<2 and D'<1
    #bright red: LOD>=2 and D'=1
    #shades of pink/red: LOD>=2 and D'<1
        #violet. A case with high LOd, bu tlow D', so we have support for one of the metrics, and the implicated markers have bright red with other snps.
    #blue: LOD<2 and D'=1
    #https://www.broadinstitute.org/haploview/ld-display
    #https://en.wikipedia.org/wiki/Linkage_disequilibrium

#hay r2 mu bajas, pero no nos importa mucho porque no vamos a predeicr la presentica de un alelo con la frecuencia de otro, eso estará bien si quisesemos computar (rellenar datos sin genotipado). Mñas indo D' vs R2 en "https://www.biostars.org/p/133204/"

######RESULTADOS: 

#Add manually from haploview the results of blocks
#open dataframe
haplotypes_blocks = data.frame(selected_chr=NA, selected_gen=NA, selected_snp=NA, selected_block=NA)

#add a row for each snp
haplotypes_blocks = rbind.data.frame(haplotypes_blocks, 
    cbind.data.frame(selected_chr="chr_20", selected_gen="ptpn1", selected_snp="rs6067472", selected_block="block_1"),
    cbind.data.frame(selected_chr="chr_20", selected_gen="ptpn1", selected_snp="rs10485614", selected_block="block_1"),
    cbind.data.frame(selected_chr="chr_20", selected_gen="ptpn1", selected_snp="rs2143511", selected_block="block_1"),
    cbind.data.frame(selected_chr="chr_20", selected_gen="ptpn1", selected_snp="rs6020608", selected_block="block_1"),
    cbind.data.frame(selected_chr="chr_20", selected_gen="ptpn1", selected_snp="rs968701", selected_block="block_1"))

#remove the first row with NAs
haplotypes_blocks = haplotypes_blocks[-1,]
haplotypes_blocks

#add variable for the combinations of chromosomes and haplotype blocks
haplotypes_blocks$selected_chr_block = interaction(haplotypes_blocks$selected_chr, haplotypes_blocks$selected_block)
#check
identical(as.character(haplotypes_blocks$selected_chr_block), paste(haplotypes_blocks$selected_chr, ".", haplotypes_blocks$selected_block, sep=""))

####### estimate frequencies of haplotypes 

#haplo.stats used: Detailed tutorial here ("https://cran.r-project.org/web/packages/haplo.stats/vignettes/manualHaploStats.pdf")
require(haplo.stats)

#extract the names of haplotype blocks with signifcan snps (under FDR<0.1)
vector_haplotype_blocks = unique(haplotypes_blocks[which(haplotypes_blocks$selected_snp %in% unique(geno_pheno_results_fdr_0.1$snp_to_test)),]$selected_chr_block)

#for each block
hap_prob_block_insertion =list()
hap_prob_block_no_trim =list()
for(i in 1:length(vector_haplotype_blocks)){

    #select the [i] block
    selected_chr_block = vector_haplotype_blocks[i]

    #extract snps of the [i] block
    tag.SNPs_block = haplotypes_blocks[which(haplotypes_blocks$selected_chr_block == selected_chr_block),]$selected_snp

    ##make a matrix of genotypes of both snps using SNPassoc
    geno_block = make.geno(data=myData_ptpn1, SNPs.sel=tag.SNPs_block) #Matrix of alleles, such that each locus has a pair of adjacent columns of alleles, and the ORDER of columns corresponds to the order of loci on a chromosome. 
    ncol(geno_block) == 2*length(tag.SNPs_block) #If there are K loci, then ncol(geno) = 2*K. Rows represent alleles for each subject.
    summaryGeno(geno_block)
    str(geno_block)

    #set the seed
    set.seed(17865)

    ## Haplotype estimation con algoritmo EM de inserción progresiva
    hap_prob_block_insertion[[i]] = haplo.em(geno_block, locus.label=tag.SNPs_block) #Para marcadores genéticos medidos en sujetos no emparentados, con linkage desconocido, haplo.em calcula estimas de maxima verosimilud (maximum likelihood) de las probabilidades de los haplotipos. A causa de que podría haber más de un par de haplotipos consistentes con lo observado, también se computan las probabilidades posteirores de los pares de haplotipos para cada individuo. Las estimaciones típicas de haplotype frequencies intentan enumerar todos los pares posibles de haplotipos al comienzo, sin embargo haplo.em da usa un algoritmo de "inserción progresiva": i) El cual va progresivamente inserta un grupo de loci (varios locus) en el haplotipo; ii) Corre una iteración, en esa iteración quita pares de haplotipos por individuo cuando la probabilidad posterior del par está por debajo de un threshold especificado previamente; iii) Vuelve a meter un nuevo grupo de loci. Este proceso se repite hasta que se han metido todos los loci en el haplotipo. El usuario puede elgir el tamaño del grupo de loci que se inserta en cada iteración (batch size). Si el batch size es igual al número total de loci, y el threshold para ir eliminando pares de haplotipos es 0, entonces el algoritmo de "inserción progresiva" se reduce a un EM típico. Este algortimo está basado en "snphap".
    print.haplo.em(hap_prob_block_insertion[[i]]) #The print methods shows the haplotypes and their estimated frequencies, followed by the final log-likelihood statistic and the lr stat for no LD, which is the likelihood ratio test statistic contrasting the lnlike for the estimated haplotype frequencies versus the lnlike under the null assuming that alleles from all loci are in linkage equilibrium. We note that the trimming by the progressive insertion algorithm can invalidate the lr stat and the degrees of freedom (df).
    summary(hap_prob_block_insertion[[i]], show.haplo=TRUE, nlines=7) #The summary method for a haplo.em object on save.em shows the list of haplotypes per subject, and their posterior probabilities. Por ejemplo: Si tienes dos snps, te pone dos columnas con el snp1 y snp2 para el haplotipo1, luego lo repite para el haplotipo 2, por tanto otras dos columnas. Si hubiese un tercer SNP habría una tercera columna. The first part of the summary output lists the subject id (row number of input geno matrix), the codes for the haplotypes of each pair, and the posterior probabilities of the haplotype pairs. The second part gives a table of the maximum number of pairs of haplotypes per subject, versus the number of pairs used in the final posterior probabilities. The haplotype codes remove the clutter of illustrating all the alleles of the haplotypes, but may not be as informative as the actual haplotypes themselves. To see the actual haplotypes, use the show.haplo=TRUE option, as in the following example.

    ## Haplotype estimation con algoritmo EM sin inserción progresiva
    #El proceso iterativo de eliminar pares de haplotimos cuya probabilidad posterior no pasa un límite dado puede ahcer que el likelihood ratio test y los grados de libertad no sean válidos. Por esta razón, repetimos la estimación del haplotipo pero con el algoritmo sin inserción progesiva, para ver si el LRtest sale igual. Esto se hace indicando que el número de loci que se inserta sea directamente el total de ellos desde la primera iteración, por tanto acaba ahí. Además hay que establacer que no hay limite de probabilidad posterior a partir del cual eliminar pares de haplotipos por sujeto. 
    hap_prob_block_no_trim[[i]] = haplo.em(geno_block, locus.label=tag.SNPs_block, control = haplo.em.control(insert.batch.size = length(tag.SNPs_block), min.posterior = 0))
    print.haplo.em(hap_prob_block_no_trim[[i]]) #Una vez hechos los bloques de SNPs he testado si las frecuencias "reales" del haplotipo se desvían de las frecuencias esperadas si los polimorfismos no estuvieran en linkage, es decir, no formasen un haplotipo. En todos los casos las frecuencias de las diferentes combinaciones se desvían significativamente de lo esperado en equilibrio (es decir, si no hubiera linkage ó haplotipo; P = 0)
    summary(hap_prob_block_no_trim[[i]], show.haplo=TRUE, nlines=7)
}
names(hap_prob_block_insertion) <- vector_haplotype_blocks
names(hap_prob_block_no_trim) <- vector_haplotype_blocks

#print the relevant results per block (only haplotpy pahse witohut trim)
for(i in 1:length(vector_haplotype_blocks)){

    #select the results of [i] block
    selected_chr_block = vector_haplotype_blocks[i]

    #print block name
    print("")
    print(paste("#############################"))
    print(paste("#########", selected_chr_block, "#############", sep=""))
    print(paste("#############################"))
    print("")

    #print the results of LD witout trimming
    print("")     
    print(paste("######### without trim #######"))
    print("")           
    print.haplo.em(hap_prob_block_no_trim[[selected_chr_block]])
    print("") 
}

#RESULTS: ALL BLOCKS DIFFERD FROM EXPECTED UNDER NO LD

#####################################################
##### haplotype associations ########
#####################################################

#Vamos a testar asociación entre fenotipos y haplotipos que contienen SNPs que de forma individual se asocian con esos fenotipos. Para cada haplotipo, se miran todas las posibles combinaciones de alelos, y para cada una de ellas se hace un conteo: Por ejemplo, cuantos individuos tienen ATAA en alguno de los cromosomas (!! tiene que estar esa combinación de alelos en el mismo cromosoma). Si se analiza con herencia dominante, se le da un valor de 1 a aquellos individuos que tenga una ó dos copias del haplotipo, es decir, que tenga ese haplotipo en uno ó los dos cromosomas. Al resto de individuos se le da cero. De esa forma tenemos una varaible para ese haplotipo, para la cual se verá si hay asociación con el fenotipo de interés. Lo mismo se hace con las otras posibles combinaciones de aleleos de esos snps o haplpitpos. Cada una de esas covariables se meteran en los modelos con el fenotipo como respuesta. En el caso de aditivo, individuos con una copia tendrán 1, y con 2 tendrán un valor de dos. Si se analiza bajo modelo recesivo, solo se le da un valor de 1 a aquellos individuos que tengan dos copias del haplotipo. Ahí está el problema del N, si solo suman aquellos individuos con las dos copias del haplotipo, solo haplotipos con una frecuencia alta van a tener un N suficiente. En nuestro caso solemos tener uno ó dos haplotipos muy frecuentes, mientras que el resto son poco comunes, para estos últimos no va a haber N bajo un modelo recesivo, así que pasamos de usarlo. 
    #RESUMEN:cada haplotipo es una variable, con 0,1,2, dependiendo del modelod eherencia, si es odminante, con que haya una copia del haplotipo en un indivudup ya es 1, si no hay ninguna copia es 0. En el caso de los recesivos tiene que haber dos copias para que haya 1. En aditivo, con dos copias es 2. De esta forma se codifica una variable por cada haplotipo, con 0,1,2, y se testa asociación con el fenotipo. En el caso del score se hace por separado con cada covaraible de haplotipo y luego con los scores especificos se obtiene un score global. 

#haplo.stats used: Detailed tutorial here ("https://cran.r-project.org/web/packages/haplo.stats/vignettes/manualHaploStats.pdf")
require(haplo.stats)

#snps ordered by position
ptpn1_position[order(ptpn1_position$pos),]

#create a list for saving the results
haplo_list_per_block = list()

#for each haplotype block
for(i in 1:length(vector_haplotype_blocks)){

    #selected the [i] block
    selected_chr_block = vector_haplotype_blocks[i]

    #selected snps inside the [i] block
    selected_snps = haplotypes_blocks[which(haplotypes_blocks$selected_chr_block==selected_chr_block),]$selected_snp
    
    ##make a matrix of genotypes of both snps using SNPassoc
    geno_block = make.geno(data=myData_ptpn1, SNPs.sel=selected_snps) #Matrix of alleles, such that each locus has a pair of adjacent columns of alleles, and the ORDER of columns corresponds to the order of loci on a chromosome. 
    ncol(geno_block) == 2*length(selected_snps) #If there are K loci, then ncol(geno) = 2*K. Rows represent alleles for each subject.
    summaryGeno(geno_block)
    #str(geno_block)
    
    #phenotypes associated with these snps according to FDR<0.1
    significant_phenotypes = unique(geno_pheno_results_fdr_0.1[which(geno_pheno_results_fdr_0.1$snp_to_test %in% selected_snps),]$selected_pheno)
  
    #create a new list for nesting it inside haplo_list_per_block
    haplo_list_per_pheno = list()

    #for each phenotype
    for(p in 1:length(significant_phenotypes)){

        #select the [p] phenotype
        selected_phenotype = significant_phenotypes[p]

        #covariables
        if(!selected_phenotype %in% c("CVi_BP", "SBP", "DBP", "TG", "TC", "LDL", "HDL", "LDL_HDL", "Apo_A1", "Apo_B", "ApoB_ApoA1", "apoB_LDL", "TG_HDL", "Insulin", "Leptin_ng_ml", "HOMA", "QUICKI",  "TC_HDL",  "risk_score_for_log")){

            #only sex, age and center for adiposity variables
            covariables_data = cbind(myData_ptpn1$CRF_sex, myData_ptpn1$CRF_age, myData_ptpn1$center)
            covariables_names = "+CRF_sex+CRF_age+center"# we use cbind instead of cbind.data.frame becauase we want the facotrs covariables transformed into numbers. In this way it is done in the example of "haplo.score": "For a binary trait, adjusted for sex and age: x <- cbind(male, age)". Inf not, we have an error in that function. It should not be a problem because myData_ptpn1$CRF_sex  and myData_ptpn1$center are considered as factors in R, so when you bind them with myData_ptpn1$CRF_age, they become into a numberic variables with as many numbers as categories. For example, the first column (CRF_sex) has two categories and the third column has 
                #length(unique(covariables_data[,1])) == length(unique(myData_ptpn1$CRF_sex))
                #length(unique(covariables_data[,3])) == length(unique(myData_ptpn1$center))
            #Similarly, rows of the first column in covariables_data that are male in myData_ptpn1 have only 1. In contrast, rows of the first column in covariables_data that are female in myData_ptpn1 have only 2. Everything is right.
                unique(covariables_data[which(myData_ptpn1$CRF_sex == "male"),1])
                unique(covariables_data[which(myData_ptpn1$CRF_sex == "female"),1])
        } else {

            #BMI in case the response is a cardiovascular factor
            covariables_data = cbind(myData_ptpn1$CRF_BMI, myData_ptpn1$CRF_sex, myData_ptpn1$CRF_age, myData_ptpn1$center)#the same explained about about myData_ptpn1$CRF_sex and myData_ptpn1$center applied here.
            covariables_names = "+CRF_BMI+CRF_sex+CRF_age+center"            
        }

        #for discrete phenotypes (obesity) 
        if(selected_phenotype %in% c("obesity", "CVi_BP")){
            family_error = "binomial"
            response = paste("as.numeric(levels(myData_ptpn1$", selected_phenotype, ")[myData_ptpn1$", selected_phenotype, "])",sep="")#Transform the discrete phenotype to numeric. No problem with haplo.score modelling because we set the family error as "binomial". With this code we transform a factor into a numeric using exactly the same levels of that factor. You select the levels with "levels(myData_ptpn1$obesity)", and then we assign to each individual ("[myData_ptpn1$obesity]") the corresponding number (level) that will be considered as a number thanks to "as.numeric". See an exmaple with obesity to check that everything goes fine ("more information here"):
                #new_obesity_cont = eval(parse(text=response))
                #unique(new_obesity_cont[which(myData_ptpn1$obesity==0)])#all 0
                #unique(new_obesity_cont[which(myData_ptpn1$obesity==1)])#all 1
        } else{ #for continous phenotypes
            family_error = "gaussian"
            response = paste("log(myData_ptpn1$", selected_phenotype, ")", sep="") #apply log transformation      
        }

        #create a new list for nesting it inside haplo_list_per_block
        haplo_list_per_haplo_score = list()
        haplo_list_per_regression = list()

        #heritage models to test
        heritage_models_test = c("dominant", "additive")

        #for each heritage model
        for(m in 1:length(heritage_models_test)){

            #selected model
            selected_model = heritage_models_test[m]

            #repeat analyses to check convergence
            for(r in 1:3){
                
                #print the combination of chr, block, phenotype and model printed
                print(" ")
                print("#######################")                
                print(toupper(paste(selected_chr_block, "-", selected_phenotype, "-", selected_model, "-repetion", r, sep="")))
                print("#######################")   
                print(" ")


                ## obtain p.val for the association of the haplotype (haplotype score)
                #This method also has an option to compute permutation p-values, which may be needed for sparse data when distribution assumptions may not be met

                #run the model for the [r] repetition
                repetition = haplo.score(eval(parse(text=response)), geno_block, trait.type=family_error, x.adj = covariables_data, locus.label=selected_snps, haplo.effect=selected_model, min.count=10, simulate=TRUE, sim.control=score.sim.control(min.sim=10000, p.threshold=0.01))

                #print the results
                print(repetition)
            
                #haplo.effect: haplo.score allows non-additive effects for scoring haplotypes. The possible effects for haplotypes are additive, dominant, and recessive. Under recessive effects, fewer haplotypes may be scored, because subjects are required to be homozygous for haplotypes. Furthermore, there would have to be min.count such persons in the sample to have the recessive effect scored. Therefore, a recessive model should only be used on samples with common haplotypes.

                #min.count: We restrict the analysis to get scores for haplotypes with a minimum sample count using min.count=5. For the haplo.score, the skip.haplo and min.count parameters control which rare haplotypes are pooled into a common group. The min.count parameter is a recent addition to haplo.score, yet it does the same task as skip.haplo and is the same idea as haplo.min.count used in haplo.glm.control for haplo.glm. As a guideline, you may wish to set min.count to calculate scores for haplotypes with expected haplotype counts of 5 or greater in the sample. We concentrate on this expected count because it adjusts to the size of the input data. Este es mejor que skip.haplo porque se ajusta al tamaño de los datos usados. En nuestro caso, si usamos skip.haplo=0.05 (como usan los de SNPassoc en su manual), harían falta 105 counts (105 casos con un haplotipo dado) para que de no se considerare poco común, ya que 105/(2*1045)=0.05 (nº counts del haplotipo partido número máximo de counts ya que cada individuo tiene 2 cromosomas). Si por el contrario, tuviesemos solo 200 indivduos, con 20 counts sería suficiente (20/(2*200)=0.05). La frequencia no tiene en cuenta el número de individuos que tienes. Usamos un poco más del número mínimo de counts sugerido por el autor del paquete, 10 en vez de 5, para curarnos en salud. For more explanation on handling rare haplotypes, see section 5.6. 

                #simulate=TRUE: Simulated or permutation p.values: To estimate the sampling distribution of the test statistic we need many samples generated under the strong null hypothesis. If the null hypothesis is true, changing the exposure would have no effect on the outcome. By randomly shuffling the exposures we can make up as many data sets as we like. If the null hypothesis is true the shuffled data sets should look like the real data, otherwise they should look different from the real data. The ranking of the real test statistic among the shuffled test statistics gives a p-value ("http://faculty.washington.edu/kenrice/sisg/SISG-08-06.pdf"; diapo 5-6). This is exactly donde by simulate=TRUE. When simulate=TRUE, haplo.score gives simulated p-values. Simulated haplotype score statistics are the re-calculated score statistics from a permuted re-ordering of the trait and covariates and the original ordering of the genotype matrix. The simulated p-value for the global score statistic (Global sim. p-val) is the number of times the simulated global score statistic exceeds the observed, divided by the total number of simulations. Creo que cuanto mayor es el estadístico, más explicativo es, por tanto cuanto menor sea el número de simulaciones en las que el estadístico simulado sea mayor que el real, mejor será (p pequeños). Likewise, simulated p-value for the maximum score statistic (Max-stat sim. p-val) is the number of times the simulated maximum haplotype score statistic exceeds the observed maximum score statistic, divided by the total number of simulations. The maximum score statistic is the maximum of the square of the haplotypespecific score statistics, which has an unknown distribution, so its significance can only be given by the simulated p-value. Intuitively, if only one or two haplotypes are associated with the trait, the maximum score statistic should have greater power to detect association than the global statistic. The score.sim.control function manages control parameters for simulations. The haplo.score function employs the simulation p-value precision criteria of Besag and Clifford[6]. These criteria ensure that the simulated p-values for both the global and the maximum score statistics are precise for small p-values. The algorithm performs a user-defined minimum number of permutations (min.sim) to guarantee sufficient precision for the simulated p-values for score statistics of individual haplotypes. Permutations beyond this minimum are then conducted until the sample standard errors for simulated p-values for both the globalstat and max-stat score statistics are less than a threshold (p.threshold * p-value). The default value for p.threshold= 1/4 provides a two-sided 95% confidence interval for the p-value with a width that is approximately as wide as the p-value itself. Effectively, simulations are more precise for smaller p-values. POR DEFECTO NÚMERO MINIMO DE SIMULACIONES ES 1000, pero se puede cambiar.

                #score.sim.control: In simulations for haplo.score, employ the simulation p-value precision criteria of Besag and Clifford (1991).  The criteria ensures both the global and the maximum score statistic simulated p-values be precise for SMALL p-values.  First, perform min.sim simulations to guarantee sufficient precision for the score statistics on individual haplotypes.  Then continue simulations as needed until simulated p-values for both the global and max score statistics meet precision requirements set by p.threshold. The p.threshold is a paremeter used to determine p-value precision from Besag and Clifford (1991).  For a p-value calculated after min.sim simulations, continue doing simulations until the p-value's sample standard error is less than p.threshold * p-value.  The dafault value for p.threshold = 1/4 corresponds approximately to having a two-sided 95% confidence interval for the p-value with a width as wide as the p-value itself. Therefore, simulations are more precise for smaller p-values. Additionally, since simulations are stopped as soon as this criteria is met, p-values may be biased high. RESUMEN CON EJEMPLO: Empieza a correr simulaciones hasta el número mínimo, a partir de ahí para cuando el SE de los p.values sea menor que p.value*threshold. Es decir, cogemos el p.value global de cada simulación y con todos hasta la última se calcula el SE, parando cuando el SE es menor que el p.value de la simulación actual multiplicado por el threshold. Por ejemplo: Llego a 2000 simulaciones que es el mínimo, el último p.value es 0.095, y he pueso un p.threshold de 0.25 (i.e. 1/4), entonces el limite sería 0.095*0.25=0.02, cuando el SE de los p.values globales a lo largo de todas las simulaciones sea menor que 0.02, se para. Lógicamente, si el p.value es muy bajo, el SE de los p.values tendrá que ser mucho menor, porque 0.05*0.25=0.01, no se parará hasta que la SE de los p.values sea 0.01. Esto tiene sentido, cuando más bajo sea un p.value más exacto tiene que ser, cuando sea más alto, es decir, cuanto más alejado esté de 0.05 menos importante es que sea poco exacto. Así se ahorra tiempo de computación. 
                #En nuestro caso vamos a poner 2000 como número mínimo de simulaciones. Y 0.01 (1/100) como p.threshold, para que los p.values sea más fiable aun siendo grandes, pero realmente no haría falta, solo en aquellos casos rondando el 0.05. 

                #PARA TESTAR ASOCIACION con el rasgo de subgrupos dentro del haplotipo (una selección de solo unos SNPs) se puede usar haplo.score.slide ("https://cran.r-project.org/web/packages/haplo.stats/vignettes/manualHaploStats.pdf"). 

                #PARA CALCULAR EL R2 DE LO QUE EXPLICA EL HAPLOTIPO: haplo.power.qt ("https://cran.r-project.org/web/packages/haplo.stats/vignettes/manualHaploStats.pdf")

                #Si se diese el caso de que el p.value global sea signifciativo pero no haya efecto de los haplotipos, cambia el valor de epsilon en control: Mira "Score Statistic Dependencies: the eps.svd parameter" en el manual.
            }

            #fit haplo.score and save it
            set.seed(3456)
            haplo_list_per_haplo_score[[m]] = haplo.score(eval(parse(text=response)), geno_block, trait.type=family_error, x.adj = covariables_data, locus.label=selected_snps, haplo.effect=selected_model, min.count=10, simulate=TRUE, sim.control=score.sim.control(min.sim=10000, p.threshold=0.01))#details about the function inside the latter loop

            #run regression
            haplo_list_per_regression[[m]] = haplo.glm(paste(response, "~geno_block", covariables_names, sep=""), data=myData_ptpn1, family=family_error, locus.label=selected_snps, na.action="na.geno.keep", control = haplo.glm.control(haplo.effect=selected_model, haplo.min.count=10, glm.c=glm.control(maxit=2000)))# we increased the number of iterations to get convergence

            #haplo.effect: Usamos el modelo de dominancia el que sale significativo: 

            #na.action: This is applied to the model.frame. The default value of na.action=na.geno.keep will keep observations with some (but not all) missing alleles, but exclude observations missing any other data (e.g., response variable, other covariates, weight). The EM algorithm for ambiguous haplotypes accounts for missing alleles. Similar to the usual glm, na.fail creates an error if any missing values are found, and a third possible alternative is na.exclude, which deletes observations that contain one or more missing values for any data, including alleles.

            #The issue of deciding which haplotypes to use for association is critical in haplo.glm. By default it will model a rare haplotype effect so that the effects of other haplotypes are in reference to the baseline effect of the one common happlotype. The rules for choosing haplotypes to be modeled in haplo.glm are similar to the rules in haplo.score: by a minimum frequency or a minimum expected count in the sample. Two control parameters in haplo.glm.control may be used to control this setting: haplo.freq.min may be set to a selected minimum haplotype frequency, and haplo.min.count may be set to select the cut-off for minimum expected haplotype count in the sample. The default minimum frequency cut-off in haplo.glm is set to 0.01. More discussion on rare haplotypes takes place in section 6.7.4.

            #haplo.freq.min indica a partir de frecuencia se considera un haplotipo poco frecuente y por tanto se mete en la categoria de raro. haplo.min.info indica a patir de que frecuencia deja de tenerse en cuenta un haplotipo para calcualr los SE. Por tanto, si haplo.freq.min = haplo.min.info, quiere decir que todos los haplotipos includios en raros quedan fuera de las estimaciones del SE para esa categoria (rare) y por tanto saldrá NA. Normalmente, lo que pasará es que haya un pequeño grupo de haploitpos entre haplo.freq.min  y haplo.min.info, por tanto, las estimacioes de la categoria rare pueden ser muy problemáticas, así que tener muuucho cuidado. USAMOS MEJOR haplo.min.count, para homogeneizar con haplo.score

            #haplo.min.count : The minimum number of expected counts for a haplotype from the sample to be included in the model.  The count is based on estimated haplotype frequencies.  Suggested minimum is 5. haplo.min.count=5 sería una frequencia mínima de 5/(2*1045)=0.0024. Este es mejor que haplo.freq.min porque se ajusta al tamaño de los datos usados. En nuestro caso, si usamos haplo.freq.min=0.05 (como usan los de SNPassoc en su manual), harían falta 105 counts (105 casos con un haplotipo dado) para que de no se considerare poco común, ya que 105/(2*1045)=0.05 (nº counts del haplotipo partido número máximo de counts ya que cada individuo tiene 2 cromosomas). Si por el contrario, tuviesemos solo 20 indivduos, con 20 counts sería suficiente (20/(2*200)=0.05). La frequencia no tiene en cuenta el número de individuos que tienes. Usamos un poco más del número mínimo de counts sugerido por el autor del paquete, 10 en vez de 5, para curarnos en salud.  

            #The haplotype chosen for the baseline in the model is the one with the highest frequency. Sometimes the most frequent haplotype may be an at-risk haplotype, and so the measure of its effect is desired. To specify a more appropriate haplotype as the baseline in the binomial example, choose from the list of other common haplotypes, fit.bin$haplo.common. To specify an alternative baseline, such as haplotype 77, use the control parameter haplo.base and haplotype code. See section 6.7.2 Selecting the Baseline Haplotype of the manual for further information. 

            #reset the set seed ("https://stackoverflow.com/questions/40399937/how-to-reset-a-seed-value")
            rm(.Random.seed, envir=globalenv())#in this way we avoid to get the same results during repetitions to check convergence for the next models/phenotypes
        } 

        #indicate the name of the [m] model
        names(haplo_list_per_haplo_score) <- heritage_models_test 
        names(haplo_list_per_regression) <- heritage_models_test 

        #bind results from haplo score and regression into one single list
        haplo_list_per_haplo_score_regression = list(haplo_list_per_haplo_score, haplo_list_per_regression)

        #add names (first haplo score, second regression)
        names(haplo_list_per_haplo_score_regression) <- c("haplo_score", "regression")

        #check that the order of the heritage models is correct for regression (in the case of haplo.score we can see it in the output)
        eval(parse(text=paste("haplo_list_per_haplo_score_regression$regression$additive$control$haplo.effect", sep=""))) == "additive" 
        eval(parse(text=paste("haplo_list_per_haplo_score_regression$regression$dominant$control$haplo.effect", sep=""))) == "dominant"

        #save these results into the list of pheno
        haplo_list_per_pheno[[p]] = haplo_list_per_haplo_score_regression
    }

    #set names of the list
    names(haplo_list_per_pheno) <- significant_phenotypes

    #save results in the block list
    haplo_list_per_block[[i]] <- haplo_list_per_pheno  
}
  
#set names of the final list
names(haplo_list_per_block) <- vector_haplotype_blocks


###save significant results

#load allele names to interpretation
allele_names = read.csv("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/data/snps/alleles_ptpn1_v2.csv", header=TRUE)
    #IMPORTANT NOTE: I have matched the allele names of HELENA and ncbi (at 10/09/2019). Now, they all matched, but I have noted that UCP alleles have changed in ncbi. For example, alleles that I changed to match ncbi, now are in ncbi exactly as original HELENA. If you check this for these SNPs, and you see changes no panic! The important thing is you are using the complementary chain and the first allele is always the major. Indeed, in many cases both options (i.e., both chains) are included as synonimous in ncbi (i.e., HGVS).

#check that there is a correspondence between ncbi and helena alleles (first major, second minor)
check_helena_ncbi_names=NULL
for(i in 1:nrow(allele_names)){

    #select the [i] row
    selected_row = allele_names[i,]

    #extract helena and ncbi names
    helena_name = as.vector(selected_row$helena)
    ncbi_name = as.vector(selected_row$ncbi)

    #if both names are equal
    if(helena_name == ncbi_name){

        #correct
        check_helena_ncbi_names = append(check_helena_ncbi_names, "CORRECT")

    } else {#if they are different transform the helena name to the corresponding alleles of ncbi to check that they have the adequate order (major is the first and minor the second)

        #extract major and minor alleles
        helena_major = strsplit(helena_name, split="/")[[1]][1]
        helena_minor = strsplit(helena_name, split="/")[[1]][2]

        #copy helena name
        new_helena_major =  helena_major
        new_helena_minor =  helena_minor

        #convert major
        #If T is included in helena name, convert T into A
        if(grepl("T", helena_major)){
            new_helena_major = gsub("T", "A", new_helena_major)
        } 
        #If C is included in helena name, convert C into G
        if(grepl("C", helena_major)){
            new_helena_major = gsub("C", "G", new_helena_major)                
        }
        #If G is included in helena name, convert G into C        
        if(grepl("G", helena_major)){
            new_helena_major = gsub("G", "C", new_helena_major)                
        }
        #If A is included in helena name, convert A into T               
        if(grepl("A", helena_major)){
            new_helena_major = gsub("A", "T", new_helena_major)                
        }  

        #convert minor
        #If T is included in helena name, convert T into A
        if(grepl("T", helena_minor)){
            new_helena_minor = gsub("T", "A", new_helena_minor)
        } 
        #If C is included in helena name, convert C into G
        if(grepl("C", helena_minor)){
            new_helena_minor = gsub("C", "G", new_helena_minor)                
        }
        #If G is included in helena name, convert G into C        
        if(grepl("G", helena_minor)){
            new_helena_minor = gsub("G", "C", new_helena_minor)                
        }
        #If A is included in helena name, convert A into T               
        if(grepl("A", helena_minor)){
            new_helena_minor = gsub("A", "T", new_helena_minor)                
        }  

        #bind converted major and minor alleles
        new_helena_name = paste(new_helena_major, new_helena_minor, sep="/")

        #if the new name is equal to ncbi name, perfect, but if not, we have a problem
        if(ncbi_name == new_helena_name){
            #correct
            check_helena_ncbi_names = append(check_helena_ncbi_names, "CORRECT")
        } else {

            #incorrect
            check_helena_ncbi_names = append(check_helena_ncbi_names, "INCORRECT")            
        }
    }
}
#there is no incorrect case
!"INCORRECT" %in% check_helena_ncbi_names

#reorder with ucp_pos (same positon in chormosome)
allele_names = allele_names[match(ptpn1_position$snp, allele_names$snp),]
allele_names

##create a variable for allele interpretation
#change colnames of allele_names to merge with haplotypes_blocks and have alleles and blocks for interpretation
colnames(allele_names)[which(colnames(allele_names) == "snp")] <- "selected_snp"
#select thos snps from allele_names that are included in any haplotype block
allele_names_in_haplos = allele_names[which(allele_names$selected_snp %in% haplotypes_blocks$selected_snp),]
#ensure that haplotypes_blocks and allele_names have snps in the same order to bind them
allele_names_in_haplos = allele_names_in_haplos[match(haplotypes_blocks$selected_snp, allele_names_in_haplos$selected_snp),]
#bind them
allele_interpretation = cbind.data.frame(haplotypes_blocks, allele_names_in_haplos)
#check
summary(allele_interpretation$selected_snp == haplotypes_blocks$selected_snp)
#select rows that you are interested
allele_interpretation = allele_interpretation[,c(1,2,3,4,5,8,9)]
allele_interpretation

##create a variable with the same notation for allele names than in myData_ptpn1 (1 and 2)
#in the helena database the first number (1) is for the first letter in alphabetic order. Therefore A will be always 1
#create a new variable
allele_interpretation$helena_number <- NA
#for each snp
for(i in 1:nrow(allele_interpretation)){
    
    #select the [i] row
    selected_row = allele_interpretation[i,]
    
    #select the helena allele names
    helena_name = selected_row$helena

    #copy the helena name
    helena_name_new = helena_name

    #extract the allele names according to helena
    helena_alleles = strsplit(as.vector(helena_name), split="/")[[1]]

    #always A is the first, because number are afabetically ordered in genotype levels in Helena
    if("A" %in% helena_alleles){
        helena_name_new = gsub("A", "1", helena_name_new)
        helena_name_new = gsub(helena_alleles[!"A" == helena_alleles], "2", helena_name_new)                                  
    } else { #if not, C first
        if("C" %in% helena_alleles){ 
            helena_name_new = gsub("C", "1", helena_name_new)
            helena_name_new = gsub(helena_alleles[!"C" == helena_alleles], "2",  helena_name_new)                                                  
        } else { #if not G first (T always will be the last)
            if("G" %in% helena_alleles){
                helena_name_new = gsub("G", "1", helena_name_new)
                helena_name_new = gsub(helena_alleles[!"G" == helena_alleles], "2", helena_name_new)               
            }
        }            
    }

    allele_interpretation[i,]$helena_number <- helena_name_new
}
allele_interpretation

#check the conversion from allele name to number
allele_check = NULL
for(i in 1:nrow(allele_interpretation)){
    
    #select the [i] row
    selected_row = allele_interpretation[i,]
    
    #select the helena allele names
    helena_name = selected_row$helena

    #extract the allele names according to helena
    helena_alleles = strsplit(as.vector(helena_name), split="/")[[1]]

    #extract the helena number
    helena_number = selected_row$helena_number

    #always A is the first, because number are afabetically ordered in genotype levels in Helena
    if("A" %in% helena_alleles){
        helena_number = gsub("1", "A", helena_number)
        helena_number = gsub("2", helena_alleles[!"A" == helena_alleles], helena_number)                                  
    } else { #if not, C first
        if("C" %in% helena_alleles){ 
            helena_number = gsub("1", "C", helena_number)
            helena_number = gsub("2", helena_alleles[!"C" == helena_alleles], helena_number)                                                  
        } else { #if not G first (T always will be the last)
            if("G" %in% helena_alleles){
                helena_number = gsub("1", "G", helena_number)
                helena_number = gsub("2", helena_alleles[!"G" == helena_alleles], helena_number)               
            }
        }            
    }

    #check that the change was properly
    allele_check = append(allele_check, selected_row$helena==helena_number)
}
#check all the check across snps
summary(allele_check)

#for each block
singificant_results_per_blocks = list() #empty list for saving results per block
for(i in 1:length(vector_haplotype_blocks)){

    #selected the [i] block
    selected_chr_block = vector_haplotype_blocks[i]

    #selected snps inside the [i] block
    selected_snps = haplotypes_blocks[which(haplotypes_blocks$selected_chr_block==selected_chr_block),]$selected_snp

    #phenotypes associated with these snps according to FDR<0.1
    significant_phenotypes = unique(geno_pheno_results_fdr_0.1[which(geno_pheno_results_fdr_0.1$snp_to_test %in% selected_snps),]$selected_pheno)

    #for each significan pheno
    singificant_results_per_pheno = list() #empty list for saving results per pheno
    for(p in 1:length(significant_phenotypes)){
        
        #select the [p] phenotype
        selected_phenotype = significant_phenotypes[p]
        
        #heritage models to test
        heritage_models_test = c("dominant", "additive")
        
        #for each heritage model
        singificant_results_per_model = list() #empty list for saving results per model
        for(m in 1:length(heritage_models_test)){

            #select the [m] model
            selected_model = heritage_models_test[m]

            #extract score results for the [i] block, [p] pheno and [m] model. Score is the first element always ([[1]])
            result_score = eval(parse(text=paste("haplo_list_per_block$", selected_chr_block, "$", selected_phenotype, "$haplo_score$", selected_model,sep="")))

            #extract the corrected pval
            pval_score = as.vector(result_score$score.global.p.sim)

            #if the corrected pval is significant
            results_score_interval = list()#empty list for saving score pval and intervals
            if(pval_score<0.05){

                #save the pval as the second elment in the list
                results_score_interval[[1]] <- pval_score

                #load the regression for the [i] block, [p] pheno and [m] model. Regression is the second element always ([[2]])
                result_regresssion = eval(parse(text=paste("haplo_list_per_block$", selected_chr_block, "$", selected_phenotype, "$regression$", selected_model,sep="")))
                #result_regresssion
                #summary(result_regresssion, show.all.haplo=TRUE)

                #calculate confidence intervals and save them as the second elment in the list
                results_score_interval[[2]] = intervals(result_regresssion) #The rare haplotypes, those with expected counts less than haplo.min.count=5 (equivalent to having frequencies less than haplo.freq.min = 0.0113636363636364) in the above example), are pooled into a single category labeled geno.glm.rare. The haplotype chosen as the baseline category for the design matrix (most frequent haplotype, is the default) is labeled as haplo.base; more information on the baseline may be found in section 6.7.2.

                #Although we see the levels of the correction factors and the corresponding estimates compare to the reference level, the confidence intervals of haplotypes do not change if we change the reference levels of these factors because we are not including interaction between haplotype and the correction variable, like in regular GLMs. Therefore, the relationship between haplotypes and phenotype is the same across levels of the factors used to correct. I have checked this in some cases and it is true.

                #Para ver si los p.values de cada nivel son significagivos se aplica FDR. No hacemos Tukey porque requiere más de 6 niveles y que estén balanceados (no es el caso con tantas diferencias de count entre haplotipos), sino es overconservaitve. Además tendría que calcularlo a mano. https://rpubs.com/Joaquin_AR/236898. Tampoco me vale BF ni holms-BF porque asumen independencia de los test, en las comparaciones multiples de los niveles de un factor dentro de un modelo no hay indepedencia, si sin nivel es muy diferente, es más probable que también sea muy diferente del resto (dependencia positiva), por tanto BF y BF-holms serían overconservaitve, asumen más tests de los que hay. Por eso usamos FDR, que tiene en cuenta todo esto. Eso sí, con un nivel estricto del 0.05.  

                #extract again intervals and save them into an object to extract and correct pvals
                extract_intervals = results_score_interval[[2]]

                #extract variables names in intervals
                names_interval = row.names(extract_intervals[1:nrow(extract_intervals),])

                #select the last row to include in pval correction
                #If geno.block.rare exist
                if("geno_block.rare" %in% names_interval){

                    #we select the row before this. We don't use the pval of rare haplotypes. #SOLO los p.values de los niveles de haplotipo, nuestra familia de test incluye todos los niveles del factor de interés, es como cuando haces comparaciones mltiples entre los niveles de especie. Por eso no incluímos los tests de las covariables... TAMPOCO se incluye la p del nivel raro: Some datasets there may be only a few haplotypes between haplo.min.info and haplo.freq.min, and may yield misleading results for the rare haplotype coefficient. For this reason, we recommend that any inference made on the rare haplotypes be made with caution, if at all (6.7.4 Rare Haplotypes and haplo.min.info del tutorial).
                    last_row = abs(1-which(names(extract_intervals[,1]) == "geno_block.rare"))
                } else{

                    #if not and the response is adiposity, the first row after haplos is CRF_sexmale
                    if(!selected_phenotype %in% c("CVi_BP", "SBP", "DBP", "TG", "TC", "LDL", "HDL", "LDL_HDL", "Apo_A1", "Apo_B", "ApoB_ApoA1", "apoB_LDL", "TG_HDL", "Insulin", "Leptin_ng_ml", "HOMA", "QUICKI",  "TC_HDL",  "risk_score_for_log")){
                        last_row = abs(1-which(names(extract_intervals[,1]) == "CRF_sexfemale"))
                    } else{

                        #if the response is ECV, then we have as first predictor CRF_BMI
                        last_row = abs(1-which(names(extract_intervals[,1]) == "CRF_BMI"))                        
                    }
                }

                #extract the interesting pvals
                pvals = extract_intervals[1:last_row, "P-val"]

                #remove intercept pval
                pvals = pvals[-which(names(pvals) == "(Intercept)")]

                #calculate and save corrected pvals as the third object of the list
                results_score_interval[[3]] = p.adjust(pvals, method="BH")

                #extract allele names of implicated snps
                results_score_interval[[4]] = allele_interpretation[which(allele_interpretation$selected_chr_block %in% selected_chr_block),]

                #set the names of the list with pvals and intervals (score always first, second interval, third corrected pvals of intervals)
                names(results_score_interval) <- c("score_pval", "intervals", "adjusted_pvals", "allele_interpretation")
            }

            #save pval score and intervals as an element in the list per model
            singificant_results_per_model[[m]] <- results_score_interval
        }   

        #add model names to the list
        names(singificant_results_per_model) <-  heritage_models_test

        #save these results in the list per pheno
        singificant_results_per_pheno[[p]] <-  singificant_results_per_model    
    }
    
    #add pheno names to the list
    names(singificant_results_per_pheno) <-  significant_phenotypes

    #add these results as an element to the list per blocks
    singificant_results_per_blocks[[i]] <- singificant_results_per_pheno
}

#add block names to the list
names(singificant_results_per_blocks) <- vector_haplotype_blocks

#print results
singificant_results_per_blocks


#RESULTADOS
#$chr_20.block_1$CRF_hip$dominant ->FDR >0.05 entre niveles
#$chr_20.block_1$CRF_hip$additive ->FDR >0.05 entre niveles
#$chr_20.block_1$CRF_Body_fat_PC$dominant ->FDR >0.05 entre niveles
#$chr_20.block_1$CRF_Body_fat_PC$additive
#$chr_20.block_1$FMI$dominant ->FDR >0.05 entre niveles
#$chr_20.block_1$FMI$additive 

###Conclusiones: Solo nos quedamos con los casos que tiene P global < 0.05 y FDR<0.05 entre alguno de los niveles

    #chr_20.block_1$CRF_Body_fat_PC: TACCA tiene valores más altos de % grasa que AATTG (global p=0.0144; difference between groups 0.05; 95CI = 0.09 - 0.02; P = 0.0054; FDR=0.02721732 under additive model).
   
    #chr_20.block_1$FMI: TACCA tiene valores más altos de FMI que AATTG (global p=0.015; difference between groups 0.07; 95CI = 0.13 - 0.02; P = 0.0046; FDR=0.02279285 under additive model).

        #to see the pvalue in additive, which is very very low
        #singificant_results_per_blocks$chr_8.block_4$TG_HDL$additive$intervals["2221", "P-val"]
 
#me parece bien aprovechar esos resultados. Se me ha ocurrido un enfoque para el paper que nos puede permitir justificar aun más el uso de resultados que solo son significativos con FDR < 0.1. Podemos decir que primero buscamos asociaciones genotipo-fenotipo que fueran significativas por debajo de 0.1, y que de ahí seleccionamos los haplotipos-fenotipos utilizados en los análisis de haplótipo. Es decir, solo seleccionamos aquellos haplotipos de los que forman parte SNPs asociados con el fenotipo correspondiente con FDR < 0.1. Sería como hacer una fase inicial de "screening" ó exploración, en la cual estaría totalmente justificado subir el threshold porque luego se van a seguir haciendo análisis y por tanto un mayor riesgo de falsos positivos es aceptable. En cualquier caso, me parece bien meter esos resultados dejando claro eso que comentas, que si se baja a 0.05 se pierde significación, pero aun así pueden ser resultados relevantes.  

#Os paso los resultados de los análisis de haplotipo. Como os acabo de comentar, he mirado asociación entre haplotipos-fenotipos que estén correlacionados en los análisis preliminares con FDR < 0.1. Si no os parece bien el enfoque corro análisis entre todos los haplotipos y fenotipos, aunque no creo que salte nada nuevo. El paquete con el que estoy corriendo los análisis de haplotipo solo me deja usar modelo aditivo, dominante y recesivo. El recesivo no lo veo recomendable ya que requiere que los individuos sean homocigotos para cada haplotipo, así que los haplotipos con poca frecuencia no se pueden analizar, y nosotros tenemos muchos poco frecuentes. Por tanto solo he usado aditivo y dominante.

#He usado haploview para formar los haplotipos. El algoritmo que he usado se basa en la idea de que zonas del genoma muy cercanas y que se heredan juntas no recombinan entre si (es poco probable), por tanto, dentro de un haplotipo no debería haber recombinación, mientras que si podría haberla entre haplotipos. Por tanto, el algoritmo busca en nuestros datos conjuntos de SNPs contiguos que no tienen signos de recombinación, esos conjuntos ó bloques serían haplotipos. Una vez hechos los bloques de SNPs he testado si las frecuencias "reales" del haplotipo se desvían de las frecuencias esperadas si los polimorfismos no estuvieran en linkage, es decir, no formasen un haplotipo. En todos los casos las frecuencias de las diferentes combinaciones se desvían significativamente de lo esperado en equilibrio (es decir, si no hubiera linkage ó haplotipo; P = 0).

#En todos los análisis se ajusta por centro, sexo, edad y BMI (BMI solo para marcadores cardiovasculares). Primero testamos la asociacion global del haplotipo con el fenotipo obteniendo un score, además se realiza un permutation test (explciar los de la reordenacion de filas proporcion de simulaciones que tienen score más alto que le observado...) del que se un p.value. Luego corremos regressiones para ver entre que niveles del haplotipo hay diferencias. 


#Creo que con esto podría quedar cerrado el paper. Los análisis de interacciones gen*gen y gen*PA supongo que cambiarán un poco con el nuevo FDR, aunque creo recordar que sobre todo en AF salían cosas raras y dado todo lo que tenemos ya, lo dejaría en todo lo caso como otra opción para un posible "helena 3". Llegado el caso, si hubiera cosas interesantes se podría explorar la opción interacciones con haplotipos (ahora mismo no sé hacerlo pero creo que es posible). 

    #Se pueden hacer interacciones con haplotipos, PERO ten en cuenta que si ya tienes problemas de N con la mayoría de haplotipos (suele haber dos muy frecuentes y luego otros dos menos), imagina diviendo en dos ó más grupos. Así que creo que tenemos poco N para hacer esto. La interacción más fácil de hacer tanto por logistica como por N sería de la un factor de acitvidad física ó dieta. Las gen*gen tendría que tunear haplo.stat (MIRA 8.6 Creating Haplotype Effect Columns: haplo.design del tutorial de haplstat): https://cran.r-project.org/web/packages/haplo.stats/vignettes/manualHaploStats.pdf. 




##################################################
##### Gene - gene interaction across SNPs ########
##################################################

#In general, we have too much imbalance between groups, because we have too many groups for our sample size. Maybe we could control from the minimum sample size just like we did for crude associations and physical activity, but not sure how many association would have enough sample size. 




#########################################################################
################ GENE - ENVIRONMENT (PHYISCAL ACTIVITY) #################
#########################################################################

#pheno to model
pheno_to_model_interact = unique(geno_pheno_results_fdr_0.1$selected_pheno)

#models
heritage_models = c("dominant", "recessive", "overdominant", "codominant", "additive")

#function to model
#NOT REVISED FOR CONTINUOUS ENVIRONMENTAL VARIABLES
#to debug run the next loop with some value for p and then type: pheno_selected=selected_pheno; family=family; transformation=transformation; env_variable="PA_factor"; env_var_factor=TRUE; correction_variables = paste("+", control_variables, sep="")
siginficant_interactions_env = function(pheno_selected, family, transformation, env_variable, env_var_factor, correction_variables){

    #SE function to plot erro bars in bar plots
    se <- function(x) sd(x)/sqrt(length(x)) 

    #for discrete phenotypes (obesity) 
    if(family=="binomial"){
        family_error = "binomial"
    }

    #for continuous phenotypes
    if(family=="gaussian"){
        family_error = "gaussian"           
    }

    #apply transformation if apply
    if(transformation == "raw"){
        
        #if there is no transformation, the phenotype name is the same
        response = pheno_selected
    } else {

        #if there is transformation and it is not squared
        if(!transformation == "^2"){

            #add the transformation before the phenotype name
            response = paste(transformation, "(", pheno_selected, ")", sep="")
        } else {

            #add the transformation after the phenotype name
            response = paste("(", pheno_selected, ")", transformation, sep="")
        }           
    }

    #data frame for saving results
    final_output_results = data.frame(pheno_selected=NA, selected_model=NA, snps=NA, min_n=NA, p.val=NA, BF=NA, FDR=NA, d2_mine=NA, adjust_d2_mine=NA, d2_mumin=NA, adjust_d2_mumin=NA)
    
    #create a vector with the heritage models considered    
    heritage_models = c("dominant", "recessive", "overdominant", "codominant", "additive")

    #for each heritage model 
    for(j in 1:length(heritage_models)){

        #select the heritage model
        selected_model = heritage_models[j]

        #open empty vectors to save results
        d2_mine=NULL
        adjust_d2_mine=NULL
        d2_mumin=NULL
        adjust_d2_mumin=NULL
        snps = NULL
        min_n=NULL        
        p.val = NULL

        #calculate p.vale of interaction between the environmental variable and each snp under the [j] heritage model      
        for (i in 1:length(labels(myData_ptpn1))){ 
            
            #select the snp
            snp = labels(myData_ptpn1)[i]

            #extract genotype data of this SNP
            geno_data = eval(parse(text=paste(selected_model, "(na.omit(myData_ptpn1$", snp, "))", sep="")))#we consider the whole dataset because we want to extract all the levels included in this cohort for the variable. If any of this levels has low number of individuals for a given conbination with the environmental variable, this will detected in in the next lines (sample_size_per_geno). If there is only one level across the whoe cohort, we don't want to run the analysis.

            #extract genotype levels
            levels_genotypes = unique(geno_data)

            #number of levels of the [i] snp
            n_levels_snp = length(levels_genotypes)

            #calculate the number of individuals with each genotype and with data for the [p] phenotype
            if(env_var_factor){ #if the environmental variable is a factor, we will also consider the levels of this variable. 

                #levels of the environmental variable
                env_variable_data = eval(parse(text=paste("na.omit(myData_ptpn1$", env_variable, ")", sep="")))#we consider the whole dataset because we want to extract all the levels included in this cohort for the variable. If any of this levels has low number of individuals for a given conbination with the SNP, this will detected in in the next lines (sample_size_per_geno). If there is only one level across the whoe cohort, we don't want to run the analysis.

                #extract genotype levels
                levels_env_variable = sort(unique(env_variable_data))

                #number of levels of the environmental variable
                n_levels_env = length(levels_env_variable)

                #if the family error is NOT binomial and hence the phenotype is continuous
                if(!family_error == "binomial"){

                    #levels of the discrete phenotype: It is zero because the phenotype is continuous
                    n_levels_discre_pheno = 0

                    #for each level of the environmental variable
                    sample_size_per_geno = NULL                
                    for(e in 1:length(levels_env_variable)){                       
                        
                        #select the [e] level
                        selected_level_env = levels_env_variable[e]

                        #for each genotype
                        for(l in 1:length(levels_genotypes)){
                            
                            #select the [l] genotype
                            selected_level_snp = levels_genotypes[l]

                            #select those rows with the [l] level of the genotype and without na for the select phenotype
                            subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", snp, ")=='", selected_level_snp, "' & myData_ptpn1$", env_variable, "=='", selected_level_env, "' & !is.na(myData_ptpn1$", pheno_selected, ")),]", sep="")))

                            #save the number of rows
                            sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))
                        }
                    } 
                } else { #if not, and then the phenotype is discrete

                    #extract the levels of the phenotype
                    levels_factor = sort(unique(eval(parse(text=paste("myData_ptpn1$", selected_pheno, sep="")))))

                    #levels of the discrete phenotype
                    n_levels_discre_pheno = length(levels_factor)

                    #for each level of the environmental variable
                    sample_size_per_geno = NULL                
                    for(e in 1:length(levels_env_variable)){                       

                        #select the [e] level
                        selected_level_env = levels_env_variable[e]

                        #for each genotype
                        for(l in 1:length(levels_genotypes)){
                            
                            #select the [l] genotype
                            selected_level_snp = levels_genotypes[l]

                            #for each level of the factor phenotype
                            for(f in 1:length(levels_factor)){

                                #select the [f] level of the factor
                                selected_level_factor = levels_factor[f]
                            
                                #select those rows with the [l] level of the genotype and without na for the select phenotype
                                subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", snp, ")=='", selected_level_snp, "' & myData_ptpn1$", env_variable, "=='", selected_level_env, "' & myData_ptpn1$", pheno_selected, "==", selected_level_factor, "),]", sep="")))

                                #save the number of rows
                                sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))
                            }
                        }
                    }                    
                }             
            } else { #NOT REVISED for non-factor environmental variable - NOT REVISED

                #number of levels of the environmental variable: We have no env variable factor
                n_levels_env = 0

                #if the variable is continuous
                if(!family_error == "binomial"){

                    #levels of the discrete phenotype: It is zero because the phenotype is continuous
                    n_levels_discre_pheno = 0

                    #for each genotype
                    sample_size_per_geno = NULL                 
                    for(l in 1:length(levels_genotypes)){
                        
                        #select the [l] genotype
                        selected_level_snp = levels_genotypes[l]

                        #select those rows with the [l] level of the genotype and without na for the select phenotype
                        subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", snp, ")=='", selected_level_snp, "' & !is.na(myData_ptpn1$", env_variable, ") & !is.na(myData_ptpn1$", pheno_selected, ")),]", sep="")))

                        #save the number of rows
                        sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))
                    }  
                } else {#if not and thus is a discrete variable

                    #extract the levels of the phenotype
                    levels_factor = sort(unique(eval(parse(text=paste("myData_ptpn1$", selected_pheno, sep="")))))

                    #levels of the discrete phenotype
                    n_levels_discre_pheno = length(levels_factor)

                    #for each genotype
                    sample_size_per_geno = NULL                 
                    for(l in 1:length(levels_genotypes)){
                        
                        #select the [l] genotype
                        selected_level_snp = levels_genotypes[l]

                        #for each level of the discrete phenotype
                        for(f in 1:length(levels_factor)){

                            #select the [f] level of the factor
                            selected_level_factor = levels_factor[f]

                            #select those rows with the [l] level of the genotype and without na for the select phenotype
                            subset_geno_no_na = eval(parse(text=paste("myData_ptpn1[which(", selected_model, "(myData_ptpn1$", snp, ")=='", selected_level_snp, "' & !is.na(myData_ptpn1$", env_variable, ") & myData_ptpn1$", pheno_selected, "=='", selected_level_factor, "'),]", sep="")))
                            
                            #save the number of rows
                            sample_size_per_geno =  append(sample_size_per_geno, nrow(subset_geno_no_na))
                        }
                    }
                }                             
            }
                        
            #If we are analyzing SNPs with MAF < 0.1, we only want the dominant model, avoid the other models. Moreover, if the snp has only one genotpye in general or in a level of the environmental variable (in case is a factor), or in general all individuals have the same genotype skip this model
            if(selected_model %in% c("recessive", "overdominant", "codominant", "additive") & !snp %in% snps_maf_low_0.9 | TRUE %in% (sample_size_per_geno < 10)  | n_levels_snp==1 | n_levels_env == 1 | n_levels_discre_pheno == 1){

                #save all as NA
                d2_mine=append(d2_mine, NA)
                adjust_d2_mine=append(adjust_d2_mine, NA)
                d2_mumin=append(d2_mumin, NA)
                adjust_d2_mumin=append(adjust_d2_mumin, NA)
                min_n=append(min_n, NA)
                p.val = append(p.val, NA)
                snps = append(snps, snp) 
            } else { #if not, we can model 
                
                #test significant reductions of deviance with and without interaction
                model1 = glm(formula(paste(response, "~", selected_model, "(", snp, ")*", env_variable, correction_variables, sep="")), data=myData_ptpn1, family=family_error)
                model2 = glm(formula(paste(response, "~", selected_model, "(", snp, ")+", env_variable, correction_variables, sep="")), data=myData_ptpn1, family=family_error)
                    #we do not need to make a subset of the data specifying that we only want rows with data for the SNP, because in this case, the second model (nested or null) also includes the SNP. The difference between the two models is the interaction.  
                
                #extract the p-value
                results = anova(model1, model2, test="Chi")


                ##calculate the R2. WE WANT the R2 of the interaction after adjusting by the controlling variables and also the SNP and the environmental variable.
                #first using the modified function of Nick Zimmermann
                d2_nick = Dsquared_mod(model=model1, null_model = model2, adjust=FALSE)
                adjust_d2_nick = Dsquared_mod(model=model1, null_model = model2, adjust=TRUE, n_levels_predictor_test=1, predictor_is_factor=FALSE) #we calculate the adjusted R2 and for that we need to set TRUE for adjust. We set the number of levels of the predictor as 1, because the difference between models is only the interaction, but both the environmental variable and the SNP are maintained in the nested model. Indeed, if you check the difference in degrees of freedom between models, in this case it is 1. We indicate that the predictor is not a factor. This is TRUE for genotypes, when we remove a SNP from the model, if you have 3 levels, the decrease in df is 3-1=2 (number of genotypes respect to the reference genotype). In this case, we are not removing a factor, df is 1. 

                #second with the MuMIn package. This statistic is is one of the several proposed pseudo-R^2's for nonlinear regression models. It is based on an improvement from _null_ (intercept only) model to the fitted model, and calculated as R^2 = 1 - exp(-2/n * logL(x) - logL(0)) where logL(x) and logL(0) are the log-likelihoods of the fitted and the _null_ model respectively. ML estimates are used if models have been fitted by REstricted ML (by calling ‘logLik’ with argument ‘REML = FALSE’). Note that the _null_ model can include the random factors of the original model, in which case the statistic represents the ‘variance explained’ by fixed effects. For OLS models the value is consistent with classical R^2. In some cases (e.g. in logistic regression), the maximum R_LR^2 is less than one.  The modification proposed by Nagelkerke (1991) adjusts the R_LR^2 to achieve 1 at its maximum: Radj^2 = R^2 / max(R^2) where max(R^2) = 1 - exp(2 / n * logL(0)) . ‘null.fit’ tries to guess the _null_ model call, given the provided fitted model object. This would be usually a ‘glm’. The function will give an error for an unrecognised class.
                d2_mumin_raw = r.squaredLR(object=model1, null=model2, null.RE=FALSE)
                    #object: a fitted model object (a.k.a. the full model)
                    #null: a fitted null model. This is the model to compare with. In my case, the simpler and nested model without the SNP under study.
                    #null.RE: logical, should the null model contain random factors?  Only used if no _null_ model is given, otherwise omitted, with a warning.
                    #the rest of arguments are for fitting the null model, but we have already fit it. 
                    #the result is: ‘r.squaredLR’ returns a value of R_LR^2, and the attribute ‘"adj.r.squared"’ gives the Nagelkerke's modified statistic. Note that this is not the same as nor equivalent to the classical ‘adjusted R squared’. IT IS NOT THE SAME ADJUST THAN THE ZIMMERMAN ADJUSTMENT
                #extract d2 and adjust d2 from the mumin results
                d2_mumin_pkg = d2_mumin_raw[1]
                adjust_d2_mumin_pkg = attributes(d2_mumin_raw)$adj.r.squared

                #save the p.val and the snp name
                d2_mine=append(d2_mine, d2_nick)
                adjust_d2_mine=append(adjust_d2_mine, adjust_d2_nick)
                d2_mumin=append(d2_mumin, d2_mumin_pkg)
                adjust_d2_mumin=append(adjust_d2_mumin, adjust_d2_mumin_pkg)
                min_n=append(min_n, min(sample_size_per_geno))                
                p.val = append(p.val, results[,which(names(results)=="Pr(>Chi)")][2])
                snps = append(snps, snp)
            }
        }
  
        #calculate FDR
        pval_to_FDR<-p.val
        FDR = p.adjust(pval_to_FDR, method="BH")

        #indicate if there is significance after bonferroni correction
        bonferroni_significance = 0.05/length(na.omit(p.val))
        BF = NULL
        for (i in 1:length(pval_to_FDR)){
            selected_pval = pval_to_FDR[i]
            if(!is.na(selected_pval)){
                if(selected_pval < bonferroni_significance){
                    BF = append(BF, "YES")
                } else { 
                    BF = append(BF, "NO")                
                } 
            } else {
                BF = append(BF, NA) 
            }       
        }  

        #bind all results
        final_results = cbind.data.frame(rep(pheno_selected, length(snps)), rep(selected_model, length(snps)), snps, min_n, p.val, BF, FDR, d2_mine, adjust_d2_mine, d2_mumin, adjust_d2_mumin) 
        colnames(final_results)[1] <- "pheno_selected"
        colnames(final_results)[2] <- "selected_model"
        colnames(final_results)[3] <- "snps"
        colnames(final_results)[4] <- "min_n"
        colnames(final_results)[5] <- "p.val"
        colnames(final_results)[6] <- "BF"                
        colnames(final_results)[7] <- "FDR"
        colnames(final_results)[8] <- "d2_mine"
        colnames(final_results)[9] <- "adjust_d2_mine"
        colnames(final_results)[10] <- "d2_mumin"
        colnames(final_results)[11] <- "adjust_d2_mumin"

        #save it
        final_output_results = rbind.data.frame(final_output_results, final_results)                 
    }

    #remove first row with NA
    final_output_results = final_output_results[-which(rowSums(is.na(final_output_results)) == ncol(final_output_results)),]

    #save it
    return(final_output_results)
} #NOT REVISED FOR CONTINUOUS ENVIRONMENTAL VARIABLES
 
#with FDR<0.05
interact_PA_results = data.frame(pheno_selected=NA, selected_model=NA, snps=NA, min_n=NA, p.val=NA, BF=NA, FDR=NA, d2_mine=NA, adjust_d2_mine=NA, d2_mumin=NA, adjust_d2_mumin=NA)
#for each phenotype from significant phenotypes (pheno_to_model_interact)
for(p in 1:length(pheno_to_model_interact)){

    #select the [p] phenotype
    selected_pheno = pheno_to_model_interact[p]

    #if the [p] phenotype is not one of the binomial phenotypes
    if(!selected_pheno %in% c("CVi_BP", "obesity")){

        #the family for the glm is gaussian
        family="gaussian"

        #and the transformation is log
        transformation="log"
    } else{#if not and hence the phenotype is binomial

        #the family for the glm is binomial
        family="binomial"

        #and there is no transformation (i.e., raw)
        transformation="raw"
    }

    #select control variables: BMI is a control variable only for CVD variables
    if(!selected_pheno %in% c("CVi_BP", "SBP", "DBP", "TG", "TC", "LDL", "HDL", "LDL_HDL", "Apo_A1", "Apo_B", "ApoB_ApoA1", "apoB_LDL", "TG_HDL", "Insulin", "Leptin_ng_ml", "HOMA", "QUICKI",  "TC_HDL",  "risk_score_for_log")){
        
        #not add BMI as a covariable, because the phenotype is an adiposity marker
        control_variables = "CRF_sex+CRF_age+center"
    } else{
        
        #add also BMI as a covariable
        control_variables = "CRF_BMI+CRF_sex+CRF_age+center"         
    }

    #run the analyses
    interact_PA_results = rbind.data.frame(interact_PA_results, siginficant_interactions_env(pheno_selected=selected_pheno, family=family, transformation=transformation, env_variable="PA_factor", env_var_factor=TRUE, correction_variables = paste("+", control_variables, sep="")))
}

#remove the first row with NA
interact_PA_results = interact_PA_results[-which(rowSums(is.na(interact_PA_results)) == ncol(interact_PA_results)),]

#check we have all the interaction we should have
nrow(interact_PA_results) == length(pheno_to_model_interact) * length(heritage_models) * length(labels(myData_ptpn1))

#see the results
interact_PA_results


####### select those interactions for wich genotypes and phenotypes are previously associated according to FDR<0.1 ########
#Combinations of phenotypes*snps that interact with AF
significan_snp_pheno_interact = interaction(interact_PA_results$pheno_selected, interact_PA_results$snps)
#check
summary(significan_snp_pheno_interact == paste(interact_PA_results$pheno_selected, ".", interact_PA_results$snps, sep=""))

#Combinations of phenotypes*snps that are significant for crude associaitons
significan_snp_pheno_crude #created at the end of geno-pheno asocciation crudes (no haplotype)

#select those interactions
interact_PA_results_clean = interact_PA_results[which(significan_snp_pheno_interact %in% significan_snp_pheno_crude),]
interact_PA_results_clean

#number of interaction below an FDR of 0.05
length(which(interact_PA_results_clean$FDR<0.05)) #this is exactly the number of interactions in the figures of the main text of the previous versions of the manuscript. 

#check that these interactions implie snps and phenotypes previously associated
summary(interaction(interact_PA_results_clean$pheno_selected, interact_PA_results_clean$snps) %in% significan_snp_pheno_crude)

#check that the interactions not included are not inside the significan_snp_pheno_crude
#extract not included interactios
discarded_interactions = interact_PA_results[which(!significan_snp_pheno_interact %in% significan_snp_pheno_crude),]
#extract their pheno*snp combinations and see if they are inside significan_snp_pheno_crude
summary(!paste(discarded_interactions$pheno_selected, ".", discarded_interactions$snps, sep="") %in% significan_snp_pheno_crude)


##Comparisons R2
#R2 mine and with mumin: They are very similar.
summary(interact_PA_results_clean$d2_mine - interact_PA_results_clean$d2_mumin) 
plot(interact_PA_results_clean$d2_mine, interact_PA_results_clean$d2_mumin)
cor.test(interact_PA_results_clean$d2_mine, interact_PA_results_clean$d2_mumin, method="spearman") #rho=0.999
#there is a case that it is a little bit different:
interact_PA_results_clean[which(abs(interact_PA_results_clean$d2_mine - interact_PA_results_clean$d2_mumin) == max(na.omit(abs(interact_PA_results_clean$d2_mine - interact_PA_results_clean$d2_mumin)))),]
    #It is with obesity, but there is other association with obesity that have more similar R2 with both methods, so I do not think this is a question about logistic. In addition, the correction for non-linear GLMs is done in adjusted mumin. The difference is very small, from 0.009939082745 to 0.009320662592. 0.009939082745-0.009320662592 = 0.000618420153. (0.000618420153/0.009939082745)*100=6.22% of difference. It is only one case and very small difference. The difference in the supple data will be 0.93% instead of 0.99.
        #which((interact_PA_results_clean$d2_mine - interact_PA_results_clean$d2_mumin) == max(abs(na.omit(interact_PA_results_clean$d2_mine - interact_PA_results_clean$d2_mumin))))
dev.off()

#R2 adjusted is not equal with both methods.
summary(interact_PA_results_clean$adjust_d2_mine - interact_PA_results_clean$adjust_d2_mumin) 
plot(interact_PA_results_clean$adjust_d2_mine, interact_PA_results_clean$adjust_d2_mumin)
cor.test(interact_PA_results_clean$adjust_d2_mine, interact_PA_results_clean$adjust_d2_mumin, method="spearman") #rho=0.999
    #The non-adjusted values are similar between methods, but not when we adjust. Note that the method for adjusting used in the MuMIn function is different from the typical adjust based on the number of parameters and sample size. Mumin uses the Nagelkerke adjustment (see below). I have made the adjust manually modifying the function of Zimmermann.
dev.off()

#For OLS models (linear regression?) the value is consistent with classical R^2. In some cases (e.g. in logistic regression), the maximum R_LR^2 is less than one.  The modification proposed by Nagelkerke (1991) adjusts the R_LR^2 to achieve 1 at its maximum: Radj^2 = R^2 / max(R^2). So you are basically adjusting the R2 by the maximum R2 present in your data. I see the point because the R2 max is not 1 in logistic. I see the point, you could underestimate the R2 for logistic models, but I do not fully understand how the maximum R2 can be established for a given model comparison.

#there are differences between adjusted and non-adjusted in mumin, but this is not caused by the logistic
plot(interact_PA_results_clean$adjust_d2_mumin, interact_PA_results_clean$d2_mumin, col="red")
cor.test(interact_PA_results_clean$adjust_d2_mumin, interact_PA_results_clean$d2_mumin)
dev.off()

#you can see here how when considering only obesity (factor with 2 levels, that is, logistic), there is a perfect correlation between adjusted and un-adjusted in mumin. 
plot(interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected == "obesity"),]$adjust_d2_mumin, interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected == "obesity"),]$d2_mumin, col="red")
#cor.test(interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected == "obesity"),]$adjust_d2_mumin, interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected == "obesity"),]$d2_mumin, col="red")
dev.off()

#Given that I am not 100% sure if this a correct way to adjust and R2 does NOT seem underestimated for obesity (logistic), I am going to use the traditional pseudo R2 of glms, which is similar with the formula of Zimmermamnn and the function of MuMIn. In addition, my d2 and the mumin d2 are not so different respect to adjusted d2 using the traditional approach
summary(interact_PA_results_clean$d2_mine-interact_PA_results_clean$adjust_d2_mine)
plot(interact_PA_results_clean$d2_mine, interact_PA_results_clean$adjust_d2_mine)
summary(interact_PA_results_clean$d2_mumin-interact_PA_results_clean$adjust_d2_mine)
plot(interact_PA_results_clean$d2_mumin, interact_PA_results_clean$adjust_d2_mine)
    #There are almost no differences. Only one case (row 3) with a little bit of difference, which is the same showing some differences between mummin and my R2 without adjusting.
dev.off()

#Given that both R2 (mine and mumin) are VERY similar and that mumin is more reproducible because it is in a published package, we will use the R2 of Mumin. 

#Summary: We will use the classical R^2 calculated with the mumin package.


## use these results to create the supplementary data 2
#set the folder to save
folder_to_save_supple_data_2 = "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/supple_data"
system(paste("mkdir -p ", folder_to_save_supple_data_2, sep=""))
    #p: no error if existing, make parent directories as needed

#select the rows and columns we are interested
suppl_data_2 = interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected != "risk_score_for_log"), which(colnames(interact_PA_results_clean) %in% c("pheno_selected", "selected_model", "snps", "min_n", "p.val", "FDR", "d2_mumin"))]
    #We remove all rows belonging to the CVD risk score because this variable was finally not used in the manuscript. 
        #IMPORTANT: If we include risk_score_for_log, then we should change to risk_score in the supple data.     
    #We select the columns that includes the variables selected for the supplementary dataset 1

#convert R2 from 0-1 to percentage
suppl_data_2$d2_mumin = (suppl_data_2$d2_mumin*100)/1
    #For example: If over 1, we have 0.5, over 100 we would have X. X being (100*0.5)/1=50 -> 50% 

#change columns names
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "pheno_selected")] <- "phenotype"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "selected_model")] <- "heritage_model"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "snps")] <- "snp"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "min_n")] <- "min_sample_size"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "p.val")] <- "p_value"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "FDR")] <- "fdr"
colnames(suppl_data_2)[which(colnames(suppl_data_2) == "d2_mumin")] <- "r2_percentage"

#set the names of the files for saving
suppl_data_2_file_name = "suplementary_data_2_v2.csv"
suppl_data_2_file_name_zip = "suplementary_data_2_v2.zip"

#save the table
write.table(suppl_data_2, paste(folder_to_save_supple_data_2, "/", suppl_data_2_file_name, sep=""), col.names=TRUE, row.names=FALSE, sep=",")

#compress the text file and remove it after compression
system(paste("cd ", folder_to_save_supple_data_2, "; rm ", suppl_data_2_file_name_zip, "; zip ", suppl_data_2_file_name_zip, " ", suppl_data_2_file_name, " ; rm ", suppl_data_2_file_name, sep=""))
    #we could save directly as ".gz" using gzfile() around the file path with write.table. I avoid this option because this could give problems with the reviewers and readers, because they have to use a third party software to open it in windows.


## compare this version of the supplementary with the first one.

#Initially, I made some changes in this version for recoding the a problematic SNP respect allele names between HELENA and ncbi. But after thinking, I discover a more minimalistic way to solve this problem without touching the dataset.

#read the file of the previous version
suppl_data_2_previous_version = read.table(paste(folder_to_save_supple_data_2, "/suplementary_data_2_v1.txt.gz", sep=""), sep="\t", header=TRUE)

#read the file of the last version
suppl_data_2_last_version = read.table(unz(paste(folder_to_save_supple_data_2, "/suplementary_data_2_v2.zip", sep=""), suppl_data_2_file_name), sep=",", header=TRUE)

#select the non-numeric columns
no_numeric_cols = which(colnames(suppl_data_2_last_version) %in% c("phenotype", "heritage_model", "snp"))

#check that the non-numeric and numeric columns are the same between versions
summary(suppl_data_2_previous_version[,no_numeric_cols] == suppl_data_2_last_version[,no_numeric_cols])
summary(suppl_data_2_previous_version[,-no_numeric_cols] - suppl_data_2_last_version[,-no_numeric_cols])



###### Correlaton genotypes - AF 
#we don't include in the paper the correlation between AF and genotypes because we have filtered by imbalances between levels. If a genotype (e.g. G/G) has more individuals with high MVPA, we would have more G/G individuals in the level of high PA than in the low level. If most GG individuals have high MVPA, the level of low PA would have very few GG. Therefore, we would have a imbalance of levels. When genotypes and AF are correlated, we don't have a similar distribution sample size of genotypes between the two levels of physical activity. We have applied filters for that, eliminating cases with levels with low number of individuals. In these cases with imbalance it is difficult to test an interaction because you don't have enough individuals with the same genotype under the two levels of physical activity. In other words, if you have a great correlation between genotypes and AF, you should have great differences of sample size between genotypes across physical activity levels.

#heritage models
heritage_models = c("dominant", "recessive", "overdominant", "codominant", "additive")

#empty data frame to save results
results_cor_af_geno = data.frame(selected_snp=NA, selected_model=NA, pval=NA)

#for each snp
for(i in 1:length(labels(myData_ptpn1))){

    #select the [i] snp
    selected_snp = labels(myData_ptpn1)[i]

    #for each heritage model
    for(m in 1:length(heritage_models)){

        #select the [m] model
        selected_model = heritage_models[m]
    
        #extract data without NA for the selected snp under the selected model
        data = eval(parse(text=paste("myData_ptpn1[which(!is.na(", selected_model, "(myData_ptpn1$", selected_snp, "))),]", sep="")))

        #fit the models
        model1 = glm(paste("MVPA_mean~", selected_model, "(", selected_snp, ")", sep=""), data=data)
        model2 = glm(paste("MVPA_mean~1", sep=""), data=data)

        #calculate an anova
        results = anova(model1, model2, test="Chi")

        #extract the pval
        pval = results[which(names(results) == "Pr(>Chi)")][2,]

        #save the snp, the model and the pval
        results_cor_af_geno = rbind.data.frame(results_cor_af_geno, cbind.data.frame(selected_snp, selected_model, pval))
    }    
}
#remove the first row with NAs
results_cor_af_geno = results_cor_af_geno[-which(rowSums(is.na(results_cor_af_geno)) == ncol(results_cor_af_geno)),]

#check that all models and snps have been included
nrow(results_cor_af_geno) == length(labels(myData_ptpn1))*length(heritage_models)

#see results
results_cor_af_geno

#select significant cases
significant_cor_geno_af = results_cor_af_geno[which(results_cor_af_geno$pval < 0.05),]

#check that the significant interactions are not included in the cases of correlation between genotype and AF
#combination of snp*model with significant interactions
snp_model_significant_interac = interaction(interact_PA_results_clean$snps, interact_PA_results_clean$selected_model)
#check
snp_model_significant_interac == paste(interact_PA_results_clean$snps, interact_PA_results_clean$selected_model, sep=".")

#combination of snp*model associated with AF
snp_model_significant_assoc_af = interaction(significant_cor_geno_af$selected_snp, significant_cor_geno_af$selected_model)
#check
snp_model_significant_assoc_af == paste(significant_cor_geno_af$selected_snp, significant_cor_geno_af$selected_model, sep=".")

# significant interactions are not included in the snp*models with correlation with AF?
!unique(snp_model_significant_interac) %in% snp_model_significant_assoc_af

#see problematic cases. SNP-model combinations that show significant interactions with PA
significant_cor_geno_af[which(snp_model_significant_assoc_af %in% snp_model_significant_interac),]
    #snp_model_significant_assoc_af comes from significant_cor_geno_af, so we can use it to select rows in significant_cor_geno_af.

#CONCLUSIONS: there are 5 p.values lower than 0.05, but all of them are higher than 0.01. These are relatively high Pvalues, and only affects two snps: rs2143511 and rs6020608. rs2143511 does not show great imbalances, in the case of rs6020608 there is more imbalance (a group with 30), but this is not very bad. In addition, there are very significant and explicative associations beside these SNPs (R2 higher 1.5% in two cases), suggesting that our signals is not caused by a correlation between AF and genotypes. 




#####################################################################################
############# TEST FOR THE EXISTENCE OF ONLY ONE CAUSAL VARIANT #####################
#####################################################################################

#We are going to test whether each SNP explains variance of the phenotypes independently of the rest of SNPs. For that, we are going to calculate a full model with the 5 SNPs, and then remove each SNP and calculate a likelihood ratio test between full and nested model. In that way, we can test whether removing a SNP increases the non-explained variance. We are also going to calculate how explicative is each SNP while considering the rest of them, so we will calculate the R2 between the full model and the model without each SNP. Finally, we will show the explicative power of the full model and the haplo.glm respect to a model with only the covariates. 

#This approach is probably better than the stepwise regression because we want to know how explicate is each SNP respect to the rest of them. In the stepwise regression, the less explicative SNP is removed, checked the model, then the next less explicative SNP, and so on... until the next model is worse. 

#Using this approach, we can see whether a SNP is explaining variance of the phenotypes independently of the rest. In that case, the removal of that SNP from the full model should increase the non-explained variance (P<0.05) and the R2 of that SNP should be close to the R2 of the full and haplotype models. 

#pheno to model are those with FDR<0.1, because these are the ones used in subsequent analyses
pheno_to_model_causal_variant = unique(geno_pheno_results_fdr_0.1$selected_pheno)

#models to test
models_causal_variant = c("codominant", "dominant", "recessive", "overdominant", "additive")

#create the seed for the final dataset with results
causal_variants_results = data.frame(selected_pheno=NA, selected_model=NA, model=NA, p_value=NA, r2=NA)

#for each pheno
for(p in 1:length(pheno_to_model_causal_variant)){

    #select the [p] phenotype
    selected_pheno = pheno_to_model_causal_variant[p]

    #select the family and transformation
    if(!selected_pheno %in% c("obesity", "CVi_BP")){
        family="gaussian"
        transformation="log("
    } else{
        family="binomial"
        transformation="("
    }

    #select control variables: BMI is a control variable only for CVD variables
    if(!selected_pheno %in% c("CVi_BP", "SBP", "DBP", "TG", "TC", "LDL", "HDL", "LDL_HDL", "Apo_A1", "Apo_B", "ApoB_ApoA1", "apoB_LDL", "TG_HDL", "Insulin", "Leptin_ng_ml", "HOMA", "QUICKI",  "TC_HDL",  "risk_score_for_log")){
        control_variables = "CRF_sex+CRF_age+center"
    } else{
        control_variables = "CRF_BMI+CRF_sex+CRF_age+center"         
    }

    #open an empty data.frame to save results using as model the dataframe for saving final results
    causal_variants_results_raw = data.frame(matrix(NA, nrow=1, ncol=ncol(causal_variants_results)))
    colnames(causal_variants_results_raw) <- colnames(causal_variants_results)

    #for each model
    for(m in 1:length(models_causal_variant)){

        #select the [m] model
        selected_model = models_causal_variant[m]

        #select those SNPs that has been analyzed for the [p] phenotype and the [m] model
        snp_to_include = geno_pheno_results[which(geno_pheno_results$selected_pheno == selected_pheno & geno_pheno_results$selected_model == selected_model & !is.na(geno_pheno_results$pvals)),]$snp_to_test
            #we want SNPs analyzed, with a P-value, independently of its significance. SNPs with no P-value did not match our requeriments of sample size and were discarded in gene_pheno analyses.
            #selecting these snps, we are now sure that the analyzed snps are not disbalanced (minor homozigotes are at least 10).
            #We do not have problems of sample size for adding the 5 SNPs WITHOUT interactions between SNPs
                #Each SNP has a high genotyping success (at least 0.99), so we have genotyping data for most of individuals.
                #Sample size is around 1000 thousand, assuming that we need at least 10 individuals per group (zimmermann), we would the following sample size in the most extreme case:
                    #5 predictors (SNPs) with 3 levels (genotypes): 5*3*10 = 150 individuals
                    #1 predictor (center) with 10 levels (centers): 1*10*10 = 100 individuals
                    #1 predictor (sex) with two levels (sexes): 1*2*10 = 20 individuals
                    #2 continuous variables (age and BMI): 2*1*10 = 20 individuals
                    #TOTAL: 290 individuals. 
                #we are fine, the problem would be adding interactions between SNPs, because we are breaking down our cohort in several groups (up to 6), increasing the probability to get very small groups, specially for minor homozygotes.
            #The only problem can be that we are including more than 3-4 predictors, and then it is difficult to asses deviations from normality with residuals, and some predictors can be correlated. This can decrease our power to detect the optimal solution.

        #
        selected_snps_with_model = paste(selected_model, "(", snp_to_include, ")", collapse="+", sep="")

        is_na_snps = paste("!is.na(myData_ptpn1$", snp_to_include, ")", collapse=" & ", sep="")

        subset_snps_no_na = eval(parse(text=paste("myData_ptpn1[which(", is_na_snps, "),]", sep="")))


        model_global_1 = glm(paste(transformation, selected_pheno, ") ~ ", selected_snps_with_model, "+", control_variables, sep=""), data=subset_snps_no_na, family=family)

        model_global_2 = glm(paste(transformation, selected_pheno, ") ~ ", control_variables, sep=""), data=subset_snps_no_na, family=family)


        r2_global = r.squaredLR(object=model_global_1, null=model_global_2, null.RE=FALSE)[1]

        causal_variants_results_raw = rbind.data.frame(causal_variants_results_raw, cbind.data.frame(selected_pheno, selected_model, model="full", p_value=NA, r2=r2_global))


        #perturb?
            #maybe we can get less R2 because of multicolinearity
        
        if(selected_model %in% c("dominant", "additive")){
            r2_haplo = NULL
            for(j in 1:length(vector_haplotype_blocks)){
                selected_haplo_block = vector_haplotype_blocks[j]
                result_haplo_glm = eval(parse(text=paste("haplo_list_per_block$", selected_haplo_block, "$", selected_pheno, "$regression$", selected_model, sep="")))
                r2_haplo = append(r2_haplo, r.squaredLR(object=result_haplo_glm, null=model_global_2, null.RE=FALSE)[1])
            }

            #CHECK VERY GOOD WHERE YOU GET THE HAPLO GLM AND IF YOU CAN COMPARE IT WITH MODEL GLOBAL 2 (MODEL WITH ONLY COVARIALBES)
            
            r2_haplo = data.frame(r2_haplo)
            colnames(r2_haplo) <- vector_haplotype_blocks

        } else {

            #create an empty data.frame for R2 of each haplotype block
            empty_df_haplos = data.frame(matrix(NA, 1, length(vector_haplotype_blocks)))
            #set colnames
            colnames(empty_df_haplos) <- vector_haplotype_blocks

            r2_haplo = empty_df_haplos
                #empty_df_haplos comes form before the loop
        }

        causal_variants_results_raw = rbind.data.frame(causal_variants_results_raw, cbind.data.frame(selected_pheno, selected_model, model="haplo", p_value=NA, r2=as.numeric(r2_haplo)))


        #for each model               
        for(k in 1:length(snp_to_include)){
    
            #select the [k] snp
            selected_snp = snp_to_include[k]


            selected_snps_without_tested = paste(selected_model, "(", snp_to_include[which(selected_snp != snp_to_include)], ")", collapse="+", sep="")

            model_global_3 = glm(paste(transformation, selected_pheno, ") ~ ", selected_snps_without_tested, "+", control_variables, sep=""), data=subset_snps_no_na, family=family)
            

            p_value_without_snp = anova(model_global_1, model_global_3, test="Chi")$"Pr(>Chi)"[2]

            r2_without_snp = r.squaredLR(object=model_global_1, null=model_global_3, null.RE=FALSE)[1]


            causal_variants_results_raw = rbind.data.frame(causal_variants_results_raw, cbind.data.frame(selected_pheno, selected_model, model=selected_snp, p_value=p_value_without_snp, r2=r2_without_snp))
        }        
    }

    causal_variants_results_raw = causal_variants_results_raw[-which(rowSums(is.na(causal_variants_results_raw)) == ncol(causal_variants_results_raw)),]

    causal_variants_results = rbind.data.frame(causal_variants_results, causal_variants_results_raw)
}

causal_variants_results = causal_variants_results[-which(rowSums(is.na(causal_variants_results)) == ncol(causal_variants_results)),]

#check number of rows combinations

causal_variants_results[which(causal_variants_results$selected_pheno %in% c("CRF_Body_fat_PC", "FMI") & causal_variants_results$selected_model %in% c("additive", "dominant")),]

#NO PVALUE for full y haplo because these have been already calculated previously. add?
#r2 of each SNP is different because is the R2 after accounting for the rest of SNPs

#el 11 a veces es significativo, pero tenendo un R2 muy bajo en comparación con el full model

#estos resultados sugieren que cada SNP no explica mucha variabilidad de forma independiente de los fenotipos, porque al eliminar uno a uno del full model, no hay una reducción significativa del poder explicativo (deviance?). Ademas, los R2 de cada SNP por separado, son muy bastantes más bajos que los del full and haplotype model. Todo esto apoya nuestros resutlados sobre la existencia de un haplotipo de este gen que se asocia con adiposidad en nuestra cohorte. 

#Al hablar del último SNP en la discu o tal vez en los resultados de haplotipos, se podría mencionar que de forma individual e indepndiente cada SNP parece aportar poco, y ya explicas en detalle como lo has mirado en el cuarto supple data. 

#podrías meter todos los supples de data en un unico zip, realmente son excels que pesan poco. 





##############################################################################
############# Change risk_score names in results to plot #####################
##############################################################################

#change the name of the risk_score_for_log to risk_score, to match the name of the variable that will be used for plotting. We change for all the results at the end because all the analyses have to be done. If for example, we change the name is geno-pheno, the haplotype code will try to model with log(risk_score), which has negative values. 

#In the scripts for figures, then the p-values of risk_score_for_log will be extracted by calling risk_score. The values (with positive and negative) for plotting will be obtained from the original variable. 


##change in geno-pheno
geno_pheno_results[which(geno_pheno_results$selected_pheno == "risk_score_for_log"),]$selected_pheno <- "risk_score"
#check
!"risk_score_for_log" %in% geno_pheno_results$selected_pheno


##change in interaction (we remove non-significant interactions so it's possible that risk_score is not included)
if("risk_score_for_log" %in% interact_PA_results_clean$pheno_selected){
    interact_PA_results_clean[which(interact_PA_results_clean$pheno_selected == "risk_score_for_log"),]$pheno_selected <- "risk_score"
}
#check
!"risk_score_for_log" %in% interact_PA_results_clean$pheno_selected




#################################
##### LOAD ENVIRONMENT ##########
#################################
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/rdata/analysis_v2.RData")
require(SNPassoc)
require(genetics)