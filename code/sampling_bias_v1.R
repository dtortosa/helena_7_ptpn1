#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
    #https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
    #https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file




#########################################################################################################
######### SCRIPT FOR ANALYZING THE POSSIBLE EXISTENCE OF SAMPLING BIAS WITH PHYSICAL ACTIVITY ###########
#########################################################################################################

#In this script, we test the existence of sampling bias in the subset of individuals with physical activity. 




##############################################################
######### DIFFERENCES RESPECT TO PREVIOUS VERSIONS ###########
##############################################################

#Version 1




########################################
############ DATA PREPARATION ##########
########################################

### load the environment with the analyses run
load("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/rdata/analysis_v2.RData")
require(SNPassoc)
require(genetics)
require(MuMIn)


### select only those individuals with PA data

#save myData_ptpn1 in other object
myData_ptpn1_no_subset = myData_ptpn1

#make the subset and overwrite myData_ptpn1 with the subset
myData_ptpn1 = myData_ptpn1[which(!is.na(myData_ptpn1$PA_factor)),]

#check
identical(myData_ptpn1, myData_ptpn1_no_subset[which(!is.na(myData_ptpn1_no_subset$MVPA_mean)),])


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




#######################################################################
############ TEST DIFFERENCES SEX, AGE, BMI BETWEEN GROUPS ############
#######################################################################

##NOT DONE


##we are going to compare the 600 versus the 1057.

#check the groups
str(myData_ptpn1_no_subset) #full dataset
str(myData_ptpn1) #only PA. myData_ptpn1 is used for the PA factor subset because it is easier to recycle scripts of the tables of the main text. 

#some people say you have to compare the subset with the excluded individuals (300 in our case), but I am not really sure about that, the sample size is small? Importantly, the reviewer asked whether the 600 individuals are representative of the whole cohort: "It is important when working with a subset of the data in this way that these individuals are representative of the whole sample ". We want to know if the small group is similar the rest of the cohort.

#IMPORTANT: The distribution of age, sex and bmi are similar between groups. this is important for wilcoxon test, BUT the problem is that the groups are NOT independent!! This is fatal. We cannot compare the groups. See below.
	#https://stats.stackexchange.com/questions/181846/is-a-t-test-between-a-dataset-and-its-subset-meaningful

#compare the distribution of the variables between groups
par(mfcol=c(2,3))
plot(density(myData_ptpn1$CRF_BMI), main="BMI subset PA")
plot(density(myData_ptpn1_no_subset$CRF_BMI), main="BMI whole cohort")
plot(density(myData_ptpn1$CRF_age), main="Age subset PA")
plot(density(myData_ptpn1_no_subset$CRF_age), main="Age whole cohort")
plot(myData_ptpn1$CRF_sex, main="Sex subset PA")
plot(myData_ptpn1_no_subset$CRF_sex, main="Sex whole cohort")
	#the distributions are similar and for Age, seems normal.


##test selection

#The Wilcoxon Rank Sum test (aka Mann-Whitney) works with unequal sample sizes. The original paper (referenced below) did some analyses with DIFFERENT SAMPLE SIZES and showed its consistency and asymptotic normality (see table I, n = 8 on page 54). They also go on to show its robustness for small sample sizes. Reference: Mann, Henry B.; Whitney, Donald R. (1947). "On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other". Annals of Mathematical Statistics 18 (1): 50–60.
    #https://stats.stackexchange.com/questions/368881/wilcoxon-rank-sum-test-unequal-sample-sizes

#The Wilcoxon-Mann-Whitney test is widely used in all disciplines, probably nearly as much as the ubiquitous t-test.  Despite its lower power, it is often favoured over the t-test because of the misconception that no assumptions have to be met for the test to be valid. In fact the basic assumptions of the two tests (namely that BOTH SAMPLES ARE RANDOM SAMPLES AND ARE MUTUALLY INDEPENDENT) are identical. Not surprisingly, therefore, we find similar misuse as with the t-test concerning these aspects. But the biggest problems come where the assumptions do differ from those of the t-test - namely the distribution of the data.
    #At least the first assumption is met, as the two groups are randomly selected.
    #I think, MW test should be used only for medians? not sure...

#If the test is to be used to compare arithmetic means,  the two distributions must be BOTH SYMMETRICAL AND IDENTICAL apart from location. If the test is to be used to compare medians, the two distributions must be considered as identical apart from location. Yet we give a number of examples where distributions were clearly different, yet a significant result was still assumed to indicate either a difference in means or a difference in medians. In such situations the test is still valid to test for dominance of one distribution over another - but few researchers seem to be aware that that is what is being tested. If it is actually MEDIANS AND/OR DISTRIBUTIONS that are being compared, then reporting mean and standard error is clearly inappropriate for this test. Although most medical researchers are now aware of this, the practice is still widespread in other disciplines.
    #Both distributions should be similar except in the position.
    #In our case, they have a similar distribution more or less. See previous plots. 
        #https://influentialpoints.com/Training/Wilcoxon-Mann-Whitney_U_test_use_and_misuse.htm
        #https://influentialpoints.com/Training/Two-sample_t-test_use_and_misuse.htm

#It seems that the t test has a little bit more of assumptions than mann-whitney, although both have assumptions. I think we are fine with the assumptions of the MW test, see above.

#for sex we can use KW: Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups. It’s recommended when the assumptions of one-way ANOVA test are not met. This tutorial describes how to compute Kruskal-Wallis test in R software.
    #http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r




#################################################
############ SUPPLE TABLE S1, S2, S3 ############
#################################################

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
table_1_raw = cbind.data.frame(MAF, HWE)
rownames(table_1_raw) <- row.names(summary_snps)

#check that the HWE and MAF correspond to the same snps
table_1_raw$MAF == 1-(summary_snps$major.allele.freq/100)
table_1_raw$HWE == summary_snps$HWE

#reorder columns
table_1_raw = table_1_raw[,c("HWE", "MAF")]


### add alleles and MAF per center to table 1

#vector with centers
vector_centers = unique(myData_ptpn1$center)

#create an empty data.frame for the MAF of each center
empty_df_centers = data.frame(matrix(NA, nrow(table_1_raw), length(vector_centers)))
#set colnames
colnames(empty_df_centers) <- vector_centers

#create an empty data.frame for the major and minor alleles
empty_df_alleles = data.frame(matrix(NA, nrow(table_1_raw), 2))
#set colnames
colnames(empty_df_alleles) <- c("Major allele", "Minor allele")

#bind all DFs
table_1 = cbind.data.frame(empty_df_alleles, table_1_raw, empty_df_centers)

#add the alleles and the MAF per center
for(i in 1:nrow(table_1)){ #for each row (SNP) of table 1

    #select the [i] row
    selected_row = table_1[i,]

    #select the snp of the [i] row
    selected_snp = row.names(selected_row)

    #extract the helena and ncbi snps from [i] row
    helena_name = alleles[which(alleles$snp == selected_snp),]$helena
    ncbi_name = alleles[which(alleles$snp == selected_snp),]$ncbi

    #list alleles from helena and ncbi
    list_allels_helena = strsplit(as.vector(helena_name), split="/")[[1]]
    list_allels_ncbi = strsplit(as.vector(ncbi_name), split="/")[[1]]

    #select major and minor allele (they are ordered, major always first) for helena
    helena_major = list_allels_helena[1]
    helena_minor = list_allels_helena[2]
        #IMPORTANT: these two lines assumes that the first allele is the major and the second is the minor. I had to create the "alleles" table in that way. BE SURE.

    #select major and minor allele (they are ordered, major always first) for helena
    ncbi_major = list_allels_ncbi[1]
    ncbi_minor = list_allels_ncbi[2]
        #IMPORTANT: these two lines assumes that the first allele is the major and the second is the minor. I had to create the "alleles" table in that way. BE SURE.

    #save them in the corresponding columns
    table_1[i, which(colnames(table_1) == "Major allele")] <- ncbi_major
    table_1[i, which(colnames(table_1) == "Minor allele")] <- ncbi_minor

    #add MAF for each center
    for(j in 1:length(vector_centers)){

        #select the [j] center
        selected_center = vector_centers[j]

        #print the name of the center
        print("###############################")
        print(paste("STARTING CENTER: ", selected_center, sep=""))
        print("###############################")

        #subset myData_ptpn1 by center
        subset_center = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$center == '", selected_center, "'),]", sep="")))
        #check
        print("Check that we select rows of the selected center")
        print(summary(subset_center$center == selected_center))

        #extract allele frequencies of the [i] SNP for the [j] center
        allele_freqs = summary(subset_center[, which(colnames(subset_center) == selected_snp)])$allele.freq

        #remove the cases with NA
        allele_freqs = allele_freqs[which(!is.na(allele_freqs[,2])),]

        #extract the allele names (numbers)
        alleles_helena_numeric = row.names(allele_freqs)

        #copy the allele numbers into a new object
        alleles_helena = alleles_helena_numeric

        #change the numbers to alleles following the notation of helena
        #always A is the first, because number are afabetically ordered in Helena
        if("A" %in% list_allels_helena){
            alleles_helena = gsub("1", "A", alleles_helena)
            alleles_helena = gsub("2", list_allels_helena[!"A" == list_allels_helena], alleles_helena)
        } else { #if not, C first
            if("C" %in% list_allels_helena){ 
                alleles_helena = gsub("1", "C", alleles_helena)
                alleles_helena = gsub("2", list_allels_helena[!"C" == list_allels_helena], alleles_helena)
            } else { #if not G first (T always will be the last)
                if("G" %in% list_allels_helena){
                    alleles_helena = gsub("1", "G", alleles_helena)
                    alleles_helena = gsub("2", list_allels_helena[!"G" == list_allels_helena], alleles_helena)               
                }
            }            
        }
        
        #if the helena and ncbi name are not the same
        if(helena_name != ncbi_name){

            #copy the allele helena names in a new object
            alleles_helena_final = alleles_helena

            #convert the minor allele according to helena notation into "X"
            alleles_helena_final = gsub(helena_minor, "X", alleles_helena_final)
                #in this way, we avoid problems in palindromic SNPs. For example, if the SNP is T/A according to HELENA notation, but in ncbi is the inverse (major is minor (A/T)), then you will have problems. If you change A to T, then you have T/T, when you ask for changing T to A, the result will be A/A. In contrast, if you change the minor according to helena (A) to X (T/X), then you can change the major according to helena (T) to A, leaving A/X. Now, you can change X to T, so you have A/T. In this way, you have change from T/A to A/T. 

            #change the helena major to the correct major allele, which is indicated by ncbi.
            alleles_helena_final = gsub(helena_major, ncbi_major, alleles_helena_final)

            #change now the X, which was the helena minor, to the true minor allele according to ncbi
            alleles_helena_final = gsub("X", ncbi_minor, alleles_helena_final)
        } else { #if not, and hence the alleles are exactly the same

            #nothing to do
            alleles_helena_final = alleles_helena
        }

        #save in one dataframe the different versions of alleles names
        conversion_allele_table = cbind.data.frame(alleles_helena_numeric, alleles_helena, alleles_helena_final)
            #here the order does not mean major or minor, that was considered when the major and minor alleles according to ncbi and helena were obtained at the start of the loop. There, we used for that "alleles" object, where the first allele is the minor and the second is the major in each case.

        #select the row of allele_freqs whose name is "1", its new name will be extracted from "conversion_allele_table", selecting the final allele name of the allele "1"
        row.names(allele_freqs)[which(row.names(allele_freqs) == "1")] <- as.vector(conversion_allele_table[which(conversion_allele_table$alleles_helena_numeric == "1"),]$alleles_helena_final)

        #select the row of allele_freqs whose name is "2", its new name will be extracted from "conversion_allele_table", selecting the final allele name of the allele "2"
        row.names(allele_freqs)[which(row.names(allele_freqs) == "2")] <- as.vector(conversion_allele_table[which(conversion_allele_table$alleles_helena_numeric == "2"),]$alleles_helena_final)

        #select the frequency (second column) of the minor allele
        allele_freq_minor = allele_freqs[which(row.names(allele_freqs) == ncbi_minor), 2]

        #convert the MAF into frequency up to 1 and round to 2 decimals
        allele_freq_minor = round(allele_freq_minor/100, 2)

        #select the row of the [i] SNP and the column of the [j] center and add the corresponding MAF
        table_1[which(rownames(table_1) == selected_snp), which(colnames(table_1) == selected_center)] <- allele_freq_minor
    }
}

#check calculating the MAF with a alternative way
#this way automatically selects the allele with lower frequency in each subset by center, so if the minor/major allele changes a ina subset, the result should be different respect to table 1
for(i in 1:length(vector_centers)){

    #select the [j] center
    selected_center = vector_centers[i]

    #print the name of the center
    print("###############################")
    print(paste("STARTING CENTER: ", selected_center, sep=""))
    print("###############################")

    #subset myData_ptpn1 by center
    subset_center = eval(parse(text=paste("myData_ptpn1[which(myData_ptpn1$center == '", selected_center, "'),]", sep="")))
    #check
    print("Check that we select rows of the selected center")
    print(summary(subset_center$center == selected_center))

    #calculate summary of all SNPs for the subset of the [i] center
    summary_center = summary(subset_center)

    #calculate the new MAF
    new_MAF_center = 1-(summary_center[, which(names(summary_center) == "major.allele.freq")]/100)
    #round
    new_MAF_center = round(new_MAF_center, 2)
    #set names
    names(new_MAF_center) <- row.names(summary_center)
    #convert to DF
    new_MAF_center = data.frame(new_MAF_center)

    #merge the new MAF and table 1 by rows (SNPs)
    merged_table = merge(table_1, new_MAF_center, by="row.names")

    #check that the new MAFs is equal to the MAFs in the corresponding center
    #extract the MAF previously calculated for the [i] center
    previous_MAF = eval(parse(text=paste("merged_table[,which(colnames(merged_table) == '", selected_center, "')]", sep="")))
    print(paste("Checking differences in MAF between approaches: ", selected_center, sep=""))
    print(abs(merged_table$new_MAF_center - previous_MAF) %in% c("0", "0.01"))
        #the difference has to be zero or 0.01.
        #It a difference of 0.01 is possible due rounding with different numbers of decimals between approaches
            #round(0.445, 2) gives 0.44, while round(0.4454, 2) gives 0.45, this is the case of Zaragoza for rs2143511, for example.
}

#results of the check
    #for PTPN1, Dortmund, Heraklion and Vienna show differences in the MAF for the rs968701.
    #This SNP has a frequency very similar between alleles (almost 0.5), so it is possible that in one of the subsets by center, the minor becomes major. 
    #Indeed, I have compared the new version of table 1 with the previous version. In that version, the minor in each center was selected by frequency, independently if the minor across the whole cohort was similar to the minor in the subset. 
        #The frequency of Dortmund for the MAF of rs968701 in the previous version is 0.48, while in the new version is 0.52. 
        #The frequency of Heraklion for the MAF of rs968701 in the previous version is 0.48, while in the new version is 0.52. 
        #The frequency of Vienna for the MAF of rs968701 in the previous version is 0.45, while in the new version is 0.55.
    #It seems that the calculations are ok, just the minor allele for rs968701 has become the major in these three subsets. Because of this, in the new version we have these cases with the minor having a frequency higher than 0.5.

#reorder the table following the order in the chromosome 
table_1 = table_1[match(ptpn1_snps$snp, row.names(table_1)),]
row.names(table_1) == ptpn1_snps$snp


### check that the alleles names are ok
#extract the allele names from the table 1 (combined major and minor) along with snp names
alleles_from_table_1 = cbind.data.frame(row.names(table_1), paste(table_1[,which(colnames(table_1) == "Major allele")], "/", table_1[,which(colnames(table_1) == "Minor allele")], sep=""))
colnames(alleles_from_table_1) <- c("snp", "alleles_from_ncbi")

#merge these alleles names with the original allele names from alleles df
df_check_alleles = merge(alleles_from_table_1, alleles)

#check that the colums of alleles names are similar
identical(as.vector(df_check_alleles$alleles_from_ncbi), as.vector(df_check_alleles$ncbi))


### add the gene name to the table. This makes more sense when SNPs of multiple genes are included.
#loop for each gene
ptpn1s = unique(ptpn1_snps$gen)
#for each gene name
for(i in 1:length(ptpn1s)){

    #select the [i] gene
    ptpn1_selected = ptpn1s[i]

    #select the position of the first snp of the [i] gene
    position_to_add = which(row.names(table_1) == ptpn1_snps[which(ptpn1_snps$gen==ptpn1_selected),][1,]$snp)

    #create the row to be added
    row_to_add = data.frame(matrix(NA, 1, ncol(table_1)))
    colnames(row_to_add) <- colnames(table_1) #add the same names that in the data.frame

    #convert to upper case
    row.names(row_to_add) <- toupper(ptpn1_selected)

    #add the gene name
    if(position_to_add == 1){# if the position is the 1

        #bind first the gene name
        table_1 = rbind(row_to_add, table_1)
    } else{ #if the gene name have to be include in the middle of the table
        
        #add the gene name in the middle of the table, after the prior gene and before the snps of this gene [i]        
        table_1 = rbind(table_1[1:(position_to_add-1),], row_to_add, table_1[(position_to_add):nrow(table_1),])        
    }
}


### add the MAF according to 1000 Genomes Project

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




################################################
##### CONVERT TABLES TO LATEX AND COMPILE ######
################################################

#load the required package
require(xtable)

#path to save tex table
path_tex_table = "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/helena_study/helena_7/results/tables/"

#name tex table
name_tex_table = "tables_latex_sampling_v1.tex"

#name doc table
name_doc_table = "tables_latex_sampling_v1.odt"

#remove the previous file with tables in latex
system(paste("cd ", path_tex_table, "; rm ", name_tex_table, "; rm ", name_doc_table, sep=""))

#convert table 1 to a latex table
print.xtable(xtable(table_1_1KGP, caption="Supple Table S1", label=NULL, align="llccccccccccccccc", digits=2, display=c("s", "s", "s", "s", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_1_1KGP)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date())
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
print.xtable(xtable(table_2, caption="Supple Table S2", label=NULL, align="llcccccc", digits=c(0,0,0,2,0,2,0,2), display=c("s", "s", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_2)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #digits: We select the number of digits because we want each column to have a different number of digits. 
        #Numeric vector of length equal to one (in which case it will be replicated as necessary) or to the number of columns of the resulting table plus 1 for the column with the row.names. These row names will be then removed setting FALSE the argument "include.rownames" from print.xtable. 
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...
    #Rest of the argument in the line of the table 1

#convert table 3 to a latex table
print.xtable(xtable(table_3, caption="Supple Table S3", label=NULL, align="llcccccc", digits=2, display=c("s", "s", "f", "f", "f", "f", "f", "f")), type="latex", file=paste(path_tex_table, name_tex_table, sep=""), append=TRUE, floating=TRUE, table.placement="ht", caption.placement="top", caption.width=NULL, latex.environments="center", hline.after=c(-1,0,nrow(table_3)), NA.string="", include.rownames=FALSE, comment=TRUE, timestamp=date(), sanitize.text.function=function(x) {x})
    #sanitize.text.function:
        #All non-numeric entries (except row and column names) are sanitized in an attempt to remove characters which have special meaning for the output format. If sanitize.text.function is not NULL, it should be a function taking a character vector and returning one, and will be used for the sanitization instead of the default internal function. Default value is NULL.
            #in this way we can remove slash, etc...
    #Rest of the argument in the line of the table 1

#compile to odt with pandoc
system(paste("cd ", path_tex_table, "; pandoc -s ", name_tex_table, " -o ", name_doc_table, sep=""))

#From the ".odt" file, you have to copy each table to excel, using advanced options for pasting ("text"). In that way, you can get the +- and other symbols correctly pasted. Once you have the first version in excel with the column names correctly displayed (squared, etc...), you can just paste the content (phenotype values, frequencies...).




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

#save the previous geno_pheno_results
geno_pheno_results_no_subset = geno_pheno_results

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

#compare the R2 of the associations with the full dataset with the subset
kruskal.test(geno_pheno_results$d2_mumin, geno_pheno_results_no_subset$d2_mumin)