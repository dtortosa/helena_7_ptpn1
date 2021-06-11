myData_lipg$PA_included_excluded <- NA

myData_lipg[which(!is.na(myData_lipg$PA_factor)),]$PA_included_excluded <- "included"
myData_lipg[which(is.na(myData_lipg$PA_factor)),]$PA_included_excluded <- "excluded"



wilcox.test(myData_lipg[which(myData_lipg$PA_included_excluded == "excluded"),]$CRF_age, myData_lipg[which(myData_lipg$PA_included_excluded == "included"),]$CRF_age)

wilcox.test(myData_lipg[which(myData_lipg$PA_included_excluded %in% c("excluded", "included")),]$CRF_age, myData_lipg[which(myData_lipg$PA_included_excluded == "included"),]$CRF_age)



wilcox.test(myData_lipg[which(myData_lipg$PA_included_excluded == "excluded"),]$CRF_BMI, myData_lipg[which(myData_lipg$PA_included_excluded == "included"),]$CRF_BMI)

wilcox.test(myData_lipg[which(myData_lipg$PA_included_excluded %in% c("excluded", "included")),]$CRF_BMI, myData_lipg[which(myData_lipg$PA_included_excluded == "included"),]$CRF_BMI)



chisq.test(table(myData_lipg$PA_included_excluded, myData_lipg$CRF_sex))

kruskal.test(CRF_age ~ PA_included_excluded, data = myData_lipg)
kruskal.test(CRF_BMI ~ PA_included_excluded, data = myData_lipg)
kruskal.test(CRF_sex ~ PA_included_excluded, data = myData_lipg)

SI COMPARAMOS LOS 1045 CON LOS 698 NO HAY DIFERENCIAS, PERO SI COMPARAMOS LOS 698 CON LOS 300 RESTANTES SI LAS HAY. SI COPARAMOS LOS 698 CON EL RESTO DE LOS 3000 TAMBIÉN HAY DIFERENCIAS. HAY QUE PENSAR SOBRE ESTO. 

CON LO PRIMEOR SERIA SUFICIENTE, EL PROBLEMA ES QUE MANWITNE COSNIDERA QUE AMBAS MUESTRAS SON INDEPENDEIENTES... Y REALMENTE NO LO SON... NO PUEDO USAR KW AL NO PODER CREAR UN FACTOR. BUSCAR PRUEBA PARA MUESTRAS NO INDEPENDEIENTES?  mira codifo tablas lipg

realmente lo que nos interesa saber es si el total es igual que ese subsample, así que no debería estar mal. Mira a ver si lo de man-witnee es grave, porque la asuncion más importnate, de normalidad, creo que la cumplimos (mira lipg code tables). En el caso de los 3000 y los 1047 son un tercio, y fueron seleccionados al azar. Es al azar y la muestra es grande, aquí tenemos una muestra más pequeña y no es al azar, si no los que tienen datos....



#####################################################################
##### COMPARISON HELENA 1057 AND SELECTED SUBSET FOR INTERACTION ####
#####################################################################

#select those rows with data for PA factor
myData_ptpn1_pa_subset = myData_ptpn1[which(!is.na(myData_ptpn1$PA_factor)),]
#check
identical(myData_ptpn1_pa_subset, myData_ptpn1[which(!is.na(myData_ptpn1$MVPA_mean)),])

#the whole dataset with genetic data is "myData_ptpn1", while is the subset used in the analyses. "myData" only has 1057 individuals, similar to "myData_ptpn1", while "helena" has 3552 rows. 
str(myData_ptpn1)
str(myData_ptpn1_pa_subset)


##we are going to create a new variable in "helena" to differentiate between those individuals included and excluded for the final dataset.

#the reviewer asked for comparing between the included and the excluded: Are there significant differences found between those who remained for statistical analysis and those who were eventually excluded?

#we are going to compare age, sex and BMI between excluded and included. 

#open a new variable
helena$included_excluded <- NA

#those individuals that were removed from helena to create myData_ptpn1 are considered as excluded. They were remove using the variable "rows_to_remove_lipg"
helena[rows_to_remove_lipg,]$included_excluded <- "excluded"

#the rest of individuals are those included in myData_ptpn1, so we selected all rows except those included in "rows_to_remove_lipg"
helena[-rows_to_remove_lipg,]$included_excluded <- "included"

#convert the variable into a factor 
helena$included_excluded <- as.factor(helena$included_excluded)
summary(helena$included_excluded) #2495 excluded and 1057 included.

#check that three variables are identical in the included and myData_ptpn1, which is the dataset resulting from excluding 2 thirds of individuals without genotyping data
identical(helena[which(helena$included_excluded == "included"),]$CRF_BMI, myData_ptpn1$CRF_BMI)
identical(helena[which(helena$included_excluded == "included"),]$CRF_sex, myData_ptpn1$CRF_sex)
identical(helena[which(helena$included_excluded == "included"),]$CRF_age, myData_ptpn1$CRF_age)

#the same test but now selection the included individuals by removing those excluded in the final analysis
identical(helena[-which(helena$included_excluded == "excluded"),]$CRF_BMI, myData_ptpn1$CRF_BMI)
identical(helena[-which(helena$included_excluded == "excluded"),]$CRF_sex, myData_ptpn1$CRF_sex)
identical(helena[-which(helena$included_excluded == "excluded"),]$CRF_age, myData_ptpn1$CRF_age)

#compare the distribution of the variables between excluded and included
par(mfcol=c(2,3))
plot(density(helena[which(helena$included_excluded == "included"),]$CRF_BMI))
plot(density(helena[which(helena$included_excluded == "excluded"),]$CRF_BMI))
plot(density(helena[which(helena$included_excluded == "included"),]$CRF_age))
plot(density(helena[which(helena$included_excluded == "excluded"),]$CRF_age))
plot(helena[which(helena$included_excluded == "included"),]$CRF_sex)
plot(helena[which(helena$included_excluded == "excluded"),]$CRF_sex)


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

##modeling

#with GLM
model_sex = glm(helena$included_excluded ~ helena$CRF_sex, family="binomial")
model_age = glm(helena$included_excluded ~ helena$CRF_age, family="binomial")
model_bmi = glm(helena$included_excluded ~ helena$CRF_BMI, family="binomial")
    #ESTO NO VALE, HAABRÍA QUE MIRAR BIEN RESIDUALES DEL BINOMIAL
        #http://www.sthda.com/english/articles/36-classification-methods-essentials/148-logistic-regression-assumptions-and-diagnostics-in-r/


#Mann-witney testing if any of the distributions has a higher value than the other
#age
paste("Median Age of whole HELENA cohort", median(helena$CRF_age))
paste("Median Age of subset", median(myData_ptpn1$CRF_age))
wilcox.test(helena[which(helena$included_excluded == "excluded"),]$CRF_age, helena[which(helena$included_excluded == "included"),]$CRF_age)
#BMI
paste("Median BMI of whole HELENA cohort", median(helena$CRF_BMI))
paste("Median BMI of subset", median(myData_ptpn1$CRF_BMI))
wilcox.test(helena[which(helena$included_excluded == "excluded"),]$CRF_BMI, helena[which(helena$included_excluded == "included"),]$CRF_BMI)
#sex
#cannot be analyzed because it is a factor. We can calculate a contingency table and the apply a chi-square text, but then we have to met the assumptions
chisq.test(table(helena$included_excluded, helena$CRF_sex))
    #https://www.datacamp.com/community/tutorials/contingency-tables-r


#Kruskal-Wallis
#The Kruskal-Wallis non-parametric test is used for comparing ordinal or non-Normal variables for more than two groups, and is a generalisation of the Mann-Whitney U test
kruskal.test(CRF_age ~ included_excluded, data = helena)
kruskal.test(CRF_BMI ~ included_excluded, data = helena)
kruskal.test(CRF_sex ~ included_excluded, data = helena)
    #p-values of BMI and age are very similar to those obtained with the MW test, which makes sense. We will use these results. 
