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