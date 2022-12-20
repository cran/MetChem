## ---- echo=FALSE, results='asis'----------------------------------------------
Feature=c(
 "Cluster1, median [IQR]", "Cluster2, median [IQR]",  "Cluster3, median [IQR]",  "Cluster4, median [IQR]" , "Cluster5, median [IQR]" , "Cluster6, median [IQR]", 
 "Cluster7, median [IQR]")

MYC=c(  "4.15 [2.003 6.629]"  ,"-2.57 [-4.056 -1.642]","3.82 [1.703 5.572]", "-3.325 [-5.853 -1.232]","3.54 [2.873 4.312]",     "3.708 [3.185 4.399]"  ,  "4.353 [3.232 5.841]" )

WT=c("-4.379 [-6.376 -2.853]","2.423 [1.833 3.68]","-3.801 [-4.361 -2.586]","2.934 [2.356 3.926]" , "-3.64 [-4.329 -3.347]" , "-3.74 [-4.569 -3.505]","-4.211 [-4.687 -3.718]")

pvalue=c(  "3.66e-05","6.01e-05","7.66e-05","6.1e-03" ,"3.66e-05","3.66e-05","3.66e-05")

FDR=c( "6.40e-05","8.41e-05","8.93e-05","6.10e-03","6.40e-05","6.40e-05","6.40e-05")

da=data.frame(Feature=Feature,MYC=MYC,WT=WT,'p-value'=pvalue,FDR=FDR,check.names = FALSE)
  
knitr::kable(da, align = "lcccc")

## ---- echo=FALSE, results='asis'----------------------------------------------
pval=c(2.41e-13,1.48e-11,1.40e-10,9.02e-09,7.67e-08,1.28e-05,2.53e-05,4.14e-05,4.79e-05,1.16e-04,3.24e-04,3.24e-04,3.24e-04,7.27e-04)
nam=c( "Heteroaromatic compound","Azacycle",
  "Aromatic heteromonocyclic compound","Imidazolyl carboxylic acid derivative",
 "Aralkylamine","Azole",
  "Benzenoid","Histidine or derivatives",
  "Aromatic homomonocyclic compound","1-hydroxy-2-unsubstituted benzenoid",
 "Indole","Pyrrole",
 "Substituted pyrrole","Imidazole"  )
da=data.frame(Substituents=nam,'p-value'=pval,check.names = FALSE)
knitr::kable(da, align = "lc")

