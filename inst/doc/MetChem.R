## ---- echo=FALSE, results='asis'----------------------------------------------
Feature=c(
 "Cluster1, median [IQR]", "Cluster2, median [IQR]",  "Cluster3, median [IQR]",  "Cluster4, median [IQR]" , "Cluster5, median [IQR]" , "Cluster6, median [IQR]", 
 "Cluster7, median [IQR]", "Cluster8, median [IQR]")

MYC=c( "4.15 [2.003 6.629]","-2.57 [-4.056 -1.642]","3.482 [1.47 5.657]","-3.325 [-5.853 -1.232]","3.438 [2.472 4.082]","3.973 [3.794 4.45]","3.256 [2.382 4.548]","-2.279 [-2.99 -1.479]")

WT=c( "-4.379 [-6.376 -2.853]","2.423 [1.833 3.68]","-3.779 [-4.305 -2.423]","2.934 [2.356 3.926]","-3.376 [-3.839 -3.064]","-4.192 [-4.657 -3.695]","-3.406 [-3.94 -2.726]","2.413 [1.805 3.118]" )

pvalue=c( "3.66e-05","6.01e-05","7.66e-05","6.1e-03","3.66e-05","3.66e-05","3.66e-05","4.69e-05")

FDR=c(  "7.32e-05","8.01e-05","8.75e-05","6.10e-03","7.32e-05","7.32e-05","7.32e-05","7.51e-05")

da=data.frame(Feature=Feature,MYC=MYC,WT=WT,'p-value'=pvalue,FDR=FDR,check.names = FALSE)
  
knitr::kable(da, align = "lcccc")

## ---- echo=FALSE, results='asis'----------------------------------------------
pval=c( 5.56e-11,1.74e-09,1.81e-08,2.94e-07,2.40e-06,9.82e-06,2.20e-05,6.85e-05,1.05e-04,2.26e-04,2.26e-04,2.26e-04)
nam=c( "Heteroaromatic compound","Azacycle",
       "Aromatic heteromonocyclic compound","Imidazolyl carboxylic acid derivative",
       "Aralkylamine","Benzenoid",
       "Aromatic homomonocyclic compound","1-hydroxy-2-unsubstituted benzenoid",
       "Azole","Indole",
       "Pyrrole","Substituted pyrrole"  )
da=data.frame(Substituents=nam,'p-value'=pval,check.names = FALSE)
knitr::kable(da, align = "lc")

