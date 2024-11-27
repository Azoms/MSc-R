#setting working directory
getwd()
setwd("C:/Users/Shalom/OneDrive - Brunel University London/2209539_script")







#Question 1
#Saving pheno_covar.txt into a data frame

df<-read.table("data/pheno_covar.txt",header = TRUE, sep = " ")
head(df)      #viewing the first 6 observations


#Creating a scatter plot of age vs trait 1 
plot(df$age, df$trait1, main= "Scatterplot of Age vs Trait 1",
     xlab = "Age", ylab="Trait 1")
abline(lm(df$trait1~df$age), col = "blue")  #adding a line of best fit



#creating a plot ggplot of trait 1 per sex
library(ggplot2)  #loading package into R

ggplot(data=df, aes(x= factor(sex), y = trait1))+
  geom_bar(stat = "identity", fill ="red") +
  ggtitle("Variation of Trait1 with sex") +
  xlab("Sex") +
  ylab ("Trait1") +
  scale_x_discrete(labels = c('Females', 'Males'))








#Question 2.

#creating and saving the phenotype data frame into a file
pheno.df<-df[, 1:3]
write.table(pheno.df,"Results/pheno.txt", row.names = FALSE, quote = FALSE)

#creating a covariate a data frame
covar.df<-df[, (-3)]

#saving the co variate data frame into a file
write.table(covar.df, file= "Results/covar.txt", row.names = FALSE, quote = FALSE)








#Question 3
#running system command to test for plink

system("plink/plink",wait=TRUE)

#installing and loading packages to help view file formats into R 
library(data.table)
library(BEDMatrix)
library(qqman)

#exploring the data
fam <- read.table("data/BB5707.fam",sep=" ",header=FALSE)
print(head(fam))

bim <- read.table("data/BB5707.bim",sep="\t",header=FALSE)
print(head(bim))

bed <- BEDMatrix("data/BB5707.bed")
print(bed[1:4,1:5])

#Running a GWAS analysis
system("plink/plink --bfile data/BB5707 --pheno Results/pheno.txt  --covar Results/covar.txt --linear --ci 0.95 --adjust --out Results/Bb5707_GWAS_pheno")
 

#finding independent statistically significant SNPs
system("plink/plink --bfile data/BB5707 --clump Results/Bb5707_GWAS_pheno.assoc.linear --clump-p1 1e-05 --clump-r2 0.1 --clump-kb 300 --out Results/2209539.assoc.linear_clumped")


#read clumped file
GWASresults_clumped <- read.table("Results/2209539.assoc.linear_clumped.clumped", header=T)
GWASresults_clumped[1:9 ,1:7]    #view clumped file







#Question 4
#Creating a data frame with the file

df1<-read.table("data/differential_expression.txt",header = TRUE, sep = "")


#Stratifying the data in an added column (Significance1) based on p-values

Diff.expr<-df1      #creating a new dataframe
Diff.expr$Significance1 <- "NS"
Diff.expr<-within(Diff.expr, Significance1[significance=="Significant"]<- "Significant")
Diff.expr<-within(Diff.expr, Significance1[significance=="Borderline" ]<-"Borderline")


#generating a vulcan plot with NS-Gray, Borderline-Blue, Significant-Red

ggplot(Diff.expr,aes(x=log2FoldChange,y=pvalue,col=Significance1))+geom_point()+
  scale_color_manual(values=c("Blue","Gray","Red"))+
  labs(x="Log2Fold Change", y="P-Value")

#to give a better spread of data, plot the -log10 of Pvalues

ggplot(Diff.expr,aes(x=log2FoldChange,y=-log10(pvalue),col=Significance1))+geom_point()+
  scale_color_manual(values=c("Blue","Gray","Red"))+
  labs(x="Log2Fold Change", y="-Log10 (P-Value")







#Question 5
#Generating a basic allele score with clumped file
#reading the GWAS association file

GWASresults <- read.table("Results/Bb5707_GWAS_pheno.assoc.linear", 
                          header=T) 

#order file by P-value
GWASresults <- GWASresults[order(GWASresults$P),]

#display first 15 SNPs
GWASresults [1:15,]

#create a new data frame with independent associated SNPs
#from GWAS association result dataframe

GWASresults.df<- GWASresults[1:13,c("SNP","A1","BETA")] 
clumpedSNPs<- GWASresults.df[c(1,3,4:6,8,9,12,13),]  #new data frame


#writing the data frame into a text file 
write.table(clumpedSNPs,"Results/clumpedSNPs.txt", row.names = FALSE, quote = FALSE)


#generating a score profile with independent SNPs
system ("plink/plink --bfile data/BB5707 --score Results/clumpedSNPs.txt  --out Results/2209539_score")







#Question 6

#saving the score_profile into a data frame
SP.df<-read.table("Results/2209539_score.profile",header = T, stringsAsFactors = F)


#extracting trait1 from the data frame (df)
trait1<-df[, 3]

SP.df$trait1<-trait1 #adding trait1 column into score profile data frame


#running a correlation on allele score and trait1

fit1<-lm(SCORE~trait1, data=SP.df) # simple linear regression
summary(fit1)







