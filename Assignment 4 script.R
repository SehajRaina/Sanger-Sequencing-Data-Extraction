
library(BiocManager)
library(sangerseqR)
library(seqinr)


SeqStats<-read.csv("Data/BarcodePlateStats.csv") #quality control of sequence statistics

GoodQuality<-subset(SeqStats,Ok=="TRUE")

# get full names including folder path
filesName<-list.files("Data", pattern=".ab1", full.names = TRUE)

# then keep only the basename (file names) matching dataframe column
filesName2 <- filesName[basename(filesName) %in% GoodQuality$Chromatogram]
basename(filesName2)
# then read the data
ITS<-lapply(filesName2, read.abif) #read

ITSseq <- lapply(ITS, sangerseq) #extract

SeqX <- lapply(ITSseq, makeBaseCalls) #call


SeqX1<-list() #blank list
SeqX2<-list() #blank list


for(i in 1:length(SeqX)){
    SeqX1[i]<-SeqX[[i]]@primarySeq
    SeqX2[i]<-gsub(".*seq","",SeqX1[[i]])  #only get the primary sequence and put it in a list
}



plzWORK <- paste(" ", data2$new, data1$V1, sep="\n")

cat(plzWORK)

write.table(plzWORK, file="omg", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

df <- read.table("omg", header = FALSE)



getwd()
library(seqinr)
write.fasta(sequences=SeqX2,names=basename(filesName2),nbchar = 80, as.string = TRUE, file.out="DNA")
