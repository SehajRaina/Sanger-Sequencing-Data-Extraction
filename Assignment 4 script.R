###Assignment 4 
###Sehaj Raina


#1. Load packages 

library(BiocManager)
library(sangerseqR)
library(seqinr) #for exporting data into Fasta file 


#2. Quality control of sequence statistics

SeqStats<-read.csv("Data/BarcodePlateStats.csv") #load data into object 
GoodQuality<-subset(SeqStats,Ok=="TRUE") #subset rows that are OK aka. TRUE for the quality check 

filesName<-list.files("Data", pattern=".ab1", full.names = TRUE) #get full names including folder path of all sanger sequencing files
filesName2<-filesName[basename(filesName) %in% GoodQuality$Chromatogram] #Only keep the sanger sequencing files whos basename (file names) match in GoodQuality$Chromatogram column


#3. Read, extract and call the data using lapply function to apply function to all elements of a list 

ITS<-lapply(filesName2, read.abif) #read
ITSseq<-lapply(ITS, sangerseq) #extract
SeqX<-lapply(ITSseq, makeBaseCalls) #call


#4. Extract Primary Sequence and place in a list 

SeqX1<-list() #blank list
SeqX2<-list() #blank list


for(i in 1:length(SeqX)){
    SeqX1[i]<-SeqX[[i]]@primarySeq #only get the primary sequence and put it in SeqX1 list
    SeqX2[i]<-gsub(".*seq","",SeqX1[[i]]) #use regular expression on SeqX1 to remove everything before the sequence itself and put it into new SeqX2 list
}


#5. Export primary sequence data to a FASTA file 

DNAlength<-paste("length=",nchar(SeqX2),";", sep="") #create object with length of each sequence, FASTA style 
gene<-paste(basename(filesName2), DNAlength, "type=DNA", sep=" ") #create object stringing each file name to it's respective sequence length from DNAlength, FASTA style
  
write.fasta(sequences=SeqX2,names=gene,nbchar = 80, as.string = TRUE, file.out="DNA") #assign names from gene object to each respective sequence and output into FASTA format file
