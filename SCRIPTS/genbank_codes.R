#install.packages("ape")
#install.packages("seqinr")
library(ape)
library(seqinr)

##trnl_trnf
setwd("D:/Unifesp/Mestrado/analises/art_especies/dantas2")
list.files()

#apos criar a tabela .xls com as especies, salvar como formato .csv

data1 <- read.table(file = "dantas_trnl_trnf.csv", header=T, sep = ";", dec = ".")

linhas <- nrow(data1)

data <- data.frame()

for (i in 1:linhas){
  y <- data1[i,1]
  x <- read.GenBank(y)
  aa <- data1[i,2]
  write.dna(x, file = file.path(paste0('d.fasta')), format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)
  xx <- read.fasta(file = "d.fasta", seqtype = "DNA", as.string=T, forceDNAtolower = F)
  z <- paste(attr(x, "description"), names(x)) 
  data[i,1] <- attr(x, "species")
  data[i,2] <- attr(x, "description")
  write.fasta(sequences = xx, names = z, file.out = file.path(paste0(aa,"_",y,'.fasta')))
  }

write.table(data, file = "dados_seq.txt", quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ",", row.names = TRUE, qmethod = c("escape", "double"))


##ycf1
setwd("D:/Unifesp/Mestrado/analises/art_especies/dantas2/ycf1")
list.files()

#apos criar a tabela .xls com as especies, salvar como formato .csv

data1 <- read.table(file = "dantas_ycf1.csv", header=T, sep = ";", dec = ".")

linhas <- nrow(data1)

data <- data.frame()

for (i in 1:linhas){
  y <- data1[i,1]
  x <- read.GenBank(y)
  aa <- data1[i,2]
  write.dna(x, file = file.path(paste0('d.fasta')), format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)
  xx <- read.fasta(file = "d.fasta", seqtype = "DNA", as.string=T, forceDNAtolower = F)
  z <- paste(attr(x, "description"), names(x)) 
  data[i,1] <- attr(x, "species")
  data[i,2] <- attr(x, "description")
  write.fasta(sequences = xx, names = z, file.out = file.path(paste0(aa,"_",y,'.fasta')))
}

write.table(data, file = "dados_seq.txt", quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ",", row.names = TRUE, qmethod = c("escape", "double"))
