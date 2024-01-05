#Extract only K70 from Hesin_diag file. 
K70<-subset(hesin_diag,  startsWith(as.character(hesin_diag$diag_icd10), 'K70'))

            #remove duplicates
K70<- as.data.frame(K70[!duplicated(K70$eid),])

# create a new column that K70 equals 1
K70$K70<-1

#select only the ID and K70
K70<-select(K70, c(eid_1, K70))
