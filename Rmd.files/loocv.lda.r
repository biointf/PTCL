loocv.lda = function(data) {

library(MASS)
y = t(data)
cl.orig = factor(a)

y.red <- y [ cl.orig=="ALK-" |  cl.orig=="PTCL" , ]
cl.orig.red = as.matrix ( cl.orig [ cl.orig=="ALK-" |  cl.orig=="PTCL" ] )

#### with ALK+
#y.red <- y [ cl.orig=="ALK-" |  cl.orig=="PTCL" |  cl.orig=="ALK+", ]
#cl.orig.red = as.matrix ( cl.orig [ cl.orig=="ALK-" |  cl.orig=="PTCL" |  cl.orig=="ALK+" ] )
#cl.orig.red [ cl.orig.red == "ALK-" | cl.orig=="ALK+" ] = "ALCL"

perm.mother = rownames(y.red)
perm.son = combn (perm.mother, length(perm.mother)-1) 

# x.lda = x.red [, c("CD30.dCt" , "GZMB.dCt" , "SNFT.dCt" , "TMOD.dCt" ) ]
x.lda = y.red [, c("TNFRSF8" , "SNFTnew" , "TMOD" ) ]

output <- cbind(perm.mother, NA, NA, NA, NA)

for (i in 1:length(perm.mother)) {
  train <- x.lda [ perm.son[,i], ] 
  test <- x.lda [ ! ( rownames(x.lda) %in% perm.son[,i]) , ]
  cl <- cl.orig.red [which(rownames(x.lda)%in%perm.son[,i])]
  z <- lda(train, cl)
  p <- predict(z,test)$class
  output  [ output[,1] == rownames(test) , 2  ] = as.character(p)
  output  [ output[,1] == rownames(test) , 3  ] = z$scaling [1,1]
  output  [ output[,1] == rownames(test) , 4  ] = z$scaling [2,1]
  output  [ output[,1] == rownames(test) , 5  ] = z$scaling [3,1]
  
}

colnames(output) = c("true","LOOCV.predicted", "CD30", "BATF3","TMOD")
output

}