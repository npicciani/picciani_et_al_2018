library(corHMM)
library(phytools)

############### opsin mapping ################
#############################################

#reading the trees and the trait file
tree<-read.tree("rep4_cnidops.treefile")
opsinexp<-read.csv("rep4_eyes.csv",head=TRUE,sep=',')


#estimating the ancestral states
## running in an orderly way to later choose the model

recon_ARD <- rayDISC(tree,opsinexp,charnum=1,model="ARD",node.states="marginal")
recon_ER <- rayDISC(tree,opsinexp,charnum=1,model="ER",node.states="marginal")
recon_SYM <- rayDISC(tree,opsinexp,charnum=1,model="SYM",node.states="marginal")

write.csv(recon_ARD$states, file="statesARD-opsinexp.csv" )
write.csv(recon_ER$states, file="statesER-opsinexp.csv" )
write.csv(recon_SYM$states, file="statesSYM-opsinexp.csv" )
