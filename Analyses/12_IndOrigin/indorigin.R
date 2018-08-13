library(indorigin)
library(ape)

###reading the ufboot trees

tree<-read.tree(file="iqtree_selected_no_long_unpar_rep10_boot.ufboot",keep.multi = TRUE)
tree.num=length(tree)

###rooting then dropping the outgroups
rootedtrees <- root.multiPhylo(phy=tree,outgroup=c("Amphimedon_queenslandica","Trichoplax_adhaerens","Crassostrea_gigas","Strongylocentrotus_purpuratus"),resolve.root = TRUE)
pruned.trees<-list();class(pruned.trees)<-"multiPhylo"
pruned.trees<-lapply(rootedtrees,drop.tip, tip=c("Amphimedon_queenslandica","Trichoplax_adhaerens","Crassostrea_gigas","Strongylocentrotus_purpuratus"))
class(pruned.trees)<-"multiPhylo"
eyes<-read.csv("eyes.csv",header=TRUE)
tip.num=dim(eyes)[1]
eyesVec=numeric(tip.num)
tipNames=as.character(eyes[,1])
eyesVec=eyes[,2]
names(eyesVec)=tipNames

#######******** rate gain:loss=1:100

########1gains
Indorigin1gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=1)
bf1gains_1_100<-getBF(Indorigin1gains_1_100)
posterior1gains_1_100<- getPostProb(Indorigin1gains_1_100)
prior1gains_1_100<-getPriorProb(Indorigin1gains_1_100)
bf1gains_1_100
posterior1gains_1_100
prior1gains_1_100

########2gains
Indorigin2gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=2)
bf2gains_1_100<-getBF(Indorigin2gains_1_100)
posterior2gains_1_100<- getPostProb(Indorigin2gains_1_100)
prior2gains_1_100<-getPriorProb(Indorigin2gains_1_100)
bf2gains_1_100
posterior2gains_1_100
prior2gains_1_100

########3gains
Indorigin3gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=3)
bf3gains_1_100<-getBF(Indorigin3gains_1_100)
posterior3gains_1_100<- getPostProb(Indorigin3gains_1_100)
prior3gains_1_100<-getPriorProb(Indorigin3gains_1_100)
bf3gains_1_100
posterior3gains_1_100
prior3gains_1_100

########4gains
Indorigin4gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=4)
bf4gains_1_100<-getBF(Indorigin4gains_1_100)
posterior4gains_1_100<- getPostProb(Indorigin4gains_1_100)
prior4gains_1_100<-getPriorProb(Indorigin4gains_1_100)
bf4gains_1_100
posterior4gains_1_100
prior4gains_1_100

########5gains
Indorigin5gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=5)
bf5gains_1_100<-getBF(Indorigin5gains_1_100)
posterior5gains_1_100<- getPostProb(Indorigin5gains_1_100)
prior5gains_1_100<-getPriorProb(Indorigin5gains_1_100)
bf5gains_1_100
posterior5gains_1_100
prior5gains_1_100

########6gains
Indorigin6gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=6)
bf6gains_1_100<-getBF(Indorigin6gains_1_100)
posterior6gains_1_100<- getPostProb(Indorigin6gains_1_100)
prior6gains_1_100<-getPriorProb(Indorigin6gains_1_100)
bf6gains_1_100
posterior6gains_1_100
prior6gains_1_100

########7gains
Indorigin7gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=7)
bf7gains_1_100<-getBF(Indorigin7gains_1_100)
posterior7gains_1_100<- getPostProb(Indorigin7gains_1_100)
prior7gains_1_100<-getPriorProb(Indorigin7gains_1_100)
bf7gains_1_100
posterior7gains_1_100
prior7gains_1_100

########8gains
Indorigin8gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=8)
bf8gains_1_100<-getBF(Indorigin8gains_1_100)
posterior8gains_1_100<- getPostProb(Indorigin8gains_1_100)
prior8gains_1_100<-getPriorProb(Indorigin8gains_1_100)
bf8gains_1_100
posterior8gains_1_100
prior8gains_1_100

########9gains
Indorigin9gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=9)
bf9gains_1_100<-getBF(Indorigin9gains_1_100)
posterior9gains_1_100<- getPostProb(Indorigin9gains_1_100)
prior9gains_1_100<-getPriorProb(Indorigin9gains_1_100)
bf9gains_1_100
posterior9gains_1_100
prior9gains_1_100

########10gains
Indorigin10gains_1_100=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=100, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=10)
bf10gains_1_100<-getBF(Indorigin10gains_1_100)
posterior10gains_1_100<- getPostProb(Indorigin10gains_1_100)
prior10gains_1_100<-getPriorProb(Indorigin9gains_1_100)
bf10gains_1_100
posterior10gains_1_100
prior10gains_1_100


#######******** rate gain:loss=1:10
########1gains
Indorigin1gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=1)
bf1gains_1_10<-getBF(Indorigin1gains_1_10)
posterior1gains_1_10<- getPostProb(Indorigin1gains_1_10)
prior1gains_1_10<-getPriorProb(Indorigin1gains_1_10)
bf1gains_1_10
posterior1gains_1_10
prior1gains_1_10

########2gains
Indorigin2gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=2)
bf2gains_1_10<-getBF(Indorigin2gains_1_10)
posterior2gains_1_10<- getPostProb(Indorigin2gains_1_10)
prior2gains_1_10<-getPriorProb(Indorigin2gains_1_10)
bf2gains_1_10
posterior2gains_1_10
prior2gains_1_10


########3gains
Indorigin3gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=3)
bf3gains_1_10<-getBF(Indorigin3gains_1_10)
posterior3gains_1_10<- getPostProb(Indorigin3gains_1_10)
prior3gains_1_10<-getPriorProb(Indorigin3gains_1_10)
bf3gains_1_10
posterior3gains_1_10
prior3gains_1_10

########4gains
Indorigin4gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=4)
bf4gains_1_10<-getBF(Indorigin4gains_1_10)
posterior4gains_1_10<- getPostProb(Indorigin4gains_1_10)
prior4gains_1_10<-getPriorProb(Indorigin4gains_1_10)
bf4gains_1_10
posterior4gains_1_10
prior4gains_1_10

########5gains
Indorigin5gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=5)
bf5gains_1_10<-getBF(Indorigin5gains_1_10)
posterior5gains_1_10<- getPostProb(Indorigin5gains_1_10)
prior5gains_1_10<-getPriorProb(Indorigin5gains_1_10)
bf5gains_1_10
posterior5gains_1_10
prior5gains_1_10

########6gains
Indorigin6gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=6)
bf6gains_1_10<-getBF(Indorigin6gains_1_10)
posterior6gains_1_10<- getPostProb(Indorigin6gains_1_10)
prior6gains_1_10<-getPriorProb(Indorigin6gains_1_10)
bf6gains_1_10
posterior6gains_1_10
prior6gains_1_10

########7gains
Indorigin7gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=7)
bf7gains_1_10<-getBF(Indorigin7gains_1_10)
posterior7gains_1_10<- getPostProb(Indorigin7gains_1_10)
prior7gains_1_10<-getPriorProb(Indorigin7gains_1_10)
bf7gains_1_10
posterior7gains_1_10
prior7gains_1_10

########8gains
Indorigin8gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=8)
bf8gains_1_10<-getBF(Indorigin8gains_1_10)
posterior8gains_1_10<- getPostProb(Indorigin8gains_1_10)
prior8gains_1_10<-getPriorProb(Indorigin8gains_1_10)
bf8gains_1_10
posterior8gains_1_10
prior8gains_1_10

########9gains
Indorigin9gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=9)
bf9gains_1_10<-getBF(Indorigin9gains_1_10)
posterior9gains_1_10<- getPostProb(Indorigin9gains_1_10)
prior9gains_1_10<-getPriorProb(Indorigin9gains_1_10)
bf9gains_1_10
posterior9gains_1_10
prior9gains_1_10

########10gains
Indorigin10gains_1_10=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=10, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=10)
bf10gains_1_10<-getBF(Indorigin10gains_1_10)
posterior10gains_1_10<- getPostProb(Indorigin10gains_1_10)
prior10gains_1_10<-getPriorProb(Indorigin9gains_1_10)
bf10gains_1_10
posterior10gains_1_10
prior10gains_1_10


#######******** rate gain:loss=1:1

########1gain
Indorigin1gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=1)
bf1gains_1_1<-getBF(Indorigin1gains_1_1)
posterior1gains_1_1<- getPostProb(Indorigin1gains_1_1)
prior1gains_1_1<-getPriorProb(Indorigin1gains_1_1)
bf1gains_1_1
posterior1gains_1_1
prior1gains_1_1

########2gains
Indorigin2gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=2)
bf2gains_1_1<-getBF(Indorigin2gains_1_1)
posterior2gains_1_1<- getPostProb(Indorigin2gains_1_1)
prior2gains_1_1<-getPriorProb(Indorigin2gains_1_1)
bf2gains_1_1
posterior2gains_1_1
prior2gains_1_1

########3gains
Indorigin3gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=3)
bf3gains_1_1<-getBF(Indorigin3gains_1_1)
posterior3gains_1_1<- getPostProb(Indorigin3gains_1_1)
prior3gains_1_1<-getPriorProb(Indorigin3gains_1_1)
bf3gains_1_1
posterior3gains_1_1
prior3gains_1_1

########4gains
Indorigin4gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=4)
bf4gains_1_1<-getBF(Indorigin4gains_1_1)
posterior4gains_1_1<- getPostProb(Indorigin4gains_1_1)
prior4gains_1_1<-getPriorProb(Indorigin4gains_1_1)
bf4gains_1_1
posterior4gains_1_1
prior4gains_1_1

########5gains
Indorigin5gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=5)
bf5gains_1_1<-getBF(Indorigin5gains_1_1)
posterior5gains_1_1<- getPostProb(Indorigin5gains_1_1)
prior5gains_1_1<-getPriorProb(Indorigin5gains_1_1)
bf5gains_1_1
posterior5gains_1_1
prior5gains_1_1

########6gains
Indorigin6gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=6)
bf6gains_1_1<-getBF(Indorigin6gains_1_1)
posterior6gains_1_1<- getPostProb(Indorigin6gains_1_1)
prior6gains_1_1<-getPriorProb(Indorigin6gains_1_1)
bf6gains_1_1
posterior6gains_1_1
prior6gains_1_1

########7gains
Indorigin7gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=7)
bf7gains_1_1<-getBF(Indorigin7gains_1_1)
posterior7gains_1_1<- getPostProb(Indorigin7gains_1_1)
prior7gains_1_1<-getPriorProb(Indorigin7gains_1_1)
bf7gains_1_1
posterior7gains_1_1
prior7gains_1_1

########8gains
Indorigin8gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=8)
bf8gains_1_1<-getBF(Indorigin8gains_1_1)
posterior8gains_1_1<- getPostProb(Indorigin8gains_1_1)
prior8gains_1_1<-getPriorProb(Indorigin8gains_1_1)
bf8gains_1_1
posterior8gains_1_1
prior8gains_1_1

########9gains
Indorigin9gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=9)
bf9gains_1_1<-getBF(Indorigin9gains_1_1)
posterior9gains_1_1<- getPostProb(Indorigin9gains_1_1)
prior9gains_1_1<-getPriorProb(Indorigin9gains_1_1)
bf9gains_1_1
posterior9gains_1_1
prior9gains_1_1

########10gains
Indorigin10gains_1_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=1, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=10)
bf10gains_1_1<-getBF(Indorigin10gains_1_1)
posterior10gains_1_1<- getPostProb(Indorigin10gains_1_1)
prior10gains_1_1<-getPriorProb(Indorigin9gains_1_1)
bf10gains_1_1
posterior10gains_1_1
prior10gains_1_1

#######******** rate gain:loss=10:1

########1gains
Indorigin1gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=1)
bf1gains_10_1<-getBF(Indorigin1gains_10_1)
posterior1gains_10_1<- getPostProb(Indorigin1gains_10_1)
prior1gains_10_1<-getPriorProb(Indorigin1gains_10_1)
bf1gains_10_1
posterior1gains_10_1
prior1gains_10_1

########2gains
Indorigin2gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=2)
bf2gains_10_1<-getBF(Indorigin2gains_10_1)
posterior2gains_10_1<- getPostProb(Indorigin2gains_10_1)
prior2gains_10_1<-getPriorProb(Indorigin2gains_10_1)
bf2gains_10_1
posterior2gains_10_1
prior2gains_10_1

########3gains
Indorigin3gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=3)
bf3gains_10_1<-getBF(Indorigin3gains_10_1)
posterior3gains_10_1<- getPostProb(Indorigin3gains_10_1)
prior3gains_10_1<-getPriorProb(Indorigin3gains_10_1)
bf3gains_10_1
posterior3gains_10_1
prior3gains_10_1

########4gains
Indorigin4gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=4)
bf4gains_10_1<-getBF(Indorigin4gains_10_1)
posterior4gains_10_1<- getPostProb(Indorigin4gains_10_1)
prior4gains_10_1<-getPriorProb(Indorigin4gains_10_1)
bf4gains_10_1
posterior4gains_10_1
prior4gains_10_1

########5gains
Indorigin5gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=5)
bf5gains_10_1<-getBF(Indorigin5gains_10_1)
posterior5gains_10_1<- getPostProb(Indorigin5gains_10_1)
prior5gains_10_1<-getPriorProb(Indorigin5gains_10_1)
bf5gains_10_1
posterior5gains_10_1
prior5gains_10_1

########6gains
Indorigin6gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=6)
bf6gains_10_1<-getBF(Indorigin6gains_10_1)
posterior6gains_10_1<- getPostProb(Indorigin6gains_10_1)
prior6gains_10_1<-getPriorProb(Indorigin6gains_10_1)
bf6gains_10_1
posterior6gains_10_1
prior6gains_10_1

########7gains
Indorigin7gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=7)
bf7gains_10_1<-getBF(Indorigin7gains_10_1)
posterior7gains_10_1<- getPostProb(Indorigin7gains_10_1)
prior7gains_10_1<-getPriorProb(Indorigin7gains_10_1)
bf7gains_10_1
posterior7gains_10_1
prior7gains_10_1

########8gains
Indorigin8gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=8)
bf8gains_10_1<-getBF(Indorigin8gains_10_1)
posterior8gains_10_1<- getPostProb(Indorigin8gains_10_1)
prior8gains_10_1<-getPriorProb(Indorigin8gains_10_1)
bf8gains_10_1
posterior8gains_10_1
prior8gains_10_1

########9gains
Indorigin9gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=9)
bf9gains_10_1<-getBF(Indorigin9gains_10_1)
posterior9gains_10_1<- getPostProb(Indorigin9gains_10_1)
prior9gains_10_1<-getPriorProb(Indorigin9gains_10_1)
bf9gains_10_1
posterior9gains_10_1
prior9gains_10_1

########10gains
Indorigin10gains_10_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=10, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=10)
bf10gains_10_1<-getBF(Indorigin10gains_10_1)
posterior10gains_10_1<- getPostProb(Indorigin10gains_10_1)
prior10gains_10_1<-getPriorProb(Indorigin9gains_10_1)
bf10gains_10_1
posterior10gains_10_1
prior10gains_10_1

#######******** rate gain:loss=100:1

########1gains
Indorigin1gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=1)
bf1gains_100_1<-getBF(Indorigin1gains_100_1)
posterior1gains_100_1<- getPostProb(Indorigin1gains_100_1)
prior1gains_100_1<-getPriorProb(Indorigin1gains_100_1)
bf1gains_100_1
posterior1gains_100_1
prior1gains_100_1

########2gains
Indorigin2gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=2)
bf2gains_100_1<-getBF(Indorigin2gains_100_1)
posterior2gains_100_1<- getPostProb(Indorigin2gains_100_1)
prior2gains_100_1<-getPriorProb(Indorigin2gains_100_1)
bf2gains_100_1
posterior2gains_100_1
prior2gains_100_1

########3gains
Indorigin3gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=3)
bf3gains_100_1<-getBF(Indorigin3gains_100_1)
posterior3gains_100_1<- getPostProb(Indorigin3gains_100_1)
prior3gains_100_1<-getPriorProb(Indorigin3gains_100_1)
bf3gains_100_1
posterior3gains_100_1
prior3gains_100_1

########4gains
Indorigin4gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=4)
bf4gains_100_1<-getBF(Indorigin4gains_100_1)
posterior4gains_100_1<- getPostProb(Indorigin4gains_100_1)
prior4gains_100_1<-getPriorProb(Indorigin4gains_100_1)
bf4gains_100_1
posterior4gains_100_1
prior4gains_100_1

########5gains
Indorigin5gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=5)
bf5gains_100_1<-getBF(Indorigin5gains_100_1)
posterior5gains_100_1<- getPostProb(Indorigin5gains_100_1)
prior5gains_100_1<-getPriorProb(Indorigin5gains_100_1)
bf5gains_100_1
posterior5gains_100_1
prior5gains_100_1

########6gains
Indorigin6gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=6)
bf6gains_100_1<-getBF(Indorigin6gains_100_1)
posterior6gains_100_1<- getPostProb(Indorigin6gains_100_1)
prior6gains_100_1<-getPriorProb(Indorigin6gains_100_1)
bf6gains_100_1
posterior6gains_100_1
prior6gains_100_1

########7gains
Indorigin7gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=7)
bf7gains_100_1<-getBF(Indorigin7gains_100_1)
posterior7gains_100_1<- getPostProb(Indorigin7gains_100_1)
prior7gains_100_1<-getPriorProb(Indorigin7gains_100_1)
bf7gains_100_1
posterior7gains_100_1
prior7gains_100_1

########8gains
Indorigin8gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=8)
bf8gains_100_1<-getBF(Indorigin8gains_100_1)
posterior8gains_100_1<- getPostProb(Indorigin8gains_100_1)
prior8gains_100_1<-getPriorProb(Indorigin8gains_100_1)
bf8gains_100_1
posterior8gains_100_1
prior8gains_100_1

########9gains
Indorigin9gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=9)
bf9gains_100_1<-getBF(Indorigin9gains_100_1)
posterior9gains_100_1<- getPostProb(Indorigin9gains_100_1)
prior9gains_100_1<-getPriorProb(Indorigin9gains_100_1)
bf9gains_100_1
posterior9gains_100_1
prior9gains_100_1

########10gains
Indorigin10gains_100_1=testIndOrigin(inputTrees = c(pruned.trees),traitData = eyesVec,initLambda01 = .01, initLambda10=0.1, priorAlpha01=1, priorBeta01=1, priorAlpha10=1, priorBeta10=100, mcmcSize=11100, mcmcBurnin=100, mcmcSubsample=1, mcSize=10000, testThreshold=10)
bf10gains_100_1<-getBF(Indorigin10gains_100_1)
posterior10gains_100_1<- getPostProb(Indorigin10gains_100_1)
prior10gains_100_1<-getPriorProb(Indorigin9gains_100_1)
bf10gains_100_1
posterior10gains_100_1
prior10gains_100_1




