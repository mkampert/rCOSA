### R code from vignette source 'TutorialCOSA.Rnw'

###################################################
### code chunk number 1: TutorialCOSA.Rnw:14-18
###################################################
require(COSA)
Xmat <- as.matrix(iris[,1:4])
d <- cosa(X = Xmat)                                                           
print(d)


