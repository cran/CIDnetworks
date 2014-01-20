
#library(CIDnetworks)

data(dolphins)

DIC.vector <- c()
pdf ("all-in.pdf")
model.plain <- CID.Gibbs (edge.list=dolphins, burnin=200, thin=5)
par(mfrow=c(1,2))
plot(model.plain)
summary(model.plain)
DIC.vector <- c(DIC.vector, plain=model.plain$DIC)

model.lsm <- CID.Gibbs (edge.list=dolphins, components=list(LSM(2)), burnin=1000, thin=5)
par(mfrow=c(1,3))
plot(model.lsm)
summary(model.lsm)
DIC.vector <- c(DIC.vector, lsm=model.lsm$DIC)

model.sr <- CID.Gibbs (edge.list=dolphins, components=list(SR()), burnin=200, thin=5)
par(mfrow=c(1,3))
plot(model.sr)
summary(model.sr)
DIC.vector <- c(DIC.vector, sr=model.sr$DIC)

model.sbm <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3)), burnin=200, thin=5)
par(mfrow=c(1,3))
plot(model.sbm)
summary(model.sbm)
DIC.vector <- c(DIC.vector, sbm=model.sbm$DIC)


model.lsm.sr <- CID.Gibbs (edge.list=dolphins, components=list(SR(), LSM(2)), burnin=200, thin=5)
par(mfrow=c(2,2))
plot(model.lsm.sr)
summary(model.lsm.sr)
DIC.vector <- c(DIC.vector, lsm.sr=model.lsm.sr$DIC)

model.sbm.lsm <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3), LSM(2)), burnin=200, thin=5)
par(mfrow=c(2,2))
plot(model.sbm.lsm)
summary(model.sbm.lsm)
DIC.vector <- c(DIC.vector, sbm.lsm=model.sbm.lsm$DIC)

model.sbm.sr <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3), SR()), burnin=200, thin=5)
par(mfrow=c(2,2))
plot(model.sbm.sr)
summary(model.sbm.sr)
DIC.vector <- c(DIC.vector, sbm.sr=model.sbm.sr$DIC)

model.sbm.lsm.sr <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3), LSM(2), SR()), burnin=1000, thin=5)
par(mfrow=c(2,3))
plot(model.sbm.lsm.sr)
summary(model.sbm.lsm.sr)
DIC.vector <- c(DIC.vector, sbm.lsm.sr=model.sbm.lsm.sr$DIC)

model.sbm.lsm1.sr <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3), LSM(1), SR()), burnin=1000, thin=5)
par(mfrow=c(2,3))
plot(model.sbm.lsm1.sr)
summary(model.sbm.lsm1.sr)
DIC.vector <- c(DIC.vector, sbm.lsm1.sr=model.sbm.lsm1.sr$DIC)

model.sbm.lvm <- CID.Gibbs (edge.list=dolphins, components=list(SBM(3), LVM(2)), burnin=1000, thin=5)
par(mfrow=c(2,3))
plot(model.sbm.lvm)
summary(model.sbm.lvm)
DIC.vector <- c(DIC.vector, sbm.lvm=model.sbm.lvm$DIC)

dev.off()



print(DIC.vector)
