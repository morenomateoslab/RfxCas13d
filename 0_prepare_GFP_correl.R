
# input data
GFP.INFO <- "input/GFP_info.csv";
TIGER <- "input/human_GFP_tiger_effic_vs_TIGER.csv";
DEEP <- "input/human_GFP_DeepCas13_effic_vs_DeepCas13.csv";
RNATARG <- "input/human_GFP_tiger_effic_vs_rnaTargeting.csv";

algNames <- c("TIGER","DeepCas13","RNATargeting");

# GFP efficiency info
gfp.eff <- read.csv(GFP.INFO, sep="\t", h=T);
# tiger
tiger <- read.csv(TIGER, sep="\t", h=T)[,c("GuideSeq","GuideScore")];colnames(tiger)[2] <- "TIGER";
# deepcas13
deep <- read.csv(DEEP, sep="\t", h=T)[,c("GuideSeq","DeepScore")];colnames(deep)[2] <- "DeepCas13";
# rnatarg
rnatarg <- read.csv(RNATARG, sep="\t", h=T)[,c("GuideSeq","predicted_value_sigmoid")];colnames(rnatarg)[2] <- "RNATargeting";

# merge
m <- merge(gfp.eff, merge(tiger, merge(deep, rnatarg, by.x="GuideSeq", by.y="GuideSeq"), by.x="GuideSeq", by.y="GuideSeq"), by.x="GuideSeq", by.y="GuideSeq");

dir.create("OUT_0");write.table(m, file=paste("OUT_0/master_prediction_",nrow(m),"_GFP_tiling.csv",sep=""), sep="\t", row.names=F, col.names=T, quote=F);

# corr pearson
out.pears <- unlist(lapply(algNames,function(alg){cor(m$eff,m[,alg],method="pearson")}));
	
# pval pearson
pval.pears <- unlist(lapply(algNames,function(alg){cor.test(m$eff,m[,alg],method="pearson")$p.value}));

out <- as.data.frame(round(rbind(out.pears,pval.pears),3));out$Name <- c("r_pearson","pval_pearson");colnames(out)[1:3] <- algNames;
write.table(out[c(4,1:3)], file="OUT_0/GFP_corr.txt", sep="\t", row.names=F, col.names=T, quote=F);
