# helper function
extract.last <- function(x, pattern = "_"){
        extract.last.one <- function(y, pattern = "_"){tmp <- unlist(strsplit(as.character(y), pattern));return(tmp[length(tmp)]);}
        unlist(lapply(x, function(y){extract.last.one(y,pattern)}));
}

# input data
GUIDE.INFO <- "input/guide_info.csv";
TIGER <- "input/TIGER_on_our_gRNAs.csv";
DEEP <- "input/deepCas13_on_our_gRNAs.csv";
RNATARG <- "input/rnaTargeting_on_our_gRNAs.csv";
GFP.COR <- "OUT_0/GFP_corr.txt";

algNames <- c("TIGER","DeepCas13","RNATargeting");

# guide efficiency info
guide.eff <- read.csv(GUIDE.INFO, sep="\t", h=T);
# tiger
tiger <- read.csv(TIGER, sep="\t", h=T)[,c("index","Guide.Score")];colnames(tiger)[2] <- "TIGER";
# deepcas13
deep <- read.csv(DEEP, sep="\t", h=T)[,c("index","DeepScore")];colnames(deep)[2] <- "DeepCas13";
# rnatarg
rnatarg <- read.csv(RNATARG, sep="\t", h=T)[,c("index","predicted_value_sigmoid")];colnames(rnatarg)[2] <- "RNATargeting";

# merge
m <- merge(guide.eff, merge(tiger, merge(deep, rnatarg, by.x="index", by.y="index"), by.x="index", by.y="index"), by.x="index", by.y="index");

dir.create("OUT_1");write.table(m, file=paste("OUT_1/master_prediction_",nrow(m),"_gRNAs.csv",sep=""), sep="\t", row.names=F, col.names=T, quote=F);

# by cocktail
m$cocktail <- as.integer(extract.last(m$index));out.pears <- pval.pears <- list();
for(cocktail in 1:length(unique(m$cocktail))){
	tmp <- m[m$cocktail==cocktail,];

	# corr pearson
	out.pears[[paste("cocktail",cocktail,sep="_")]] <- unlist(lapply(algNames,function(alg){cor(tmp$eff,tmp[,alg],method="pearson")}));
	
	# pval pearson
	pval.pears[[paste("cocktail",cocktail,sep="_")]] <- unlist(lapply(algNames,function(alg){cor.test(tmp$eff,tmp[,alg],method="pearson")$p.value}));
}

# pearson
pearsMat <- t(matrix(unlist(out.pears), nrow=length(out.pears[[1]]), ncol=length(out.pears)));
pearsPvalMat <- t(matrix(unlist(pval.pears), nrow=length(pval.pears[[1]]), ncol=length(pval.pears)));
colnames(pearsMat) <- colnames(pearsPvalMat) <- algNames;
rownames(pearsMat) <- rownames(pearsPvalMat) <- names(out.pears);

# add GFP row
gfp <- read.csv(GFP.COR, sep="\t", h=T);
pearsMat <- rbind(pearsMat, gfp[gfp$Name=="r_pearson",algNames]);
pearsPvalMat <- rbind(pearsPvalMat, gfp[gfp$Name=="pval_pearson",algNames]);

rownames(pearsMat)[nrow(pearsMat)] <- rownames(pearsPvalMat)[nrow(pearsPvalMat)] <- "GFP_tiling";

library(pheatmap);
svg("OUT_1/r_pearson_heatmap.svg");
	pheatmap(pearsMat, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=T, fontsize_number=16, fontsize_row=16, fontsize_col=16, main="R Pearson", fontsize=20);
dev.off();

svg("OUT_1/pval_pearson_heatmap.svg");
        pheatmap(pearsPvalMat, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=T, fontsize_number=16, fontsize_row=16, fontsize_col=16, main="P-val Pearson", fontsize=20);
dev.off();
