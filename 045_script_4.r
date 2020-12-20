#

# select disease-relevant SNPs


setwd("E:/freelance/applications/045_randforest")


# selection threshold
num_mdg_snp <- 300
num_oa_go <- 20
num_oa_snp <- 200


# d1
did <- "d1"
ria_fn <- "d1_plink1_om1.5_rf_impt_agre.csv"

# d2
did <- "d2"
ria_fn <- "d2_plink1_a1e-5_rf_impt_agre.csv"


# mean_MeanDecreaseGini
ria <- read.table(ria_fn, sep=",", header=TRUE, as.is="snp")
rownames(ria) <- ria[,"snp"]
top_mdg_snps <- ria[1:num_mdg_snp, "snp"]

# GO to SNP list
go2s_ls_fn <- paste("rap", did, "go2snp.rd", sep="_")
load(go2s_ls_fn)

# OOB error rate and AUC by GO
oa_fn <- paste("rap", did, "go_oob_auc.tsv", sep="_")
oa_df <- read.table(oa_fn, sep="\t", header=TRUE, as.is="GO")
oa_df[,"1-OOB+AUC"] <- 1 - oa_df[,"OOB"] + oa_df[,"AUC"]
oa_df <- oa_df[order(oa_df[,"1-OOB+AUC"], decreasing=TRUE),]
oa_fn_w1 <- paste("rap", did, "go_oob_auc_1.tsv", sep="_")
write.table(oa_df, file=oa_fn_w1, sep="\t", col.names=TRUE, row.names=FALSE)
# select by "1-OOB+AUC" then ranked by "mean_MeanDecreaseGini"
top_oa_go = oa_df[1:num_oa_go, "GO"]
top_oa_snps_ls <- go2snp_ls[top_oa_go]
top_oa_snps_uniq <- unique(unlist(top_oa_snps_ls))
ria_top_oa <- ria[top_oa_snps_uniq,]
top_oa_snps_uniq_ranked <- ria_top_oa[order(ria_top_oa[,"mean_MeanDecreaseGini"], decreasing=TRUE), "snp"]
top_oa_snps_final <- top_oa_snps_uniq_ranked[1:num_oa_snp]

# combine SNPs selected by MeanDecreaseGini and those by "1-OOB+AUC"
top_snps <- na.omit(unique(c(top_mdg_snps, top_oa_snps_final)))
ts_fn <- paste(paste(did, "top_snps", num_mdg_snp, num_oa_go, num_oa_snp, sep="_"), ".txt", sep="")
write.table(top_snps, file=ts_fn, col.names=FALSE, row.names=FALSE, quote=FALSE)


