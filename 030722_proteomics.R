library(DEqMS)
library(limma)
library(matrixStats)
library(ggrepel)

proteomics_dataset <- as.data.frame(X20220617_JHB_proteomics_proteinonly)

norm_abundance_data <- seq(41, 49, 1)
print(norm_abundance_data)

# filter out low FDR and change rownames to protein accession

filtered_proteomics_dataset <- proteomics_dataset[proteomics_dataset$"Protein FDR Confidence: Combined" != "Low", ]

row.names(filtered_proteomics_dataset) <- filtered_proteomics_dataset$Accession

filtered_proteomics_dataset <- filtered_proteomics_dataset[norm_abundance_data]

# log-transform and remove NAs

filtered_proteomics_dataset <- log2(filtered_proteomics_dataset)
filtered_proteomics_dataset <- na.omit(filtered_proteomics_dataset)

# check median centering

boxplot(filtered_proteomics_dataset, las=2, main="Proteomics data")

# design table

conditions <- as.factor(c(rep("Ciclesonide", 3), rep("Control", 3), rep("OA", 3)))

design = model.matrix(~0+conditions)
colnames(design) = gsub("conditions", "", colnames(design))

# contrasts

contrast_list <- c("Ciclesonide-Control","OA-Control","Ciclesonide-OA")

contrast =  makeContrasts(contrasts=contrast_list, levels=design)

fit1 <- lmFit(filtered_proteomics_dataset, design)

fit2 <- contrasts.fit(fit1, contrasts = contrast)

fit3 <- eBayes(fit2)

# DEqMS analysis

count_columns = seq(41, 49, 1)

psm.count.table = data.frame(count = rowMins(as.matrix(proteomics_dataset[,count_columns])), row.names = proteomics_dataset$Accession)

fit3$count = psm.count.table[rownames(fit3$coefficients), "count"]

fit4 = spectraCounteBayes(fit3)

print(fit4)

# visualise output

VarianceBoxplot(fit4, n = 30, main = "Proteomics dataset", xlab = "PSM count")

VarianceScatterplot(fit4, main = "Proteomics dataset")

# extract results

DEqMS.results = outputResult(fit4, coef_col = 1)
head(fit4$coefficients)
head(DEqMS.results)

write.table(DEqMS.results,"DEqMS.results.Proteomics_dataset_v1.txt",sep = "\t",
            row.names = F, quote = F)

# volcano plot

DEqMS.results$log.sca.pval = -log10(DEqMS.results$sca.P.Value)

ggplot(DEqMS.results, aes(x = logFC, y = log.sca.pval)) + 
  geom_point(size = 0.5)+
  theme_bw(base_size = 16) + # change theme
  xlab(expression("log2(Ciclesonide/Control)")) + # x-axis label
  ylab(expression("-log10(P-value)")) + # y-axis label
  geom_vline(xintercept = c(-1,1), colour = "red") + # add fold change cutoffs
  geom_hline(yintercept = 3, colour = "red") + # add significance cutoffs
  geom_vline(xintercept = 0, colour = "black") + # add 0 lines
  scale_colour_gradient(low = "black", high = "black", guide = FALSE)+
  geom_text_repel(data = subset(DEqMS.results, abs(logFC) >  1 & log.sca.pval > 3),
                  aes(logFC, log.sca.pval, label = gene))
