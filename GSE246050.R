
library(dplyr)
library(readr)
library(xlsx)

mRNA <- read_csv("GSE246050.csv")

# Display first few rows of the dataframe
head(mRNA)

# filtering data with cutoff 
sig.limma <- mRNA[abs(mRNA$logFC) > 1 & mRNA$adj.P.Val < 0.05,]
View(sig.limma)


cat("The total number of dysregulated MRNA with limma model:",nrow(sig.limma))

write.csv(sig.limma,file ="DEMRNA_RNASeq_2.csv")

downregulated_mrna <- sig.limma %>%
  filter(logFC < 0)

cat("Total number of downregulated mRNA:", nrow(downregulated_mrna), "\n")

write_csv(downregulated_mrna, "downregulated_RNASeqMRNA_2.csv")

upregulated_mrna <- sig.limma %>%
  filter(logFC > 0)

cat("Total number of downregulated mRNA:", nrow(upregulated_mrna), "\n")

write_csv(upregulated_mrna, "upregulated_RNASeqMRNA_2.csv")

library(ggplot2)

vol_plot <- ggplot(data=mRNA, aes(x=logFC, y= -log10(adj.P.Val))) + geom_point()
vol_plot

vol_plot + 
  geom_hline(yintercept = -log10(0.5),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")   

# Change xlim() ----------------------------------------------------------------
# Manually specify x-axis limits   
vol_plot + 
  geom_hline(yintercept = -log10(0.5),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") + 
  xlim(-25, 25) 

# Modify scale_x_continuous() --------------------------------------------------
vol_plot + 
  geom_hline(yintercept = -log10(0.5),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-15, 15, 2)), # Modify x-axis tick intervals    
                     limits = c(-15,15)) 

# Create new categorical column ------------------------------------------------ 
diseased_vs_healthy <- mRNA %>%
  mutate(Gene_Type = case_when(logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
                               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   

# Obtain gene_type counts ------------------------------------------------------           
diseased_vs_healthy %>%
  count(Gene_Type)


cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

diseased_vs_healthy %>%
  ggplot(aes(x = logFC,
             y = -log10(adj.P.Val),
             fill = Gene_Type,    
             size = Gene_Type,
             alpha = Gene_Type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-15, 15, 2)),       
                     limits = c(-15, 15))+
  labs(title = "Volcano plot of differentially expressed mRNA in GSE246050",
       x = "log2(Fold Change)",
       y = "-log10(Adj.P.Value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 


