
library(dplyr)
library(readr)
library(xlsx)

lncRNA <- read_csv("GSE125866.csv")

# Display first few rows of the dataframe
head(lncRNA)

# filtering data with cutoff 
sig.limma <- lncRNA[abs(lncRNA$logFC) > 1 & lncRNA$adj.P.Val < 0.05,]
View(sig.limma)


cat("The total number of dysregulated lncRNA with limma model:",nrow(sig.limma))

write.csv(sig.limma,file ="DElncRNA.csv")

downregulated_lncRNA <- sig.limma %>%
  filter(logFC < 0)

cat("Total number of downregulated lncRNA:", nrow(downregulated_lncRNA), "\n")

write_csv(downregulated_lncRNA, "downregulated_lncRNA.csv")

upregulated_lncRNA <- sig.limma %>%
  filter(logFC > 0)

cat("Total number of upregulated lncRNA:", nrow(upregulated_lncRNA), "\n")

write_csv(upregulated_lncRNA, "upregulated_lncRNA.csv")

library(ggplot2)

vol_plot <- ggplot(data=lncRNA, aes(x=logFC, y= -log10(adj.P.Val))) + geom_point()
vol_plot

vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")   

# Change xlim() ----------------------------------------------------------------
# Manually specify x-axis limits   
vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") + 
  xlim(-15, 15) 

# Modify scale_x_continuous() --------------------------------------------------
vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-15, 15, 2)), # Modify x-axis tick intervals    
                     limits = c(-15,15)) 

# Create new categorical column ------------------------------------------------ 
diseased_vs_healthy <- lncRNA %>%
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
  labs(title = "Volcano plot of differentially expressed lncRNA in GSE125866",
       x = "log2(Fold Change)",
       y = "-log10(Adj.P.Value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 


