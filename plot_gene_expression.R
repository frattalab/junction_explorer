# plot_gene_expression.R

# Load necessary libraries inside the function
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
library(forcats)
library(ggplot2)
library(ggpubr)
library(reshape2)

plot_expression <- function(this_gene, rna_expression_nmd, sh_exp, sk_exp, upf1_rna, chx_rna) {

    # NMD fractionation plot
    p1 <- rna_expression_nmd %>% 
        filter(gene_name == this_gene) %>% 
        mutate(significant_change = (padj < 0.1)) %>%
        mutate(frac = case_when(grepl("ratio", condition) ~ "Nuclear vs Cytoplasmic expression",
                                grepl("cytoplasmic", condition) ~ "Cytoplasmic fraction",
                                grepl("nuclear", condition) ~ "Nuclear fraction")) %>%
        mutate(cond2 = case_when(condition %in% c("nuclear_frac", "cytoplasmic_frac") ~ "TDP-43 KD v control",
                                 condition == 'cytoplasmic_nuclear_ratio_KD' ~ "Nuclear vs Cytoplasmic TDP-43 KD",
                                 condition == 'cytoplasmic_nuclear_ratio_Control' ~ "Nuclear vs Cytoplasmic Control",
                                 condition == 'cytoplasmic_frac_inKD_afterNMDi' ~ "TDP-43 KD v TDP-43 KD NMDi",
                                 condition == 'nuclear_frac_inKD_afterNMDi' ~ "TDP-43 KD v TDP-43 KD NMDi")) %>%
        ggplot(aes(x = cond2, y = log2FoldChange, fill = significant_change)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        coord_flip() +
        ggtitle("SMG1i + Fractionation") +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = c("TRUE" = "#9E2A2B", "FALSE" = "#FFF3B0")) +
        theme(legend.position = 'bottom') +
        labs(fill = "Adj. P-value < 0.1") +
        facet_wrap(~frac, scales = 'free', nrow = 3) + 
        ylab("Log2FC") + 
        xlab(element_blank())
    
    p2 <- sk_exp %>%
        filter(symbol %in% c("TARDBP", this_gene)) %>%
        mutate(significant_change = (padj < 0.1)) %>%
        select(symbol, log2FoldChange, source, significant_change) %>%
        pivot_wider(names_from = 'symbol', values_from = c("log2FoldChange")) %>%
        reshape2::melt(id.vars = c("source", "TARDBP", "significant_change")) %>%
        ggplot(aes(x = TARDBP, y = value, fill = significant_change)) +
        geom_point(pch = 21, size = 4) +
        geom_hline(yintercept = 0) +
        ggtitle(paste0("SK-N-BE(2) dose curve ", this_gene)) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = c("TRUE" = "#9E2A2B", "FALSE" = "#FFF3B0")) +
        theme(legend.position = 'bottom') +
        labs(fill = "Adj. P-value < 0.1") +
        ggpubr::stat_cor(method = 'spearman', cor.coef.name = 'rho') + 
        ylab("Log2FC after TDP-43 KD") +
        xlab("TARDBP Log2FC")
    
    p3 <- sh_exp %>%
        filter(symbol %in% c("TARDBP", this_gene)) %>%
        mutate(significant_change = (padj < 0.1)) %>%
        select(symbol, log2FoldChange, source, significant_change) %>%
        pivot_wider(names_from = 'symbol', values_from = c("log2FoldChange")) %>%
        reshape2::melt(id.vars = c("source", "TARDBP", "significant_change")) %>%
        ggplot(aes(x = TARDBP, y = value, fill = significant_change)) +
        geom_point(pch = 21, size = 4) +
        geom_hline(yintercept = 0) +
        ggtitle(paste0("SH-SY5Y dose curve ", this_gene)) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = c("TRUE" = "#9E2A2B", "FALSE" = "#FFF3B0")) +
        theme(legend.position = 'bottom') +
        labs(fill = "Adj. P-value < 0.1") +
        ggpubr::stat_cor(method = 'spearman', cor.coef.name = 'rho') +
        ylab("Log2FC after TDP-43 KD") +
        xlab("TARDBP Log2FC")
    
    p4 <- upf1_rna %>%
        filter(gene_name == this_gene) %>%
        mutate(significant_change = (padj < 0.1)) %>%
        mutate(experiment = case_match(experiment, "ctrlctrl_ctrlTDP43" ~ "TDP-43 KD v Control",
                                       "ctrlTDP43_UPF1TDP43" ~ "TDP-43+UPF1 KD v TDP-43 KD",
                                       "UPF1ctrl_vs_ctrlctrl" ~ "UPF1 KD v Control")) %>%
        filter(!is.na(experiment)) %>%
        mutate(experiment = fct_relevel(experiment, "TDP-43+UPF1 KD v TDP-43 KD", "TDP-43 KD v Control", after = Inf)) %>%
        ggplot(aes(x = experiment, y = log2fold_change, fill = significant_change)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        coord_flip() +
        ggtitle("UPF1 KD i3 Cortical") +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = c("TRUE" = "#9E2A2B", "FALSE" = "#FFF3B0")) +
        theme(legend.position = 'bottom') +
        labs(fill = "Adj. P-value < 0.1", ylab = "Log2FC", xlab = element_blank()) + 
        xlab(element_blank())
    
    p5 <- chx_rna %>%
        filter(gene_name == this_gene) %>%
        mutate(significant_change = (padj < 0.1)) %>%
        mutate(experiment = case_match(experiment, "Control_TDP43KD|Control_Control" ~ "TDP-43 KD v Control",
                                       "Cycloheximide_TDP43KD|Control_TDP43KD" ~ "TDP-43 KD + CHX v TDP-43 KD",
                                       "Cycloheximide_Control|Control_Control" ~ "CHX v Control")) %>%
        filter(!is.na(experiment)) %>%
        ggplot(aes(x = experiment, y = log2fold_change, fill = significant_change)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        coord_flip() +
        ggtitle("CHX treatment SH-SY5Y") +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = c("TRUE" = "#9E2A2B", "FALSE" = "#FFF3B0")) +
        theme(legend.position = 'bottom') +
        labs(fill = "Adj. P-value < 0.1") +
        ylab("Log2FC") + 
        xlab(element_blank())
    
    # Arrange and return the combined plot
    my_plot = ggpubr::ggarrange(
        ggpubr::ggarrange(p2, p3, ncol = 2, nrow = 1),
        ggpubr::ggarrange(p4, p5, p1, ncol = 3, nrow = 1),
        nrow = 2
    )
    
    
    return(my_plot)
}