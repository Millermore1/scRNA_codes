library(tidyverse)

DESs <- read.csv("G:/sn_RNAseq_data/PRMT8KO_snRNA_seq/PV_neurons_Gqvsvector_GO_enrichment_G1.csv")
## 选择展示前几个 GO term
top_n_terms <- 10

plot_df <- DESs %>%
  mutate(
    GeneCount = str_count(Genes, "\\|") + 1,
    gene_label = str_replace_all(Genes, "\\|", ", ")
  ) %>%
  arrange(desc(Z.score)) %>%
  slice_head(n = top_n_terms) %>%
  mutate(
    Description = str_to_sentence(Description),
    y = rev(seq_len(n())),
    term_y = y + 0.22,
    gene_y = y - 0.22
  )


## 如果基因太多，可以只显示前若干个
max_genes_show <- 8

plot_df <- plot_df %>%
  mutate(
    gene_label = map_chr(str_split(gene_label, ",\\s*"), function(x) {
      if(length(x) > max_genes_show) {
        paste0(paste(x[1:max_genes_show], collapse = ", "), ", ...")
      } else {
        paste(x, collapse = ", ")
      }
    })
  )

x_max <- max(plot_df$Z.score) * 1.08

p <- ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = Z.score,
      ymin = y + 0.03,
      ymax = y + 0.47
    ),
    fill = "#FFDE7B",
    colour = NA
  ) +
  geom_text(
    aes(
      x = 0.05,
      y = term_y,
      label = Description
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.2,
    fontface = "bold",
    family = "Arial",
    colour = "black"
  ) +
  geom_text(
    aes(
      x = 0.05,
      y = gene_y,
      label = gene_label
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3,
    fontface = "italic",
    family = "Arial",
    colour = "black"
  ) +
  scale_x_continuous(
    limits = c(0, x_max),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, max(plot_df$y) + 0.65),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.4, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 10),
    plot.margin = margin(3, 3, 3, 3)
  )

p

ggsave(
  "GO_bar_with_genes_PV_Gq_G1.pdf",
  p,
  width = 120,
  height = 100,
  units = "mm",
  device = cairo_pdf
)

getwd()
