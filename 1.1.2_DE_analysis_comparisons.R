library("RColorBrewer")
library("tidyverse")
library("ggrepel")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

## Upset plot ----
a <- read.delim("0_Input_DE_list.tsv", check.names=F)
b <- list()

for(i in 1:ncol(a)){
  b[[i]] <- a[,i][is.na(a[,i])==F]
  names(b)[i] <- names(a)[i]
}

UpSetR::upset(UpSetR::fromList(b), nsets = 10, number.angles = 0, order.by = "freq", set_size.numbers_size = 15, text.scale = c(1.3, 1.3, 1.3, 1.3, 1.2, 1.3), sets.bar.color = "black", att.color = "black", main.bar.color = "black")  

## Volcano plot ----
pal2 <- colorRampPalette(brewer.pal(9, "RdBu"))(256)

ggplot(a, aes(x=-log2FC, y=-log10(padj), fill=-log2FC)) + 
  geom_point(shape=21, alpha=0.75, col="black") +
  theme_bw() +
  scale_fill_gradientn(colours=rev(pal2), limits=c(-12,12)) +
  geom_text_repel(aes(label=ifelse(-log10(padj) > 30, as.character(ID),'')),
                  size = 2.5, max.overlaps=30, box.padding = 0.25, lineheight = 0.5)+
  facet_wrap(~Comparison)+
  xlim(c(-12,12)) +
  ylim(c(0, 120))+
  labs(x="log2FC - E vs S", y="-log10(adjusted p-value)", col="log2FC", size="log10(Median level)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red")

### Dotplot -----
pal2 <- colorRampPalette(brewer.pal(9, "RdBu"))(256)

a$ID <- factor(a$ID, levels=rev(a$ID))

colorh <- c(pal2[250], pal2[5])

ggplot(a, aes(x=-q, y=ID, fill=GeneList, size=X.GeneInGOAndHitList)) + 
  geom_point(pch=21, col="black") +
  theme_bw() +
  scale_fill_manual(values = colorh)+
  facet_wrap(~GeneList, scale="free", nrow=2)+
  labs(x="-log10(adjusted p-value)", y="GO Biological Processes",
       size="N° DE genes", fill="Regulation")

### CLMM -------
library(tidyverse)
library(ordinal)
library(emmeans)

dat <- read.delim("Input.txt", header = TRUE, check.names = FALSE)

dat_long <- dat %>%
  pivot_longer(cols = c(`100 nM`, `1 µM`, `10 µM`, `100 µM`, `1 mM`), names_to = "Concentration", values_to = "Score")

dat_long$Concentration <- factor(dat_long$Concentration, levels = c("100 nM", "1 µM", "10 µM", "100 µM", "1 mM"), ordered = FALSE)

dat_long$Condition <- factor(dat_long$Condition, levels = c("E.coli", "S", "M"))

dat_long$Experiment <- factor(dat_long$Experiment)
dat_long$Score_ord <- factor(dat_long$Score, levels = c(0, 0.5, 1), ordered = TRUE)

mod_additive <- clmm(
  Score_ord ~ Condition + Concentration + (1 | Experiment),
  data = dat_long,
  link = "logit",
  Hess = TRUE,
  nAGQ = 10
)

anova(mod_additive, mod_full)

emm_condition <- emmeans(mod_full, ~ Condition | Concentration, mode = "latent")
pairs(emm_condition, adjust = "holm")

### Bayesian model -----
library("tidyverse")
library("brms")
library("emmeans")
library("posterior")

# Import
dat <- read.delim("Input.txt", header = TRUE,check.names = FALSE)

# Wide -> long
dat_long <- dat %>%
  pivot_longer(
    cols = c(`100 nM`, `1 µM`, `10 µM`, `100 µM`, `1 mM`),
    names_to = "Concentration",
    values_to = "Score"
  ) %>%
  mutate(Experiment = factor(Experiment),
         Condition = factor(Condition, levels = c("E.coli", "S", "M")),
         Concentration = factor(Concentration, levels = c("100 nM", "1 µM", "10 µM", "100 µM", "1 mM")),
    
# brms ordinal response
Score_ord = ordered(Score, levels = c(0, 0.5, 1)))

priors <- c(prior(normal(0, 2.5), class = "b"))

fit_bayes_adjusted <- brm(
  Score_ord ~ Condition + Concentration,
  data = dat_long,
  family = cumulative(link = "logit",  threshold = "flexible"),
  prior = priors, chains = 4, iter = 4000, warmup = 2000, cores = 4, 
  seed = 12345, control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(fit_bayes_adjusted)

emm_bayes <- emmeans(fit_bayes_adjusted, ~ Condition)
contr_bayes <- pairs(emm_bayes)

summary(contr_bayes, point.est = median, level = 0.95)
res <- as.data.frame(summary(contr_bayes, point.est = median, level = 0.95))

res$OR <- exp(res$estimate)
res$OR_lower <- exp(res$lower.HPD)
res$OR_upper <- exp(res$upper.HPD)