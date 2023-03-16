# Data analysis for the Osmiini provision manipulation experiment.
# Data was collected in 2021 and 2022 by Sébastien Rivest at the Rocky Mountain Biological Laboratory.
# The effect of Asteraceae pollen spines and Lupinus pollen alkaloids on survival probability... 
# were investigated  on 9 species of bees in the Osmiini tribe and their cleptoparasite (sapyga. sp).
# Bee eggs were sampled from nests collected from artificial nesting structures ("Bloc") at 12 different sites.

# Load packages.
Packages <- c("ape", "nlme", "ggtree", "tidyverse", "readxl", "blmeco", "vegan", "tidyr",
              "ggplot2", "ggthemes", "lattice", "magrittr", "ggpubr", "gridExtra", "brms",
              "rstan", "bayesplot", "tidybayes", "ggridges", "glue", "stringr", "forcats",
              "mgcv", "pracma", "extrafont", "survival", "survminer")
lapply(Packages, library, character.only = TRUE)

# Load dataset.
osmiini_survival <- read_excel("C:/Users/seb69/OneDrive/Documents/Doc/egg-transfer_figs/Last versions/Osmiini_survival_data.xlsx")osmiini_survival <- read_excel("C:/Users/seb69/OneDrive/Documents/Doc/egg-transfer_figs/Last versions/Osmiini_survival_data.xlsx")
osmiini_survival <- osmiini_survival %>% 
  mutate(Nest = paste(.$Site, .$Bloc, .$Nest, .$Year, sep = "-"), # Unique nest and bloc IDs  
         Bloc = paste(.$Site, .$Bloc, .$Year, sep = "-")) 

# Build phylogenetic tree.
# Sapyga. sp is included only for visualization. 
# We considered that their is no phylogenetic correlation between this wasp and the osmiine bee species.
Osmini_tree <- read.tree(text='(Sapyga. sp:0.00000001, (H. fulgida:0.051, (O. lignaria:0.029,
  ((O. iridis:0.026, (O. montana:0.010, O. subaustralis:0.006):0.019):0.018:,
  (O. coloradensis:0.036, (O.tersula:0.014, 
  (O. tristella:0.010, O. pusilla:0.009):0.007):0.019):0.007):0.009):0.051):0.00000001);')
plot(Osmini_tree, cex=0.5)
# Phylogenetic correlation matrix
corr_tree <- ape::vcv.phylo(Osmini_tree, corr = TRUE)

# Model of the effect of Asteraceae spines in function of host-use type
# Extract data from the Asteraceae spines experiment.
survival_aster <- osmiini_survival %>% 
  subset(Experiment == "aster") %>% 
  mutate(grouping = paste(.$Species, .$Site, .$Bloc, .$Nest, sep = "-"),
         census_1 = 2-status,
         Year = as.factor(Year)) %>%
  subset(., Species != "Sapyga.sp")

# Bayesian model
aster_mod <- brm(days.to.event|cens(census_1) ~ Diet*treatment + Year +
                   (1|gr(Species, cov = corr_tree)) + (1|Site/Bloc/Nest),
                 family = brmsfamily("cox"),
                 prior = c(set_prior("normal(0,5)", class = "b"),
                           set_prior("student_t(3,0,2.5)", class = "sd")),
                 data = survival_aster, data2 = list(corr_tree = corr_tree),
                 warmup = 2500, iter = 10000, chains = 4, thin = 8, 
                 control = list(adapt_delta = 0.99), 
                 cores = 4, seed = 1234)
summary(aster_mod) # Note that values in the model summary correspond to log hazard ratios.
# Look for the presence of autocorrelation in the chains.
plot(aster_mod)
# Verify if the chains have converged.
mcmc_acf_bar(aster_mod, regex_pars = c("sd"))
# Test hypotheses for comparisons not shown in the model summary.
hypothesis(aster_mod, c("Dietfabaceae:treatmentT = Dietgeneralist:treatmentT"))

# Model of the effect of Lupinus alkaloids in function of host-use type
# Extract data from the Lupinus alkaloids experiment.
survival_lup <- osmiini_survival %>%
  subset(Experiment == "lup") %>%
  mutate(grouping = paste(.$Species, .$Site, .$Bloc, .$Nest, sep = "-"),
         census_1 = 2-status,
         Diet = factor(.$Diet, levels = c("generalist", "cleptoparasite", "fabaceae", "aster")),
         Year = as.factor(Year)) %>%
  subset(., Species != "Sapyga.sp")
# Bayesian model
lupinus_mod <- brm(days.to.event|cens(census_1) ~ treatment*Diet + Year +
                     (1|gr(Species, cov = corr_tree)) + (1|Site/Bloc/Nest),
                   family = brmsfamily("cox"),
                   prior = c(set_prior("normal(0,5)", class = "b"),
                             set_prior("student_t(3,0,2.5)", class = "sd")),
                   data = survival_lup, data2 = list(corr_tree = corr_tree),
                   warmup = 2500, iter = 10000, chains = 4, thin = 8, 
                   control = list(adapt_delta = 0.99),
                   cores = 4, seed = 1234)
summary(lupinus_mod)
plot(lupinus_mod)
mcmc_acf_bar(lupinus_mod, regex_pars = c("sd"))
hypothesis(lupinus_mod, c("treatmentT:Dietaster = treatmentT:Dietfabaceae"))

# Plots of Kaplan-Meier survival curves in function of host-use types
# Asteraceae spine experiment
# Asteraceae specialists
# Extract data for Asteraceae specialists.
surv_aster_aster <- survival_aster[survival_aster$Diet == "aster",] 
# Compute survival curve
fit_aster_aster <- survfit(Surv(days.to.event, status == 2) ~ treatment, 
                           data = surv_aster_aster)
# Plot
plot_aster_aster <- ggsurvplot(fit_aster_aster, data = surv_aster_aster, 
                               palette = c("#f69f56","#ca2e30"),
                               linetype = c(2,1), conf.int = TRUE, xlim = c(0,59.5))
# Generalists
surv_gen_aster <- survival_aster[survival_aster$Diet == "generalist",]
fit_gen_aster <- survfit(Surv(days.to.event, status == 2) ~ treatment, 
                         data = surv_gen_aster)
plot_gen_aster <- ggsurvplot(fit_gen_aster, data = surv_gen_aster,
                             palette = c("#f69f56","#ca2e30"),
                             linetype = c(2,1), conf.int = TRUE, xlim = c(0,59.5))
# Fabaceae specialists
surv_fab_aster <- survival_aster[survival_aster$Diet == "fabaceae",]
fit_fab_aster <- survfit(Surv(days.to.event, status == 2) ~ treatment, 
                         data = surv_fab_aster )
plot_fab_aster <- ggsurvplot(fit_fab_aster, data = surv_fab_aster,
                             palette = c("#f69f56","#ca2e30"),
                             linetype = c(2,1), conf.int = TRUE, xlim = c(0,59.5))
# Combine plots 
aster_plots <- list() 
aster_plots[[1]] <- plot_aster_aster
aster_plots[[2]] <- plot_fab_aster
aster_plots[[3]] <- plot_gen_aster
arrange_ggsurvplots(aster_plots, ncol = 3)

# Lupinus alkaloids experiment
# Asteraceae specialists
surv_aster_lup <- survival_lup[survival_lup$Diet == "aster",]
fit_aster_lup <- survfit(Surv(days.to.event, status == 2) ~ treatment, 
                         data = surv_aster_lup)
plot_aster_lup <- ggsurvplot(fit_aster_lup, data = surv_aster_lup,
                             palette = c("#77b1d9","#1d6391"),
                             linetype = c(2,1), conf.int = TRUE, xlim = c(0,47))
# Generalists
surv_gen_lup <- survival_lup[survival_lup$Diet == "generalist",]
fit_gen_lup <-survfit(Surv(days.to.event, status == 2) ~ treatment, 
                      data = surv_gen_lup)
plot_gen_lup <- ggsurvplot(fit_gen_lup, data = surv_gen_lup,
                           palette = c("#77b1d9","#1d6391"),
                           linetype = c(2,1), conf.int = TRUE, xlim = c(0,47))
# Fabaceae specialists
surv_fab_lup <- survival_lup[survival_lup$Diet == "fabaceae",]
fit_fab_lup <- survfit(Surv(days.to.event, status == 2) ~ treatment, 
                       data = surv_fab_lup)
plot_fab_lup <- ggsurvplot(fit_fab_lup, data = surv_fab_lup,
                           palette = c("#77b1d9","#1d6391"),
                           linetype = c(2,1), conf.int = TRUE, xlim = c(0,47))

lup_plots <- list()
lup_plots[[1]] <- plot_aster_lup
lup_plots[[2]] <- plot_gen_lup
lup_plots[[3]] <- plot_fab_lup
arrange_ggsurvplots(lup_plots, ncol = 3)

# Comparison of the effect of Asteraceae spines and Lupinus alkaloids on larval development between species
# Effect of Asteraceae spines
# Extract data from the Asteraceae spines experiment (with Sapyga sp).
survival_aster_all <- osmiini_survival %>% 
  subset(Experiment == "aster") %>% 
  mutate(grouping = paste(.$Species, .$Site, .$Bloc, .$Nest, sep = "-"),
         census_1 = 2-status,
         Year = as.factor(Year))

aster_species_mod <- brm(days.to.event|cens(census_1) ~ treatment + Year +
                           (1 + treatment|gr(Species, cov = corr_tree)) + (1|Site/Bloc/Nest),
                         family = brmsfamily("cox"),
                         prior = c(set_prior("normal(0,5)", class = "b"),
                                   set_prior("student_t(3,0,2.5)", class = "sd")),
                         data = survival_aster_all, data2 = list(corr_tree = corr_tree),
                         warmup = 2500, iter = 10000, chains = 4, thin = 8, 
                         control = list(adapt_delta = 0.99),
                         cores = 4, seed = 1234)
summary(aster_species_mod) # Note that values in the model summary corresponds to log hazard ratio.
plot(aster_species_mod)
mcmc_acf_bar(aster_species_mod, regex_pars = c("sd"))

# Plot of hazard ratio for all species
# Sample values from the model to draw posterior distributions.
draws<-spread_draws(aster_species_mod, r_Species[Species,term], b_treatmentT) %>% 
  filter(term == "treatmentT") %>% 
  mutate(b_treatmentT = r_Species + b_treatmentT)

pheno.data <- draws %>% 
  ungroup() %>%
  mutate(Species = str_replace_all(Species, "[.]", " ")) %>% 
  mutate(Species = reorder(Species, b_treatmentT))
# Summarize data and back transform posterior values from log(hazard ratio) to hazard ratio.
pheno.data.summary <- group_by(pheno.data, Species) %>% 
  mean_qi(b_treatmentT) %>%
  mutate(b_treatmentT = exp(b_treatmentT), # Back transform posterior values.
         .lower = exp(.lower),
         .upper = exp(.upper))

colors <- c("black", "grey", "#1d6391", "#1d6391", "#f69f56", "#f69f56",
            "#ca2e30", "#f69f56", "#1d6391", "#1d6391", "#1d6391")
species <- c("NA", "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
             "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla")
names(colors) <- species
# Plot
p_aster <- ggplot(aes(x = b_treatmentT, y = factor(Species, levels = c(
  "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
  "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla"))), 
  data = pheno.data) +
  geom_errorbarh(data = pheno.data.summary,
                 aes(fill="black",x = b_treatmentT, y = factor(Species, levels = c(
                   "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
                   "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla")),
                   xmin=.lower, xmax=.upper), 
                 height = 0.15, size = 0.5) +
  geom_point(data = pheno.data.summary,
             aes(fill = Species,
                 x = b_treatmentT, y = factor(Species, levels = c(
                   "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
                   "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla"))), 
             size = 5, col = "black", shape = 21) +
  scale_fill_manual(values = colors) +
  geom_vline(xintercept = 1, color = "black", linetype="dashed", size = 1, alpha = 0.5) +
  coord_cartesian(xlim = c(0, 21)) +
  labs(x = "Hazard ratio",
       y = element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.y = element_text(face="italic", size = 10))
p_aster

# Effect of Lupinus pollen alkaloids
survival_lup_all <- osmiini_survival %>%
  subset(Experiment == "lup") %>%
  mutate(grouping = paste(.$Species, .$Site, .$Bloc, .$Nest, sep = "-"),
         census_1 = 2-status,
         Diet = factor(.$Diet, levels = c("generalist", "cleptoparasite", "fabaceae", "aster")),
         Year = as.factor(Year))

Lupinus_species_mod  <- brm(days.to.event|cens(census_1) ~ treatment + Year +
                              (1 + treatment|gr(Species, cov = corr_tree)) + (1|Site/Bloc/Nest),
                            family = brmsfamily("cox"),
                            prior = c(set_prior("normal(0,5)", class = "b"),
                                      set_prior("student_t(3,0,2.5)", class = "sd")),
                            data = survival_lup_all, data2 = list(corr_tree = corr_tree),
                            warmup = 2500, iter = 10000, chains = 4, thin = 8, 
                            control = list(adapt_delta = 0.99),
                            cores = 4, seed = 1234)
summary(Lupinus_species_mod)
plot(Lupinus_species_mod)
mcmc_acf_bar(Lupinus_species_mod, regex_pars = c("sd"))

# Plot of hazard ratio for all species
draws<-spread_draws(Lupinus_species_mod, r_Species[Species,term], b_treatmentT) %>% 
  filter(term == "treatmentT") %>% 
  mutate(b_treatmentT = r_Species + b_treatmentT)

pheno.data <- draws %>% 
  ungroup() %>%
  mutate(Species = str_replace_all(Species, "[.]", " ")) %>% 
  mutate(Species = reorder(Species, b_treatmentT))

pheno.data.summary <- group_by(pheno.data, Species) %>% 
  mean_qi(b_treatmentT)

pheno.data.summary <- group_by(pheno.data, Species) %>% 
  mean_qi(b_treatmentT) %>%
  mutate(b_treatmentT = exp(b_treatmentT), # Back transform posterior values
         .lower = exp(.lower),
         .upper = exp(.upper))

p_lup <- ggplot(aes(x = b_treatmentT, y = factor(Species, levels = c(
  "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
  "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla"))), 
  data = pheno.data) +
  geom_errorbarh(data = pheno.data.summary,
                 aes(fill="black",x = b_treatmentT, y = factor(Species, levels = c(
                   "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
                   "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla")),
                   xmin=.lower, xmax=.upper), 
                 height = 0.15, size = 0.5) +
  geom_point(data = pheno.data.summary,
             aes(fill = Species,
                 x = b_treatmentT, y = factor(Species, levels = c(
                   "Sapyga sp", "H fulgida", "O lignaria", "O montana", "O subaustralis",
                   "O iridis",  "O coloradensis", "O tersula", "O tristella", "O pusilla"))), 
             size = 5, col = "black", shape = 21) +
  scale_fill_manual(values = colors) + 
  geom_vline(xintercept = 1, color = "black", linetype="dashed", size = 1, alpha = 0.5) +
  coord_cartesian(xlim = c(0, 21)) +
  labs(x = "Hazard ratio",
       y = element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.y = element_text(face="italic", size = 10))
p_lup

# Combine plots of Lupinus and Asteraceae experiments with Osmiini phylogeny.
p2<-ggtree(Osmini_tree, branch.length="none") %>% flip(4,16)
ggarrange(p2+ theme(plot.margin = unit(c(0.8,2.9,0.8,0.4), "cm")),
          p_lup + theme(plot.margin = unit(c(0.7,1,0,-3), "cm")),
          p_aster + theme(plot.margin = unit(c(0.7,1,0,-0.8), "cm",),
                          axis.text.y = element_blank()),
          nrow =1)
