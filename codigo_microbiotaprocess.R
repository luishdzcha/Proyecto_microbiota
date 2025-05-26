#Cargar paquetes necesarios y librerías
library(MicrobiotaProcess)

BiocManager::install("MicrobiotaProcess")
library(BiocManager)
library(phyloseq)

install.packages("ggpubr")
library(ggpubr)

install.packages("ggtree")
library(ggtree)
library(dplyr)
library(ggplot2)
library(HMP16SData)
library(rstatix)
library(SummarizedExperiment)



# Cargar datos desde phyloseq
data(GlobalPatterns)  # Carga el dataset de ejemplo GlobalPatterns incluido en phyloseq
ps <- GlobalPatterns  # Asigna el dataset a un objeto llamado ps
ps  # Muestra un resumen del objeto phyloseq (número de muestras, taxones, etc.)

# Convertir a formato phyloseq (Esto ya no se hace) 
ps_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(ps)
ps_phyloseq 

# ======================
# ANÁLISIS DE DIVERSIDAD ALFA
# ======================

# Calcular métricas de diversidad alfa (Shannon, Simpson, Chao1)
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1"))
alpha_div  # Muestra las métricas calculadas para cada muestra

# Combinar las métricas con los metadatos
alpha_div <- cbind(sample_data(ps), alpha_div)

# Gráfico de boxplot por grupo (SampleType)
ggboxplot(alpha_div, x = "SampleType", y = "Shannon", 
          fill = "SampleType", palette = "jco") +
  stat_compare_means()  # Test Kruskal-Wallis para diferencias entre grupos

# INTERPRETACIÓN:
# - Shannon: Mide diversidad considerando riqueza y uniformidad (valores altos = mayor diversidad)
# - Simpson: Mide dominancia (valores cercanos a 1 = dominancia por pocas especies)
# - Chao1: Estimador de riqueza de especies (considera especies raras)
# El test Kruskal-Wallis indica si hay diferencias significativas entre los grupos

# ======================
# ANÁLISIS DE DIVERSIDAD BETA
# ======================

# Calcular matriz de distancia (Bray-Curtis)
dist_bray <- phyloseq::distance(ps, method = "bray")  # Bray-Curtis: considera abundancia

# PCoA (Análisis de Coordenadas Principales)
pcoa <- ordinate(ps, method = "PCoA", distance = dist_bray)

# Gráfico de PCoA
plot_ordination(ps, pcoa, color = "SampleType") +
  geom_point(size = 3) +  # Puntos para cada muestra
  stat_ellipse() +       # Elipses de confianza por grupo
  theme_minimal()

# INTERPRETACIÓN:
# - PCoA muestra patrones de similitud entre muestras
# - Muestras más cercanas = comunidades más similares
# - Elipses que no se solapan indican diferencias significativas entre grupos

# ======================
# GRÁFICO DE BARRAS (COMPOSICIÓN)
# ======================

# Agrupar por Filum y transformar a abundancia relativa
ps_phylum <- tax_glom(ps, taxrank = "Phylum")  # Agrupa OTUs a nivel de Phylum
ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x) * 100)  # Convierte a porcentajes

# Gráfico de barras
plot_bar(ps_rel, fill = "Phylum") +
  facet_wrap(~SampleType, scales = "free_x") +  # Paneles separados por tipo de muestra
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# INTERPRETACIÓN:
# - Muestra la composición relativa de cada phylum por tipo de muestra
# - Altura de barras = abundancia relativa (%)
# - Patrones visuales de dominancia (ej: Bacteroidetes en heces)


# ======================
# ANÁLISIS ESTADÍSTICO DETALLADO
# ======================

# Resumen estadístico por grupo
resumen_alfa <- alpha_div %>%
  group_by(SampleType) %>%
  summarise(
    Shannon_mean = mean(Shannon),
    Shannon_sd = sd(Shannon),
    Simpson_mean = mean(1-Simpson), # Convertir a diversidad (1-Simpson)
    Simpson_sd = sd(1-Simpson),
    Chao1_mean = mean(Chao1),
    Chao1_sd = sd(Chao1)
  )
print(resumen_alfa)

# Pruebas estadísticas:
# Kruskal-Wallis para Shannon
kruskal.test(Shannon ~ SampleType, data = alpha_div) %>% broom::tidy()

# Test post-hoc Dunn (comparaciones por pares)
dunn_test(alpha_div, Shannon ~ SampleType, p.adjust.method = "BH") %>% 
  filter(p.adj < 0.05) # Mostrar solo diferencias significativas

# INTERPRETACIÓN:
# - Kruskal-Wallis: Evalúa si hay diferencias globales entre grupos
# - Dunn Test: Identifica qué pares de grupos difieren significativamente
# - p.adj < 0.05 indica diferencias estadísticamente significativas

