# creating a figure to show mutation types

# Load packages
library(ggplot2)
library(dplyr)

# Read data
muts <- read.delim("data/mutations.tsv", header = TRUE, stringsAsFactors = FALSE)

# Fix missense classification
muts$MutationType[
  muts$MutationType == "Missense" & 
    (muts$PolyPhen > 0.5 | muts$SIFT < 0.05)
] <- "Damaging Missense"

muts$MutationType[
  muts$MutationType == "Missense" & 
    (muts$PolyPhen <= 0.5 & muts$SIFT >= 0.05)
] <- "Tolerated Missense"

# Create Impact category
muts <- muts %>%
  mutate(Impact = case_when(
    MutationType %in% c("Damaging Missense", "Nonsense", 
                        "Frameshift Deletion", "Frameshift Insertion", "Splice") ~ "Damaging",
    MutationType == "Tolerated Missense" ~ "Tolerated",
    TRUE ~ "Other"
  ))

# rename frameshift insertion and frameshift deletion
muts$MutationType[muts$MutationType == 'Frameshift Insertion'] <- 'Insertion'
muts$MutationType[muts$MutationType == 'Frameshift Deletion'] <- 'Deletion'

# ---- Prepare INNER ring (detailed) ----
inner <- muts %>%
  dplyr::count(MutationType) %>%
  mutate(ring = "inner")

# ---- Prepare outer ring (collapsed) ----
outer <- dplyr::count(muts, Impact)
colnames(outer)[colnames(outer) == "Impact"] <- "MutationType"
outer$ring <- "outer"

# Combine
plot_data <- bind_rows(outer, inner)

# Compute positions
plot_data <- plot_data %>%
  group_by(ring) %>%
  mutate(
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0)
  )

png(file = 'plots/mut_types.png',
    width = 5, height = 5, units = 'in', res = 1000)

# ---- Plot ----
ggplot(plot_data) +
  geom_rect(aes(
    ymin = ymin, ymax = ymax,
    xmin = ifelse(ring == "inner", 0.5, 1.5),
    xmax = ifelse(ring == "inner", 1.5, 2),
    fill = MutationType
  ), color = "white") +
  
  coord_polar(theta = "y") +
  xlim(0, 2.5) +
  theme_void() +
  
  geom_text(aes(
    x = ifelse(ring == "inner", 1, 2),
    y = (ymin + ymax)/2,
    label = paste0(MutationType, "\n", n)
  ), size = 3.5) +
  
  scale_fill_manual(values = c(
    '#aaaaaa', '#77aadd', '#dddd00', '#ffff00',
    '#ffcc66', '#ff9911', '#88cc88', '#aaccee', 'red'
  )) +
  
  labs(
    title = "KANSL1 Mutation Types",
    subtitle = "n = 25"
  ) +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      margin = margin(b = 2)
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      margin = margin(b = -50)
    ),
    
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

dev.off()

