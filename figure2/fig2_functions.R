#########################################
# Figure 2 
#########################################
library(dplyr)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(sf)
library(ggpubr)

# NOTE, TODO: you must add the path to the data directory here
data_dir <- "..."

# Data for each state
nm_hh <- read.csv(file.path(data_dir, "aux", "nm.csv"))
fl_sch <- read.csv(file.path(data_dir, "aux", "fl_county_schooldata.csv"))
fl_sch$county_fips <- as.character(fl_sch$county_fips)
ny_emp <- read.csv(file.path(data_dir, "aux", "ny.csv"))

# some color palettes 
pal <- wes_palette("Zissou1", 7, type = "continuous")
pal1 <- wes_palette("Royal2", 7, type = "continuous")

# Household sizes
# ==============================================================================
A <- nm_hh %>%
  select(c(household_id, household_size, person_race)) %>%
  distinct() %>%
  mutate(person_race = as.factor(recode(person_race,
                                        "0" = "White",
                                        "1" = "Black",
                                        "2" = "Asian",
                                        "3" = "Native American",
                                        "4" = "Pacific Islander",
                                        "5" = "Other",
                                        "6" = "Multi"))) %>%
  select(c(household_size, person_race)) %>%
  filter(person_race == "White" | person_race == "Native American") %>%
  ggplot() +
  theme_minimal() +
  xlim(0, 10) +
  geom_density(aes(x = household_size, fill = person_race, 
                   color = person_race), adjust = 6, alpha = 0.25) +
  labs(title = paste("Household Size Distribution by Race, NM"), 
       x = "Household Size",
       y = "Density",
       fill = "Race",
       color = "Race") +
  theme(legend.title=element_blank(), text = element_text(size = 8))

# Schoolgroup sizes
# ==============================================================================
# Create maps of FL
florida_counties <- counties(state = "FL", cb = TRUE, year = 2020) %>%
  st_as_sf() # convert to sf object for mapping

florida_map <- florida_counties %>%
  left_join(fl_sch, by = c("GEOID" = "county_fips"))

# Mean student:teacher ratio in FL
b1 <- florida_map %>%
  ggplot() +
  geom_sf(aes(fill = mean_schoolgroup_size), color = "white", size = 0.2) +
  scale_fill_viridis_c(name = "", na.value = "gray90") +
  theme_minimal() +
  labs(
    title = "Mean Student:Teacher Ratio, FL"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    text = element_text(size = 8)
  )

# Percent Hispanic students by county in FL
b2 <- ggplot(florida_map) +
  geom_sf(aes(fill = pct_hispanic_students), color = "white", size = 0.2) +
  scale_fill_viridis_c(name = "", na.value = "gray90") +
  theme_minimal() +
  labs(
    title = "% Hispanic Students, FL"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    text = element_text(size = 8)
  )

# Correlation plot
b3 <- florida_map %>%
  ggplot() +
  aes(x = mean_schoolgroup_size, y = pct_hispanic_students) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, col = "#00BFC4") +
  stat_cor(method = "pearson", label.x = 20) +
  labs(
    x = "Mean student:teacher ratio",
    y = "% Hispanic students"
  ) + 
  theme(text = element_text(size = 8))

B <- plot_grid(b1, b2, b3, labels = c("B1", "B2", "B3"), nrow = 1)

# NAICS codes
# ==============================================================================
# race 
naics_tracts <- function(occ) {
  naics_names <- data.frame(
    code = c("111", "112", "113", "114", "445", "541", "622", "623", "722"),
    name = c("Crop Production", "Animal Production", "Forestry and Logging", 
             "Fishing and Hunting", "Food and Beverage Stores", 
             "Professional, Scientific, and Technical Services", "Hospitals", 
             "Nursing and Residential Care Facilities",
             "Food Services and Drinking Places")
  )
  
  # merge NAICS names
  naics <- ny_emp %>%
    filter(person_naics != 0) %>%
    mutate(person_naics = as.character(person_naics), 
           person_race = as.factor(person_race)) %>%
    left_join(naics_names, by = c("person_naics" = "code")) 
  
  # calculate proportions within each tract by race
  tract_level <- naics %>%
    group_by(fips_code, person_race) %>%
    summarise(
      tract_total = n(), # total employed in this tract
      race_total = sum(person_naics == occ), # total in selected NAICS for this race in this tract
      .groups = "drop"
    ) %>%
    mutate(proportion = race_total / tract_total * 100) # proportion within the tract
  
  # mean and variability across tracts by race
  summary_data <- tract_level %>%
    group_by(person_race) %>%
    summarise(
      mean_pct = mean(proportion, na.rm = TRUE),
      sd_pct = sd(proportion, na.rm = TRUE), # standard deviation
      se_pct = sd_pct / sqrt(n()), # standard error
      .groups = "drop"
    )
  
  industry_name <- naics_names$name[naics_names$code == occ]

  # plot mean with error bars (standard error)
  ggplot(summary_data, aes(x = person_race, y = mean_pct, fill = person_race)) +
    theme_minimal() +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
      aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
      width = 0.2
    ) +
    scale_x_discrete(labels = c("Asian", "Black", "Multi", "AIAN", "Other", "HIPI", "White")) +
    scale_fill_manual(values = pal1) +
    labs(
      title = paste("% Employment in", industry_name),
      subtitle = "New York",
      x = "Race",
      y = "Mean Percent Employed (± SE)"
    ) +
    ylim(c(0, 12)) +
    guides(fill = FALSE) 
}

c1 <- naics_tracts("622") + 
  theme(plot.title = element_text(size = 8)) + labs(x = NULL, y = NULL)
c2 <- naics_tracts("623") +
  labs(
    x = NULL,
    y = NULL
  ) + 
  theme(plot.title = element_text(size = 8))

c3 <- naics_tracts("722") +
  labs(
    x = NULL,
    y = NULL
  ) + 
  theme(plot.title = element_text(size = 8))

c4 <- naics_tracts("541") +
  labs(
    x = NULL,
    y = NULL
  ) + 
  theme(plot.title = element_text(size = 8))

C <- plot_grid(c1, c2, c3, c4)

# large Figure 2
plot_grid(
  A, B, C,
  nrow = 3,
  labels = c("A", "", "C"),
  rel_heights = c(1, 1.5, 2)
)

