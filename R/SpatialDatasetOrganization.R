library(sf)
library(readr)
library(dplyr)

# csv
# Those results were sent to GoogleEathEngine
# all data
all.pts <- read_csv("../Data/CSV/Database_18_06_01.csv") %>%
  mutate(Longitude = as.numeric(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, na.fail = FALSE)
plot(st_geometry(all.pts))
write_sf(all.pts, "../Data/Vector/Database_10_06_01.shp")

# Flora ----
# data sent by e-mail on "ter, 29 de out 12:49" by Pablo

# Richness
richness <- read_csv("../Data/CSV/Dados.RF.Plant.richness.csv") %>% 
  mutate(GDef50 = GrossDef_50Km,
         phox5 = phihox_05Km) %>% select(Site, RR_var, GDef50, phox5, Longitude, Latitude)
richness.sf <- st_as_sf(richness, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(richness.sf))
write_sf(richness.sf, "../Data/Vector/FloraRichness.shp")

# Abundance
abundance <- read_csv("../Data/CSV/Dados.RF.Plant.abundance.csv") %>% 
  mutate(f03_100Km = f2003_100Km) %>% select(Site, RR_var, f03_100Km, sa_75Km, Longitude, Latitude) 
abundance.sf <- st_as_sf(abundance, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(abundance.sf))
write_sf(abundance.sf, "../Data/Vector/FloraAbundance.shp")


# Invertebrados ----
# data sent by e-mail on "ter, 29 de out 12:49" by Pablo

# Richness
InvertRichness <- read_csv("../Data/CSV/Dados.RF.Invert.richness.csv") %>% 
  select(-1)
colnames(InvertRichness)
InvertRichness.sf <- st_as_sf(InvertRichness, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(InvertRichness.sf))
write_sf(InvertRichness.sf, "../Data/Vector/InvertebradosRichness.shp")

# Abundance
InvertAbundance <- read_csv("../Data/CSV/Dados.RF.Invert.abundance.csv") %>% select(-1) 
colnames(InvertAbundance)
InvertAbundance.sf <- st_as_sf(InvertAbundance, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(InvertAbundance.sf))
write_sf(InvertAbundance.sf, "../Data/Vector/InvertAbundance.shp")

# Vertebrados ----
# data sent by e-mail on "ter, 29 de out 12:49" by Pablo

# Richness
VertRichness <- read_csv("../Data/CSV/Dados.RF.Vert.richness.csv") %>% 
  select(-1)
colnames(VertRichness)
VertRichness.sf <- st_as_sf(VertRichness, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(VertRichness.sf))
write_sf(VertRichness.sf, "../Data/Vector/VertebradosRichness.shp")

# Abundance
VertAbundance <- read_csv("../Data/CSV/Dados.RF.Vert.abundance.csv") %>% select(-1) 
colnames(VertAbundance)
VertAbundance.sf <- st_as_sf(VertAbundance, coords = c("Longitude", "Latitude"), crs = 4326)
plot(st_geometry(VertAbundance.sf))
write_sf(VertAbundance.sf, "../Data/Vector/VertAbundance.shp")
