library(readr)
library(dplyr)

bio_revisionNew <- read_csv("../Data/CSV/Biodiversity_Revision_09_06_2020.csv")
bio_revision <- read_csv("../Data/CSV/Biodiversity_Revision.csv") %>% 
  select(Site)

bio_revisionNew %>% right_join(bio_revision) %>% 
  write_csv("../Data/CSV/Biodiversity_Revision_09_06_2020_filtrado.csv")


# analise do resultado - confirmando summary
bio_revisionNew %>% right_join(bio_revision) %>% select( starts_with("percUrbArea")) %>% max()
