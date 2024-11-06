sample_table_histSA_c4$Location <- factor(sample_table_histSA_c4$Location, levels=c("Concepcion (South Chile)", "La Serena",
                                                        "San Carlos de Bariloche",
                                                        "Mendoza", "Rio Negro", "Hilario Ascasubi", "Colonia (Uruguay)", "Porto Alegre"))
sample_table_histSA_c4$historic_E_W <- factor(sample_table_histSA_c4$historic_E_W, levels=c("West", "Intermediate", "East"))

ggplot(sample_table_histSA_c4, aes(Location, MhFV_meandepth, color=historic_E_W))+
  geom_boxplot()+
  geom_point(alpha=0.7)+
  scale_colour_viridis(discrete=T)+
  theme_bw()

ggplot(sample_table_histSA_c4, aes(Location, MhFV_coverage, color=historic_E_W))+
  geom_boxplot()+
  geom_point(alpha=0.7)+
  scale_colour_viridis(discrete=T)+
  theme_bw()
