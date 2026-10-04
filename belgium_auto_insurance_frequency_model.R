library(Hmisc)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(forcats)
library(broom)
library(mgcv)      
library(sf)        
library(classInt)  


ds_o <- read.xlsx("mtpl_belgium_291024.xlsx")

#basictransformation
ds <- ds_o %>% 
  mutate_if(is.character, as.factor) %>%
  mutate(sex = factor(sex, levels = c("male", "female"))) %>%
  mutate(fuel = factor(fuel, levels = c('gasoline', 'diesel'))) %>%
  mutate(power = ifelse(power >= 111, 111, power)) %>%
  mutate(power = ifelse(power < 20, 20, power)) %>%
  mutate(power = cut(power, seq(10, 120, 10), right = FALSE, include.lowest = TRUE)) %>%
  mutate(bm = ifelse(bm > 11, 12, bm)) %>%
  mutate(bm = as.factor(bm)) %>%
  mutate(ageph_q20 = cut(ageph, breaks = quantile(ageph, probs = seq(0, 1, by = 0.05), na.rm = TRUE), include.lowest = TRUE, right = FALSE))


ds_new <- ds %>%
  mutate(
    coverage_group = fct_collapse(coverage,
                                  "TPL_Extended" = c("TPL+", "TPL++") 
    ),
    power_group = fct_collapse(power,
                               "Very_Low"     = c("[20,30)"),
                               "Low"          = c("[30,40)"),
                               "Middle_Base"  = c("[40,50)", "[50,60)", "[60,70)"),
                               "High"         = c("[70,80)", "[80,90)", "[90,100)", "[100,110)"),
                               "Very_High"    = c("[110,120]")
    ),
    age_group = fct_collapse(ageph_q20,
                             "Young_Risk"      = c("[18,25)"),
                             "Young_Moderate"  = c("[25,28)", "[28,31)"), 
                             "Adult_Base"      = c("[31,33)", "[33,35)", "[35,37)", "[37,39)",
                                                   "[39,41)", "[41,43)", "[43,46)", "[46,48)", "[48,50)"),
                             "Pre_Senior"      = c("[50,52)", "[52,55)", "[55,58)"),
                             "Senior"          = c("[58,61)", "[61,65)", "[65,68)", "[68,73)"),
                             "Elderly"         = c("[73,95]")
    ),
    bm_group = fct_collapse(bm,
                            "0"      = c("0"),
                            "1"      = c("1"),
                            "2"      = c("2"),
                            "3"      = c("3"),
                            "4-5"    = c("4", "5"),
                            "6-8"    = c("6", "7", "8"),
                            "9-10"   = c("9", "10"),
                            "BM_Top" = c("11", "12")
    )
  ) %>%
   #setting reference levels
  mutate(coverage_group = fct_relevel(coverage_group, "TPL")) %>%
  mutate(power_group = fct_relevel(power_group, "Middle_Base")) %>%
  mutate(age_group = fct_relevel(age_group, "Adult_Base")) %>%
  mutate(bm_group = fct_relevel(bm_group, "0"))

#final GLM model 
final_glm <- glm(nclaims ~ coverage_group + fuel + age_group + bm_group + power_group,
                 offset = log(expo),
                 family = poisson(link = 'log'),
                 data = ds_new)

summary(final_glm)
AIC(final_glm)

confint.default(final_glm)

ds_new %>% 
  group_by(coverage_group) %>%
  summarize(emp_freq = sum(nclaims)/sum(expo),expo=sum(expo)) %>% 
  arrange(desc(expo)) %>% 
  as.data.frame()

ds_new %>% 
  group_by(fuel) %>%
  summarize(emp_freq = sum(nclaims)/sum(expo),expo=sum(expo)) %>% 
  arrange(desc(expo)) %>% 
  as.data.frame()

ds_new %>% 
  group_by(age_group) %>%
  summarize(emp_freq = sum(nclaims)/sum(expo),expo=sum(expo)) %>% 
  arrange(desc(expo)) %>% 
  as.data.frame()

ds_new %>% 
  group_by(bm_group) %>%
  summarize(emp_freq = sum(nclaims)/sum(expo),expo=sum(expo)) %>% 
  arrange(desc(expo)) %>% 
  as.data.frame()

ds_new %>%
  group_by(power_group) %>%
  summarize(total_expo = sum(expo)) %>%
  arrange(desc(total_expo)) %>% 
  as.data.frame() 


# GAM model
model_gam <- gam(nclaims ~ coverage_group + fuel + age_group + bm_group + power_group + s(long, lat, bs = "tp"),
                 offset = log(expo),
                 family = poisson(link = "log"),
                 data = ds_new)

summary(model_gam)
AIC(model_gam)


#extracting the influence of location
pred_spatial <- predict(model_gam, type = "terms", terms = "s(long,lat)")

#dividing into 5 zones using the Fisher method
classint_fisher <- classIntervals(pred_spatial, 5, style = "fisher")

#adding a zone to the dataset
ds_new$geo_zone <- cut(pred_spatial, 
                       breaks = classint_fisher$brks, 
                       right = FALSE, 
                       include.lowest = TRUE, 
                       labels = c("Zone_1", "Zone_2", "Zone_3", "Zone_4", "Zone_5"))

#reduction to 4 zones
ds_new$geo_zone_4groups <- fct_collapse(ds_new$geo_zone,
                                        "Zone_1"   = "Zone_1",
                                        "Zone_2"   = "Zone_2",
                                        "Zone_3_4" = c("Zone_3", "Zone_4"), 
                                        "Zone_5"   = "Zone_5")

#setting the reference level
ds_new$geo_zone_4groups <- fct_relevel(ds_new$geo_zone_4groups, "Zone_3_4")

#
final_glm_spatial <- glm(nclaims ~ coverage_group + fuel + age_group + bm_group + power_group + geo_zone_4groups,
                         offset = log(expo),
                         family = poisson(link = "log"),
                         data = ds_new)

summary(final_glm_spatial)
AIC(final_glm_spatial)


#load map
belgium_shape_sf <- st_read('./4326/postaldistricts.shp')

#preparing data for prediction
post_ds <- belgium_shape_sf
post_ds <- st_centroid(post_ds)
post_ds$long <- do.call(rbind, post_ds$geometry)[,1]
post_ds$lat  <- do.call(rbind, post_ds$geometry)[,2]

post_ds$fuel           <- factor(levels(ds_new$fuel)[1],           levels = levels(ds_new$fuel))
post_ds$coverage_group <- factor(levels(ds_new$coverage_group)[1], levels = levels(ds_new$coverage_group))
post_ds$age_group      <- factor(levels(ds_new$age_group)[1],      levels = levels(ds_new$age_group))
post_ds$bm_group       <- factor(levels(ds_new$bm_group)[1],       levels = levels(ds_new$bm_group))
post_ds$power_group    <- factor(levels(ds_new$power_group)[1],    levels = levels(ds_new$power_group))
post_ds$expo           <- 1

#risk map prediction
pred_map <- predict(model_gam, newdata = post_ds, type = "terms", terms = "s(long,lat)")

#connecting to the map
ds_pred <- data.frame(
  nouveau_PO = post_ds$nouveau_PO, 
  fit_spatial = pred_map[,1]
)
ds_pred_unique <- ds_pred %>% distinct(nouveau_PO, .keep_all = TRUE)
belgium_shape_map <- left_join(belgium_shape_sf, ds_pred_unique, by = "nouveau_PO")

ggplot(belgium_shape_map) +
  geom_sf(aes(fill = fit_spatial), colour = NA) +
  ggtitle("Spatial Risk Map") +
  scale_fill_gradient(low = "#99CCFF", high = "#003366", name = "Risk (Log)") +
  theme_bw()

#map with tariff zones
classint_fisher_map <- classIntervals(ds_pred_unique$fit_spatial, n = 4, style = "fisher")

belgium_shape_map$geo_zone_map <- cut(belgium_shape_map$fit_spatial, 
                                      breaks = classint_fisher_map$brks, 
                                      include.lowest = TRUE,
                                      labels = c("Zone 1", "Zone 2", "Zone 3 (Base)", "Zone 4"))

ggplot(belgium_shape_map) +
  geom_sf(aes(fill = geo_zone_map), color = "white", size = 0.05) +
  scale_fill_brewer(palette = "OrRd", name = "Tariff zone") +
  ggtitle("Risk map (4 zonses)") +
  theme_bw()

