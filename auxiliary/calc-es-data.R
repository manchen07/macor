library(dplyr)
library(tidyr)
library(polycor)

data("summary_dat")

# --------------------------------
# Get correlation estimates
# --------------------------------

# Calculate correlation between y1 and y3 or y4
es_dat <- 
  summary_dat %>%
  dplyr::select(-ends_with("scenario")) %>% 
  mutate(cor_y1y3 = cal_rce(d_y3, n_trt, n_ctl),
         v_cor_y1y3 = cal_v_rce(d_y3, n_trt, n_ctl),
         cor_y1y4 = cal_rce(d_y4, n_trt, n_ctl),
         v_cor_y1y4 = cal_v_rce(d_y4, n_trt, n_ctl))

# Calculate biserial correlation between y2 and y4
es_dat <- 
  es_dat %>%
  mutate(d2_y4 = (m_1_y4 - m_0_y4)/sqrt(((n_1_y4-1)*sd_1_y4^2 + (n_0_y4-1)*sd_0_y4^2)/(n_1_y4 + n_0_y4 -2)),
         cor_y2y4 = dtob(d2_y4, n_0_y4, n_1_y4), # d to rb
         v_cor_y2y4 = cal_v_rb(cor_y2y4, n_0_y4, n_1_y4)) # sampling variance of rb (Soper's) 

# Calculate tetrachoric correlation between y1 (group) and y2 (dichotomous)

## From contingency table

es_dat_y2 <- 
  subset(es_dat, !is.na(ctl_y2_0)) %>% 
  dplyr::select(studyid, ctl_y2_0, ctl_y2_1, trt_y2_0, trt_y2_1)

rho <- 
  apply(as.matrix(es_dat_y2[,grep("y2",names(es_dat_y2))]), 1, pcor) %>%
  t() %>% as.data.frame() %>%
  mutate(studyid = es_dat_y2$studyid) %>%  
  dplyr::select(studyid, cor_y1y2 = pcor, v_cor_y1y2 = v_pcor)

es_dat <- 
  es_dat %>%
  left_join(rho, by = "studyid")

save(es_dat, file = "data/es_dat.RData", compress = TRUE)
