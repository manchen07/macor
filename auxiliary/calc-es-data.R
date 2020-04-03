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
  mutate(rce_y1y3 = cal_rce(d_y3, n_trt, n_ctl),
         v_rce_y1y3 = cal_v_rce(d_y3, n_trt, n_ctl),
         rce_y1y4 = cal_rce(d_y4, n_trt, n_ctl),
         v_rce_y1y4 = cal_v_rce(d_y4, n_trt, n_ctl))

# Calculate biserial correlation between y2 and y4
es_dat <- 
  es_dat %>%
  mutate(rb_y4_d = dtob(d_y4, n_y2_0, n_y2_1), # d to rb
         rb_y4_t = ttob(t_y4, n_y2_0, n_y2_1), # t to rb
         d_y4_t = cal_t_d(t_y4, n_y2_0, n_y2_1), # d to t
         rb_y4_dt = dtob(d_y4_t, n_y2_0, n_y2_1), # d to t to rb
         v_rb_y4 = cal_v_rb(rb_y4_d, n_y2_0, n_y2_1)) # sampling variance of rb (Soper's) 

# Suggestions
# Sampstat -> d -> r_pb -> r_b using `dtob`  
# t stat -> r_pb -> r_b using `ttob` 

# Calculate tetrachoric correlation between y1 (group) and y2 (dichotomous)

## From contingency table

es_dat_y2 <- 
  subset(es_dat, !is.na(ctl_y2_0)) %>% 
  dplyr::select(studyid, ctl_y2_0, ctl_y2_1, trt_y2_0, trt_y2_1)

rho <- 
  apply(as.matrix(es_dat_y2[,grep("y2",names(es_dat_y2))]), 1, pcor) %>%
  t() %>% as.data.frame() %>%
  mutate(studyid = es_dat_y2$studyid) %>%  
  dplyr::select(studyid, pcor_y1y2 = pcor, v_pcor_y1y2 = v_pcor)

## Final approximation from Divgi (1979)

es_dat <- 
  es_dat %>%
  left_join(rho, by = "studyid") %>% 
  mutate(r_tet_y1y2 = calc_tet(ctl_y2_0, trt_y2_1, ctl_y2_1, trt_y2_0)[, "r_tet"],
         v_tet_y1y2 = calc_tet(ctl_y2_0, trt_y2_1, ctl_y2_1, trt_y2_0)[, "v_tet"])

# or to d to r
es_dat <- 
  es_dat %>% 
  mutate(d_or_y1y2 = log(or) * sqrt(3) / pi, # calculate d from or
         vd_or_y1y2 = 3 * v_or / (pi^2), # variance of d from or
         r_d_or_y1y2 = dtob(d_or_y1y2, n_trt, n_ctl),
         diff = pcor_y1y2 - r_d_or_y1y2) # the difference btw tetrachoric from polychor function and the tetrachoric from OR to d to r

save(es_dat, file = "data/es_dat.RData", compress = TRUE)
