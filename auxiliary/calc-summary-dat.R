library(dplyr)
library(tidyr)

data("raw_dat")

# --------------------------------
# Calculate Summary stats
# --------------------------------

sim_dat <- 
  bind_rows(raw_dat_miss, .id = "column_label")

sim_clean <- 
  sim_dat %>% 
  dplyr::rename(studyid = column_label, y1 = group, y2 = Verbal, y3 = Stroop, y4 = Internalzing) %>%
  dplyr::select(-N) %>% 
  mutate(studyid = as.numeric(studyid),
         y1 = ifelse(y1 == "control", "ctl", "trt")) %>% 
  group_by(studyid) %>%
  mutate(y2_dic = ifelse(y2 <= median(y2), 0, 1),
         y2_dic = factor(y2_dic)) %>% 
  ungroup()

sum_y3y4 <- 
  sim_clean %>%
  pivot_longer(cols = y3:y4,
               names_to = "variables",
               values_to = "values") %>% 
  group_by(studyid, y1, variables) %>% 
  summarise(n = n(),
            m = mean(values),
            sd = sd(values)) %>% 
  pivot_longer(cols = n:sd,
               names_to = "stat",
               values_to = "val")%>% 
  unite(stat, stat, y1, variables, sep = "_") %>% 
  pivot_wider(names_from = stat,
              values_from = val)%>% 
  mutate(n_ctl = n_ctl_y3,
         n_trt = n_trt_y3) %>% 
  dplyr::select(-c(n_ctl_y3, n_ctl_y4, n_trt_y3, n_trt_y4))

sum_y2 <- 
  sim_clean %>%
  group_by(studyid, y1, y2_dic) %>% 
  summarise(n = n()) %>% 
  mutate(y1_y2_dic = ifelse(is.na(y2_dic), NA, paste0(y1, "_y2_", y2_dic, sep = ""))) %>% 
  ungroup() %>% 
  dplyr::select(studyid, y1_y2_dic, n) %>%
  filter(!is.na(y1_y2_dic)) %>% 
  pivot_wider(names_from = y1_y2_dic, values_from = n)

sum_cor <- 
  sim_clean %>%
  group_by(studyid) %>% 
  mutate(cor23 = cor(y2, y3),
         cor24 = cor(y2, y4),
         cor34 = cor(y3, y4)) %>%
  ungroup() %>% 
  group_by(studyid, y1) %>% 
  mutate(cor23_y1 = cor(y2, y3),
         cor24_y1 = cor(y2, y4),
         cor34_y1 = cor(y3, y4)) %>% 
  ungroup() %>%
  dplyr::select(studyid, y1, starts_with("cor")) %>% 
  unique() %>% 
  pivot_longer(cols = cor23_y1:cor34_y1,
               names_to = "cor",
               values_to = "val") %>% 
  unite(cor, cor, y1) %>% 
  pivot_wider(names_from = cor,
              values_from = val) %>%
  rename_all(funs(c("studyid", "cor23", "cor24", "cor34", "cor23_ctl", 
                    "cor24_ctl", "cor34_ctl", "cor23_trt", "cor24_trt", "cor34_trt")))

summary_dat <- 
  sum_cor %>% 
  left_join(sum_y3y4, by = "studyid") %>% 
  left_join(sum_y2, by = "studyid") %>% 
  dplyr::select(studyid, starts_with("n_"), starts_with("cor"), starts_with("m"), 
                starts_with("sd"), starts_with("ctl_y2"), starts_with("trt_y2"))

summary_dat <- 
  summary_dat %>% 
  mutate(n_y2_0 = ctl_y2_0 + trt_y2_0,
         n_y2_1 = ctl_y2_1 + trt_y2_1,
         d_y3 = cal_d(m_trt_y3, sd_trt_y3, n_trt, m_ctl_y3, sd_ctl_y3, n_ctl)[,"d"], v_d_y3 = cal_d(m_trt_y3, sd_trt_y3, n_trt, m_ctl_y3, sd_ctl_y3, n_ctl)[,"v_d"],
         d_y4 = cal_d(m_trt_y4, sd_trt_y4, n_trt, m_ctl_y4, sd_ctl_y4, n_ctl)[,"d"], v_d_y4 = cal_d(m_trt_y4, sd_trt_y4, n_trt, m_ctl_y4, sd_ctl_y4, n_ctl)[,"v_d"],
         t_y3 = cal_t(n_trt, m_trt_y3, sd_trt_y3, m_ctl_y3, sd_ctl_y3, n_ctl),
         t_y4 = cal_t(n_trt, m_trt_y4, sd_trt_y4, m_ctl_y4, sd_ctl_y4, n_ctl),
         or = ctl_y2_0 * trt_y2_1 / (ctl_y2_1 * trt_y2_0),
         v_or = 1/ctl_y2_0 + 1/ctl_y2_1 + 1/trt_y2_0 + 1/trt_y2_1)

# --------------------------------------------
# Introduce missingness at the scenario level
# --------------------------------------------

set.seed(20200401)
scenario <- c("summary", "d", "t")
y13_scenario <- sample(scenario, 100, replace = TRUE)
y14_scenario <- sample(scenario, 100, replace = TRUE)

scenario2 <- c("cellsize", "or")
y12_scenario <- sample(scenario2, 100, replace = TRUE)

scenario3 <- c("summary", "cor")
y24_scenario <- sample(scenario3, 100, replace = TRUE)

summary_dat <-
  summary_dat %>% 
  mutate(y13_scenario = y13_scenario,
         y14_scenario = y14_scenario,
         y12_scenario = y12_scenario,
         y24_scenario = y24_scenario)

summary_dat <- 
  summary_dat %>% 
  
  # regarding y1 and y3
  mutate_at(.vars = vars(m_ctl_y3, m_trt_y3, sd_ctl_y3, sd_trt_y3),
            .funs = list(~ ifelse(y13_scenario %in% c("d", "t"), NA, .))) %>% 
  mutate_at(.vars = vars(d_y3, v_d_y3),
            .funs = list(~ ifelse(y13_scenario == "d", ., NA))) %>% 
  mutate(t_y3 = ifelse(y13_scenario == "t", t_y3, NA)) %>% 
  
  # regarding y1 and y4
  mutate_at(.vars = vars(m_ctl_y4, m_trt_y4, sd_ctl_y4, sd_trt_y4),
            .funs = list(~ ifelse(y14_scenario %in% c("d", "t"), NA, .))) %>% 
  mutate_at(.vars = vars(d_y4, v_d_y4),
            .funs = list(~ ifelse(y14_scenario == "d", ., NA))) %>% 
  mutate(t_y4 = ifelse(y14_scenario == "t", t_y4, NA)) %>% 
  
  # regarding y1 and y2
  mutate_at(.vars = vars(ctl_y2_0, ctl_y2_1, trt_y2_0, trt_y2_1, n_y2_0, n_y2_1),
            .funs = list(~ ifelse(y12_scenario == "cellsize", ., NA))) %>% 
  mutate_at(.vars = vars(or, v_or),
            .funs = list(~ ifelse(y12_scenario == "or", ., NA))) %>% 
  
  # regarding y2 and y4
  mutate(cor24 = ifelse(y24_scenario == "cor", cor24, NA)) # did not include summary bc it might conflicts with y1 and y4

save(summary_dat, file = "data/summary_dat.RData", compress = TRUE)