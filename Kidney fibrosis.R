library(readxl)
library(tidyverse)
library(janitor)

#### 1. Data Management ####
xwalk=read_excel("CRS_postTAC wk6_G1-G3 fibrosis 020521.xlsx", sheet = "xwalk") %>% 
  bind_rows(read_excel("CRS postTAC wk12_G2n4_fibrosis.xlsx", sheet = "xwalk")) %>% 
  mutate(id = as.character(id))

## Clean Data
df_6 = read_excel("CRS_postTAC wk6_G1-G3 fibrosis 020521.xlsx", sheet = "cortex_data")%>%
  clean_names() %>% 
  mutate(time="6wks")

df_12 = read_excel("CRS postTAC wk12_G2n4_fibrosis.xlsx", sheet = "cortex_data") %>% 
  clean_names() %>% 
  mutate(time="12wks")

df= df_6 %>% 
  bind_rows(df_12) %>% 
  select(-min_thr,-max_thr) %>% 
  #mutate(location = ifelse(str_detect(label,"cortex"),
                           #"Cortex",
                           #"Medulla")) %>% 
  mutate(id = str_sub(label, 1,3)) %>% 
  left_join(xwalk) %>% 
  select(time, trt, diet, id, pct_area) %>% 
  filter(!is.na(trt)) %>% 
  group_by(trt, diet, id, time) %>%
  summarize(pct_area = mean(pct_area)) %>% 
  ungroup()  %>% 
  mutate(time= factor(time, levels= c("6wks","12wks")))%>% 
  arrange(trt, diet, time)



## Plot Cortex
df_avg =  df %>% 
 # filter(location=="Cortex") %>% 
  group_by(trt, diet, time) %>% 
  summarize(pct_area = mean(pct_area)) %>% 
  ungroup() %>% 
  mutate(time= factor(time, levels= c("6wks","12wks")))

ggplot()+
  geom_bar(data = df_avg ,
           aes(x=time, y = pct_area, fill = trt),
           alpha = 0.5,
           position="dodge",stat="identity")+
  geom_point(data = df,
             aes(x = time, y= pct_area, col = trt, group = trt),
             position=position_jitterdodge(jitter.width = 0.1),
             size = 4)+
  theme_bw()+
  theme(strip.text = element_text(size = 18),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_blank())+
  labs(y="% Fibrosis")+
 # scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~diet)

df_chow = df %>% 
  filter(diet=="chow") %>% 
  group_by(trt, time) %>% 
  mutate(trt2 = paste0(trt, row_number())) %>% 
  ungroup() %>% 
  select(-diet,-trt,-id) %>% 
  pivot_wider(names_from = trt2, values_from = pct_area) 

 #mutate(week = parse_number(as.character(time))) %>% 
 # select(-time) %>% 
  #select(week, everything())

df_HFD = df %>% 
  filter(diet=="HFD") %>% 
  group_by(trt, time) %>% 
  mutate(trt2 = paste0(trt, row_number())) %>% 
  ungroup() %>% 
  select(-diet,-trt,-id) %>% 
  pivot_wider(names_from = trt2, values_from = pct_area) 

write.csv(df_HFD,"results/HFD_kidney fibrosis.csv")
  

