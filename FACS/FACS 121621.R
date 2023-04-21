library(readxl)
library(tidyverse)
library(plotly)
library(broom)
library(cowplot)

##### Prep #####
## Clean gate to cell crosswalk
xwalk_cell_gate = read_excel("FACS manual clean_postTAC wk1.xlsx", 
                             sheet = 'xwalk') %>% 
  rowwise() %>% 
  mutate(gate = str_sub(gate,
                        str_locate(gate,"LSG")[1],
                        -1L)) 

#x= read_excel("FACS manual clean_CRS.xlsx", 
 #             sheet = "kidney")
## Create function to clean
clean_facs_heart = function(x){
  dfa =  x %>% 
    pivot_longer(-c(sample_name, id)) %>% 
    filter(!is.na(value)) %>% 
    mutate(trt = ifelse(str_detect(name, "TAC"),"TAC","Sham"),
           genotype =  case_when(str_detect(name, "LB2")~"LB2",
                                 str_detect(name, "b2fl")~"b2fl"))%>% 
    rowwise() %>% 
    rename(sample =  id) %>% 
    mutate(gate = str_sub(sample_name,
                          str_locate(sample_name,"LSG")[1],
                          -1L)) %>% 
    ungroup() %>%
    left_join(xwalk_cell_gate) %>% 
    select(sample, trt, genotype, cell,order,value) %>% 
    filter(!is.na(cell)) %>% 
    group_by(sample,trt, genotype) %>% 
    group_modify(~{
      singlet_value_tmp = .x %>% 
        filter(cell == "singlet") %>% 
        pull(value)
      cd_value_tmp =.x %>% 
        filter(cell == "CD45+") %>% 
        pull(value)
      macro_value_tmp = .x %>% 
        filter(cell == "Mac") %>% 
        pull(value)
      myeloid_value_tmp = .x %>% 
        filter(cell == "myeloid") %>% 
        pull(value)
      .x %>%
        mutate(singlet = singlet_value_tmp, 
               cd45 = cd_value_tmp,
               mac = macro_value_tmp,
               myeloid = myeloid_value_tmp
        ) 
    }) %>% 
    ungroup() %>% 
    mutate(pct_singlet = round(value/singlet*100,2),
           pct_cd45 = round(value/cd45*100,2),
           pct_mac = round(value/mac*100,2),
           pct_myeloid = round(value/myeloid*100,2)) %>% 
    filter(cell!="singlet")
  
  cell_order = dfa %>% 
    select(cell, order) %>% 
    distinct() %>% 
    arrange(order) %>% 
    pull(cell)
  
  dfa %>% 
    mutate(cell = factor(cell,levels = rev(cell_order))) %>% 
    mutate(genotype = factor(genotype, levels = c("b2fl","LB2")))%>%
    mutate(trt = factor(trt, levels = c("Sham","TAC")))
  
}


##function 2 updated 092121###
plot_pct_singlet2 = function(df1, df1_avg){
  df_avg_tmp =  df1 %>% 
    add_count(trt, genotype, cell) %>% 
    group_by(trt, genotype, cell) %>% 
    summarise(se = sd(value_to_plot, na.rm = T)/sqrt(n),
              value_to_plot = mean(value_to_plot)) %>% 
    ungroup() %>% 
    distinct()%>% 
    mutate(trt = factor(trt, levels = c("Sham","TAC")))
  
  ggplot(data = df_avg_tmp,
         aes(x=genotype, y = value_to_plot, fill = trt))+
    geom_bar(data = df_avg_tmp ,
             aes(y=value_to_plot, x = genotype, fill = trt),
             alpha = 0.5,
             color = 'black',
             width = 0.6,
             size = 1,
             position=position_dodge(width = 0.7),stat="identity")+
    geom_point(data = df1 %>% 
                 mutate(trt = factor(trt, levels = c("Sham","TAC"))),
               aes(x = genotype, y=value_to_plot, 
                   fill = trt, 
                   group = trt),
               position=position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.7),
               shape = 21,
               size = 2.5)+
    geom_errorbar(data = df_avg_tmp,
                  aes(ymin = value_to_plot-se, ymax = value_to_plot+se), 
                  width=.1,
                  position=position_dodge(.7))+
    #facet_wrap(~trt)+
    labs(x = "", y= "%")+
    scale_fill_manual(values = c('white','black'))+
    scale_y_continuous(expand = c(0,0), limits = c(0,max(df1$value_to_plot)*1.1)) +
    theme_prism()
  
}


##### Clean Data #####
df_heart = read_excel("FACS manual clean_postTAC wk1.xlsx", sheet="heart") %>% 
  clean_facs_heart 





######  1. heart Figures  ######
df_heart %>% 
  filter(!cell=="CD45+") %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2() +
  labs(y = "% of Singlets")+
  facet_wrap(~cell)

df_heart %>% 
  filter(cell%in%c("CCR2+mac","CCR2-mac")) %>% 
  mutate(value_to_plot = pct_mac) %>% 
  plot_pct_singlet2()+
  labs(y="% of mac")+
  facet_wrap(~cell)

df_myeloid= df_kidney %>% 
  filter(cell=="myeloid") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_neutrophil= df_kidney %>% 
  filter(cell=="neutrophil") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_macrophage= df_kidney %>% 
  filter(cell=="Mac") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_M2_singlet = df_kidney %>% 
  filter(cell=="M2") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_M2_mac = df_kidney %>% 
  filter(cell=="M2") %>% 
  select(sample, genotype, cell, pct_mac) %>% 
  pivot_wider(names_from = sample, values_from = pct_mac)

df_mono = df_kidney %>% 
  filter(cell=="inf.mono") %>% 
  select(sample, genotype, cell, pct_singlet) %>%
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_kidney_prism = bind_rows(df_CD45,
                           df_myeloid,
                           df_neutrophil,
                           df_macrophage,
                           df_M2_singlet,
                           df_M2_mac,
                           df_mono)
write.csv(df_kidney_prism, file= "kidney FACS prism.csv")

df_kidney %>% 
  filter(cell=="myeloid") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD45+CD11b+")+
  facet_wrap(~sex)

ggsave(filename = "kidney_myeloid.png")

df_kidney %>% 
  filter(cell=="DCs") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+CD11c")

ggsave(filename= "kidney_DCs.png")

df_kidney %>% 
  filter(cell=="Mac") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+F4_80+")+
  facet_wrap(~sex)

ggsave(filename = "kidney_mac_CD11b included.png")
 
df_kidney %>% 
  filter(cell=="neutrophil") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+Ly6G+")+
  facet_wrap(~sex)

ggsave(filename = "kidney_neutrophil.png")

df_kidney %>% 
  filter(cell=="M2") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets") +
  ggtitle("CD206+mac")

df_kidney %>% 
  filter(cell=="M2") %>% 
  mutate(value_to_plot = pct_mac) %>% 
  plot_pct_singlet2()+
  labs(y="% of Mac")+
  ggtitle("CD206+mac")

df_kidney %>% 
  filter(cell=="inf.mono") %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("Ly6C hi CD11b+")+
  facet_wrap(~sex)

####2. heart figures####

df_CD45= df_heart %>% 
  filter(cell=="CD45+") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_myeloid= df_heart %>% 
  filter(cell=="myeloid") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_neutrophil= df_heart %>% 
  filter(cell=="neutrophil") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_macrophage= df_heart %>% 
  filter(cell=="Mac") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_M2_singlet = df_heart %>% 
  filter(cell=="M2") %>% 
  select(sample, genotype, cell, pct_singlet) %>% 
  pivot_wider(names_from = sample, values_from = pct_singlet)

df_M2_mac = df_heart %>% 
  filter(cell=="M2") %>% 
  select(sample, genotype, cell, pct_mac) %>% 
  pivot_wider(names_from = sample, values_from = pct_mac)

df_mono = df_heart %>% 
  filter(cell=="inf.mono") %>% 
  select(sample, genotype, cell, pct_singlet) %>%
  pivot_wider(names_from = sample, values_from = pct_singlet)

write.csv(df_mono, file="heart mono.csv")

df_heart_prism = bind_rows(df_CD45,
                           df_myeloid,
                           df_neutrophil,
                           df_macrophage,
                           df_M2_singlet,
                           df_M2_mac,
                           df_mono)

write.csv(df_heart_prism, file="heart FACS prism.csv")

#########################
df_heart %>% 
  filter(cell == "CD45+") %>%
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2() +
  labs(y = "% of Singlets")+
  ggtitle("CD45+")


df_heart %>% 
  filter(cell=="myeloid") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD45+CD11b+")

ggsave(filename = "kidney_myeloid.png")

df_heart %>% 
  filter(cell=="Mac") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+F4_80+")

ggsave(filename = "kidney_mac_CD11b included.png")

df_heart %>% 
  filter(cell=="neutrophil") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+Ly6G+")

ggsave(filename = "kidney_neutrophil.png")

df_heart %>% 
  filter(cell=="M2") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets") +
  ggtitle("CD206+mac")

df_heart %>% 
  filter(cell=="M2") %>% 
  mutate(value_to_plot = pct_mac) %>% 
  plot_pct_singlet2()+
  labs(y="% of Mac")+
  ggtitle("CD206+mac")

df_heart %>% 
  filter(cell=="inf.mono") %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("Ly6C hi CD11b+")


#####3. Spleen figures#####
df_spleen %>% 
  filter(cell == "CD45+") %>%
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2() +
  labs(y = "% of Singlets")+
  ggtitle("CD45+")


df_spleen %>% 
  filter(cell=="myeloid") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD45+CD11b+")

df_spleen %>% 
  filter(cell=="Mac") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+F4_80+")

df_spleen %>% 
  filter(cell=="inf.mono") %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("Ly6C hi CD11b+")

df_spleen %>% 
  filter(cell=="neutrophil") %>% 
  mutate(value_to_plot=pct_singlet) %>% 
  plot_pct_singlet2()+
  labs(y="% of Singlets")+
  ggtitle("CD11b+Ly6G+")


######  1. T-test for mLN (value_to_plot) ######
## Heart genotype Results chow vs HFD
results_mLN_genotype = df_mLN %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  arrange(trt, genotype, cell) %>% 
  group_by(trt, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(genotype) %>% unique()
    values_grp1 = .x %>% filter(genotype==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(genotype==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_mLN_genotype,file = "results_mLN_genotype_pct_singlet.csv")
results_mLN_genotype_pct_myeloid = df_mLN %>% 
  filter(cell%in%c("DCs","Mac","neutrophils")) %>% 
  mutate(value_to_plot = pct_myeloid) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  arrange(trt, genotype, cell) %>% 
  group_by(trt, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(genotype) %>% unique()
    values_grp1 = .x %>% filter(genotype==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(genotype==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_mLN_genotype_pct_myeloid,file = "results_mLN_genotype_pct_myeloid.csv")
## Heart Treatment Results MI_d7 vs MI_d25
results_mLN_trt = df_mLN %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  group_by(genotype, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(trt) %>% unique()
    values_grp1 = .x %>% filter(trt==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(trt==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_mLN_trt, file = "results_mLN_trt_pct_singlet.csv")
results_mLN_trt_pct_myeloid = df_mLN %>% 
  mutate(value_to_plot = pct_myeloid) %>% 
  filter(cell%in%c("DCs","Mac","neutrophils")) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  group_by(genotype, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(trt) %>% unique()
    values_grp1 = .x %>% filter(trt==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(trt==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_mLN_trt_pct_myeloid, file = "results_mLN_trt_pct_myeloid.csv")

######  2. T-test for Heart (value_to_plot)  ######
## Heart genotype Results chow vs HFD
results_heart_genotype = df_heart %>%
  mutate(value_to_plot  = pct_singlet) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  arrange(trt, genotype, cell) %>% 
  group_by(trt, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(genotype) %>% unique()
    values_grp1 = .x %>% filter(genotype==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(genotype==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_heart_genotype,file = "results_heart_genotype_pct_singlet.csv")
## Heart Treatment Results MI_d7 vs MI_d25
results_heart_trt = df_heart %>% 
  mutate(value_to_plot = pct_singlet) %>% 
  select(sample, trt, genotype, cell, value_to_plot) %>% 
  group_by(genotype, cell) %>% 
  group_modify(~{
    groups_tmp = .x %>% pull(trt) %>% unique()
    values_grp1 = .x %>% filter(trt==groups_tmp[1]) %>% pull(value_to_plot)
    values_grp2 = .x %>% filter(trt==groups_tmp[2]) %>% pull(value_to_plot)
    model = t.test(values_grp1,values_grp2)
    p_value_tmp = tidy(model) %>% pull(p.value)
    .x %>% 
      mutate(p_value = p_value_tmp)
  }) %>% 
  ungroup() %>%
  select(-sample, -value_to_plot) %>% 
  distinct() %>% 
  arrange(p_value)
write.csv(results_heart_trt,file = "results_heart_trt_pct_singlet.csv")


