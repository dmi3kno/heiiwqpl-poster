tar_target(chocolate_nested, {
  chocolate_raw <- readxl::read_excel("data/Foodex 2 L5 dashboard.xlsx", skip = 2) 
  chocolate <- chocolate_raw %>% 
    filter(`Survey start year`>2000) %>% 
    mutate(Survey=paste0("(", `Population Group`,", ",`Survey start year`,")")) %>% 
    select(21,10,14:19) %>% 
    rename("50th percentile"="Median",
           "Ncons"="Number of consumers")
  
  chocolate %>% 
    pivot_longer(-c(Survey,Ncons), names_to = "p", values_to = "q") %>% 
    mutate(p=parse_number(p)/100) %>% 
    group_by(Survey,Ncons,q) %>% 
    summarize(p=mean(p)) %>% 
    group_by(Survey, Ncons) %>% 
    summarise(p=list(p), q=list(q)) %>% 
    mutate(mtlg_a=map2(p, q, ~qpd::approx_max_metalog(q=.y, p_grid =.x, bl=0))) %>% 
    mutate(mtlg_v=map_lgl(mtlg_a, qpd::is_metalog_valid, bl=0))
})
