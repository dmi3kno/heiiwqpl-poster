library(targets)
library(tidyverse)
library(qpd)
dbpoints_plot_data <- tar_read(chocolate_nested)  %>%
  mutate(pq_grd=map(mtlg_a, ~tibble::tibble(p_grd=qpd::make_pgrid(),
                                            q_grd=qpd::qmetalog(p_grd, .x, bl=0)))) %>%
  select(Survey, Ncons, pq_grd) %>%
  unnest(cols = c(pq_grd))

pts <- tar_read(chocolate_nested) %>% unnest(c(p,q)) %>%
  filter(Survey=="(Elderly, 2010)") %>%
  slice(1:5)

dist_df <- tibble(ps = qpd::make_pgrid(),
                .draw=seq_along(ps),
                `Myerson`=with(pts, qpd::qMyerson(ps, q[p==0.05], q[p==0.5], q[p==0.95], alpha = 0.05)),
                `J-QPDS`=with(pts, qpd::qJQPDS(ps, q[p==0.05], q[p==0.5], q[p==0.95], alpha = 0.05)),
                `Metalog`=with(pts, qpd::qmetalog(ps, qpd::fit_metalog(p[c(2:5)],q[c(2:5)], bl=0), bl=0))
                ) %>%
  pivot_longer(-c(ps,.draw), names_to = "Distribution", values_to = "qs")

ggplot(dist_df)+
  geom_line(aes(x=qs, y=ps, color=Distribution), size=0.75)+
  geom_point(data=pts,aes(y=p, x=q), color="black")+
  coord_cartesian(xlim = c(0,25))+
  hrbrthemes::theme_ipsum_rc()+
  scale_y_continuous(breaks = seq(0,1, by=0.1))+
  scale_color_brewer(palette = "Set2")+
  labs(title="Chronic Bitter Chocolate Consumption",
       subtitle="Swedish National Dietary Survey - Riksmaten elderly 2010-11",
       caption="Source: EFSA Food Consumption Database",
       x="Consumption, g/day", y="CDF")

ggsave("qpd_compare.png", width = 10, height = 6, bg="white")
