#---------------------------
#interactive plot
#---------------------------

library(pacman)
pacman::p_load(plotly)
pacman::p_load(GGally)
pacman::p_load(survival)
pacman::p_load(cowplot)
pacman::p_load(broom)
pacman::p_load_current_gh("sahirbhatnagar/casebase")
pacman::p_load(Epi)
"%ni%" <- Negate("%in%")

data(lung, package = "survival")
sf_lung <- survival::survfit(survival::Surv(time, status) ~ 1, data = lung)
p1 <- GGally::ggsurv(sf_lung, main = "Kaplan-Meier Curve for the NCCTG Lung Cancer Data")
plotly::ggplotly(p1)

lung <- transform(lung, sex = factor(sex, levels = 1:2, labels = c("Male","Female")))
sf_sex <- survival::survfit(Surv(time, status) ~ sex, data = lung)
pl_sex <- GGally::ggsurv(sf_sex, main = "Kaplan-Meier Curve for the NCCTG Lung Cancer Data Stratified by Sex")
log_rank_sex <- survival::survdiff(Surv(time, status) ~ sex, data = lung)


pl_sex_annotated <- pl_sex + ggplot2::geom_text(aes(label = sprintf("log-rank test p-value: %0.2g", 
                                                                    pchisq(log_rank_sex$chisq, df = 1, lower.tail = F)),
                                                    x = 750, y = 0.9))

plotly::ggplotly(pl_sex_annotated)

data(nickel, package = "Epi")
attach(nickel)
LL <- Lexis.diagram( age=c(10,100), date=c(1900,1990), 
                     entry.age=age1st, exit.age=ageout, birth.date=dob, 
                     fail=(icd %in% c(162,163)), lwd.life=1,
                     cex.fail=0.8, col.fail=c("green","red") )

LL[nickel$icd %in% c(162,163),"cause"] <- "lung"
LL[nickel$icd %in% c(160),"cause"] <- "nasal"
LL[nickel$icd %ni% c(160,162,163), "cause"] <- "other"

lex_plot <- ggplot(LL, aes(x=entry.date, xend=exit.date, y=entry.age, yend=exit.age)) + 
  xlab("Calendar time") +
  ylab("Age") + 
  labs(title = "Lexis Diagram of Nickel Smelting Workers in South Wales")+
  scale_y_continuous(breaks = seq(10,100,10)) +
  scale_x_continuous(breaks = seq(1900,1990,10)) +
  geom_segment(size=.4, colour="grey") +
  geom_point(aes(x = exit.date, y = exit.age, color = cause), 
             data = LL[LL$cause %in% c("lung","nasal"),]) + 
  scale_color_brewer(palette = "Set1") + theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + background_grid(major = "xy", minor = "xy",colour.major = "grey")

ggplotly(lex_plot)
