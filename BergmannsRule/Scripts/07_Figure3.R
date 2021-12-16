##############################################################
# Authors: 
# Erin Henry, Ana Benitez-Lopez (@anabenlop)
# Email: erinhenry55@gmail.com, abenitez81@gmail.com

# Script to create figure 3 in the Bergmann's rule paper
# Made 18 October 2020
# Adapted on 16 December 2021

# Figure 3: Meta-regressions (Heat Balance and Migration) ----------------------

### Migration meta-regression panel
# get migration models
mig.mods <- readRDS("Results/Bergmannsrule_results_MR_mig.rds")

# Get sample sizes for figure legend
birds <- readRDS("Results/BergmannsRule_results_correlationsBirdsMR_20211115.rds")
nrow(subset(birds,migration=="Resident"))/7 # N resident species
nrow(subset(birds,migration=="Migratory"))/7 # N migratory species

names(mig.mods)

# dataframe of models
model <- data.frame(beta = unlist(lapply(mig.mods,"[","beta")),
                    env.var = rep(c("MT","MinT","MaxT","MP","PET","NPP","NPPsd"),each=2),
                    ci.lb = unlist(lapply(mig.mods,"[","ci.lb")),
                    ci.ub = unlist(lapply(mig.mods,"[","ci.ub")),
                    n = unlist(sapply(mig.mods,"[","k")),
                    migration = rep(c("Migratory","Resident"),times=7),
                    ID = 1:14)
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)

model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")
model$env.var <- factor(model$env.var,levels=c("MT","MaxT","NPP","NPPsd"))

# make plot
p <- ggplot(model,aes(x=beta,y=migration),show.legend=T) +
  geom_vline(xintercept = 0, color = "gray80",size=.5) +
  geom_point(mapping=aes(x=beta,y=migration,color=env.var),show.legend=T,
             data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=F,size=.75, 
                 mapping=aes(x=beta,y=migration,color=env.var,xmin=ci.lb,
                             xmax=ci.ub,height=0)) +
  theme_classic2() +
  labs(x="Spearman's r",y=NULL) +
  scale_color_manual(name="Variable",
                     values=c("#4477AA","#EE6677","#228833","#AA3377")) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=8),
        plot.margin=unit(c(1,5,1,1),"lines"),
        panel.border=element_rect(fill=NA),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9))

p <- p + facet_wrap(~env.var,nrow=4) + 
  theme(strip.text.x = element_blank(), # remove facet labels
        strip.background = element_blank())

p



### Heat balance meta-regression (Panel b)
hb.mod <- readRDS("Results/BergmannsRule_results_MR_heatBalance.rds")

# Get sample sizes for figure legend
hb <- readRDS("Results/BergmannsRule_results_correlations_20211114.rds")
nrow(subset(hb,class=="reptile" |order=="Anura"))/7 # N thermoregulators
nrow(subset(hb,order=="Caudata"))/7 # N thermoconformers
rm(hb)

# dataframe of models
model <- data.frame(beta = hb.mod$beta,
                    env.var <- c("MT","MT"),
                    ci.lb = hb.mod$ci.lb,
                    ci.ub = hb.mod$ci.ub,
                    therm = c("TC","TR"))
model$beta <- transf.ztor(model$beta)
model$ci.lb <- transf.ztor(model$ci.lb)
model$ci.ub <- transf.ztor(model$ci.ub)

model <- subset(model,env.var!="MP" & env.var!="MinT" & env.var!="PET")
model$therm <- factor(model$therm,levels=c("TR","TC"))

# make plot
p2 <- ggplot(model,aes(x=beta,y=therm),show.legend=F) +
  geom_vline(xintercept = 0, color = "gray80",size=.55) +
  geom_point(mapping=aes(x=beta,y=therm,color=env.var),show.legend=F,
             data=model,size=2) + 
  geom_errorbarh(data=model,show.legend=F,size=.75, 
                 mapping=aes(x=beta,y=therm,color=env.var,xmin=ci.lb,
                             xmax=ci.ub,height=0)) +
  theme_classic2() +
  labs(x="Spearman's r",y=NULL) +
  scale_color_manual(name="Variable",
                     values=c("#4477AA")) +
  theme(axis.title=element_text(size=9),
        plot.title = element_text(size=12),
        axis.text=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        plot.margin=unit(c(1,5,1,1),"lines"),
        panel.border=element_rect(fill=NA),
        legend.spacing.x = unit(.1,'cm'),
  )
p2

### Combine (Into Fig 3.)
ggarrange(p2,p,ncol=2,nrow=1,
          labels="auto",#hjust=-5,vjust=2,
          common.legend = T,
          legend = "bottom") # save as 1000 x 400

### Save figure
ggsave(filename='Figures/Figure_3.png', 
       width=180, height=80, units = 'mm', dpi=600)
