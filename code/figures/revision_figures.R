suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(patchwork); library(svglite); library(ragg); library(grid)
})
options(stringsAsFactors=FALSE, scipen=999)

args <- commandArgs(trailingOnly=TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash="/", mustWork=TRUE) else normalizePath(".", winslash="/", mustWork=TRUE)
root <- file.path(repo_root,"data","derived")
render_root <- file.path(repo_root,"results")
out_dir <- file.path(render_root,"recreated_figures")
source_dir <- file.path(render_root,"recreated_figure_source_data")
qa_dir <- Sys.getenv("MCD_QA_DIR", unset=file.path(render_root,"qa"))
dir.create(qa_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)

pal <- c(ink="#222222",gray="#6F6F6F",light="#D9D9D9",pale="#F4F5F7",
         blue="#0072B2",sky="#56B4E9",green="#009E73",orange="#E69F00",
         vermillion="#D55E00",purple="#CC79A7",yellow="#F0E442")
dest_order <- c("CIRCULATORY","NEOPLASMS","INFECTIOUS_PARASITIC","RESPIRATORY",
                "DIGESTIVE_OTHER","DIABETES_NUTRITIONAL","EXTERNAL",
                "OTHER_GALLBLADDER_BILIARY","OTHER")
dest_cols <- c(CIRCULATORY="#0072B2",NEOPLASMS="#D55E00",INFECTIOUS_PARASITIC="#CC79A7",
               RESPIRATORY="#56B4E9",DIGESTIVE_OTHER="#009E73",DIABETES_NUTRITIONAL="#E69F00",
               EXTERNAL="#6C5B7B",OTHER_GALLBLADDER_BILIARY="#B7B7B7",OTHER="#525252")
disease_order <- c("K80","I21","C18","E10_E14","N18","A40_A41")
disease_labels <- c(K80="Cholelithiasis (K80)",I21="Acute myocardial infarction (I21)",
                    C18="Colon cancer (C18)",E10_E14="Diabetes (E10-E14)",
                    N18="Chronic kidney disease (N18)",A40_A41="Sepsis (A40-A41)")
disease_cols <- c(K80="#0072B2",I21="#D55E00",C18="#009E73",E10_E14="#E69F00",N18="#CC79A7",A40_A41="#6C5B7B")

theme_pub <- function(base_size=7,base_family="Arial") theme_classic(base_size=base_size,base_family=base_family)+
  theme(axis.line=element_line(linewidth=.35,colour=pal["ink"]),axis.ticks=element_line(linewidth=.35),
        axis.title=element_text(size=base_size),axis.text=element_text(size=base_size-.35),
        plot.title=element_text(size=base_size+.55,face="bold"),
        plot.subtitle=element_text(size=base_size-.15,colour=pal["gray"],lineheight=.95),
        plot.caption=element_text(size=base_size-.8,colour=pal["gray"],hjust=0),
        legend.title=element_text(size=base_size-.2,face="bold"),legend.text=element_text(size=base_size-.45),
        legend.key.height=unit(3,"mm"),legend.key.width=unit(4,"mm"),
        strip.text=element_text(size=base_size-.05,face="bold"),
        panel.grid.major.y=element_line(linewidth=.23,colour="#ECECEC"),panel.grid.minor=element_blank(),
        plot.margin=margin(3,4,3,4))
theme_set(theme_pub())
tag_theme <- theme(plot.tag=element_text(family="Arial",size=9,face="bold"),plot.tag.position=c(0,1))

read_rel <- function(rel){
  x <- utils::read.csv(file.path(root,sub("^tables/","",rel)),check.names=FALSE,stringsAsFactors=FALSE,encoding="UTF-8")
  names(x)[1] <- sub("^\\ufeff","",names(x)[1]); as_tibble(x)
}
write_source <- function(x,name){p<-file.path(source_dir,name);utils::write.csv(as.data.frame(x),p,row.names=FALSE,fileEncoding="UTF-8");p}
save_pub <- function(plot,stem,width_mm=183,height_mm=150){
  w<-width_mm/25.4;h<-height_mm/25.4
  svglite::svglite(file.path(out_dir,paste0(stem,".svg")),width=w,height=h,bg="white",system_fonts=list(Arial="Arial"));print(plot);dev.off()
  grDevices::cairo_pdf(file.path(out_dir,paste0(stem,".pdf")),width=w,height=h,family="Arial",bg="white");print(plot);dev.off()
  ragg::agg_png(file.path(out_dir,paste0(stem,".png")),width=w,height=h,units="in",res=300,background="white");print(plot);dev.off()
  ragg::agg_tiff(file.path(out_dir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",background="white");print(plot);dev.off()
}
blank_panel <- function(xlim=c(0,10),ylim=c(0,10)) ggplot()+coord_cartesian(xlim=xlim,ylim=ylim,clip="off")+
  theme_void(base_family="Arial")+theme(plot.title=element_text(size=7.55,face="bold"),plot.margin=margin(3,4,3,4))
add_box <- function(p,xmin,xmax,ymin,ymax,label,fill="white",colour=pal["ink"],size=2.35,fontface="plain") p+
  annotate("rect",xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill,colour=colour,linewidth=.45)+
  annotate("text",x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label,family="Arial",size=size,fontface=fontface,lineheight=.95)

annual <- read_rel("tables/k80_annual_main.csv")
std <- read_rel("tables/k80_standardized_annual.csv")
contrasts <- read_rel("tables/k80_standardized_contrasts.csv")
flow <- read_rel("tables/cohort_annual_flow.csv")
comp <- read_rel("tables/k80_composition_ucf.csv")
decomp <- read_rel("tables/k80_kitagawa_decomposition.csv")
decomp_uncertainty <- read_rel("tables/k80_kitagawa_decomposition_uncertainty.csv")
multi <- read_rel("tables/k80_multilabel_sensitivity.csv")
dest <- read_rel("tables/ucd_destination_annual.csv")
dest_con <- read_rel("tables/ucd_destination_contrasts.csv")
dest_con_simultaneous <- read_rel("tables/ucd_destination_contrasts_simultaneous.csv")
dest_other <- read_rel("tables/ucd_destination_other_ucd.csv")
cross_std <- read_rel("tables/cross_disease_standardized_annual.csv")
cross_con <- read_rel("tables/cross_disease_standardized_contrasts.csv")
race <- read_rel("tables/race_ucf_annual.csv")
diag <- read_rel("tables/interaction_model_diagnostics.csv")
emp <- read_rel("tables/interaction_empirical_marginal_sensitivity.csv")
robust_current <- read_rel("tables/robustness_matrix_current.csv")
era_evidence <- read_rel("tables/era_break_evidence.csv")
composition_sensitivity <- read_rel("tables/k80_composition_sensitivity.csv")

stopifnot(sum(annual$A_record_axis_K80)==51084,sum(annual$B_main_A_and_UCD_K80)==27514)
stopifnot(all(annual$gap_A_minus_B==annual$A_record_axis_K80-annual$B_main_A_and_UCD_K80))
stopifnot(all(abs(dest%>%group_by(year)%>%summarise(s=sum(std_probability),.groups="drop")%>%pull(s)-1)<1e-8))

# Figure 1 --------------------------------------------------------------------
total_records<-sum(flow$total_records);resident_records<-sum(flow$resident_records)
A_total<-sum(flow$A);B_total<-sum(flow$B);gap_total<-sum(flow$gap)
official_total<-sum(annual$B_official_all_UCD_K80);orphan_total<-sum(annual$orphan_UCD_K80)
p1a<-blank_panel(c(1997.5,2025.5),c(0,10))+labs(title="Data scope and coding landmarks")+
  annotate("rect",xmin=1998.7,xmax=2024.3,ymin=6.3,ymax=8.8,fill="#EAF3F8",colour=pal["blue"],linewidth=.5)+
  annotate("text",x=2011.5,y=7.55,label="NCHS Multiple Cause of Death\n26 annual public-use files, 1999-2024",family="Arial",size=2.7,fontface="bold")+
  annotate("segment",x=1999,xend=2024,y=4.5,yend=4.5,linewidth=.65)+
  annotate("point",x=c(1999,2003,2018,2020,2021,2022,2024),y=4.5,size=2,
           colour=c(pal["blue"],pal["gray"],pal["green"],pal["purple"],pal["vermillion"],pal["orange"],pal["blue"]))+
  annotate("text",x=c(1999,2003,2018,2019.35,2021,2022.7,2024),y=c(3.45,3.45,3.45,2.55,3.45,2.55,3.45),
           label=c("1999","2003\ncertificate","2018\nRace Recode 40","2020\nbridged end","2021\nRace Recode 40","2022\nMedCoder","2024"),family="Arial",size=1.7,lineheight=.85)+
  annotate("text",x=2011.5,y=.85,label="Resident deaths; K80 primary boundary; K80-K83 sensitivity only",family="Arial",size=2.05,colour=pal["gray"])
p1b<-blank_panel()+labs(title="Cohort construction")
p1b<-add_box(p1b,.6,9.4,8,9.5,sprintf("All parsed records  n=%s",comma(total_records)),pal["pale"],size=2.35)+
  annotate("segment",x=5,xend=5,y=8,yend=7.2,arrow=arrow(length=unit(1.5,"mm")),linewidth=.45)
p1b<-add_box(p1b,.6,9.4,5.6,7.2,sprintf("U.S. residents  n=%s",comma(resident_records)),"#F0F7F4",pal["green"],2.35)+
  annotate("segment",x=5,xend=5,y=5.6,yend=4.8,arrow=arrow(length=unit(1.5,"mm")),linewidth=.45)
p1b<-add_box(p1b,.6,9.4,3,4.8,sprintf("A: record-axis K80-coded deaths  n=%s",comma(A_total)),"#EAF3F8",pal["blue"],2.45,"bold")+
  annotate("segment",x=5,xend=2.7,y=3,yend=2.2,arrow=arrow(length=unit(1.5,"mm")),linewidth=.45)+
  annotate("segment",x=5,xend=7.3,y=3,yend=2.2,arrow=arrow(length=unit(1.5,"mm")),linewidth=.45)
p1b<-add_box(p1b,.2,4.8,.4,2.2,sprintf("B: A=1 AND UCD=K80\nn=%s",comma(B_total)),"#DDF1EA",pal["green"],2.15,"bold")
p1b<-add_box(p1b,5.2,9.8,.4,2.2,sprintf("Gap: A-B\nn=%s",comma(gap_total)),"#FFF0ED",pal["vermillion"],2.25,"bold")
p1c<-blank_panel()+labs(title="Entity and record axes are distinct")
p1c<-add_box(p1c,.5,9.5,7.8,9.4,"ENTITY AXIS\nPart/line structure and literal certificate entities","#F0F7F4",pal["green"],2.25,"bold")
p1c<-add_box(p1c,.5,9.5,4.7,6.5,"RECORD AXIS (PRIMARY A)\nComputer-processed record-axis K80 codes","#EAF3F8",pal["blue"],2.25,"bold")
p1c<-add_box(p1c,.5,9.5,1.3,3.6,sprintf("OFFICIAL-UCD K80 RECONCILIATION CONTROL\nPrimary B: A=1 AND UCD=K80  n=%s\nAll official UCD K80=%s; outside A=%s",comma(B_total),comma(official_total),comma(orphan_total)),"#FFF4E6",pal["orange"],1.95,"bold")+
  annotate("text",x=5,y=.55,label="No public individual row identifiers are displayed",family="Arial",size=2,colour=pal["gray"])
p1d<-blank_panel()+labs(title="Two denominators must not be conflated")+
  annotate("text",x=5,y=8.6,label="Conditional UCF = B / A",family="Arial",size=3.6,fontface="bold",colour=pal["blue"])+
  annotate("text",x=5,y=7.2,label="Conditional UCD complement = (A-B) / A",family="Arial",size=2.9,fontface="bold",colour=pal["vermillion"])
p1d<-add_box(p1d,.7,9.3,4.3,6.1,"POPULATION RATE\nK80 deaths / U.S. population",pal["pale"],pal["gray"],2.35)
p1d<-add_box(p1d,.7,9.3,1.7,3.5,"CONDITIONAL RECORDED UCD FRACTION\nUCD=K80 among record-axis K80-coded deaths","#EAF3F8",pal["blue"],2.2)+
  annotate("text",x=5,y=.55,label="Recorded attribution; not clinical diagnosis or causal proof",family="Arial",size=2,colour=pal["gray"])
write_source(bind_rows(
  flow%>%transmute(panel="A/B",year,metric="annual_file_flow",total_records,resident_records,A,B,gap),
  annual%>%transmute(panel="C",year,metric="axis_reconciliation",total_records=NA_real_,resident_records=NA_real_,A=A_record_axis_K80,B=B_main_A_and_UCD_K80,gap=gap_A_minus_B,B_official=B_official_all_UCD_K80,orphan=orphan_UCD_K80)
),"Figure1_source_data.csv")
fig1<-((p1a|p1b)/(p1c|p1d))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig1,"Figure1_study_universe_estimands",183,150)

# Figure 2 --------------------------------------------------------------------
std_primary<-std%>%filter(scheme=="primary_2018_2024")
counts_long<-annual%>%select(year,A=A_record_axis_K80,B=B_main_A_and_UCD_K80,gap=gap_A_minus_B)%>%pivot_longer(c(A,B,gap),names_to="series",values_to="count")
p2a<-ggplot(counts_long,aes(year,count,colour=series))+geom_line(linewidth=.72)+geom_point(size=1)+
  scale_colour_manual(values=c(A=unname(pal["blue"]),B=unname(pal["green"]),gap=unname(pal["vermillion"])),labels=c(A="Record-axis K80-coded (A)",B="A=1 & UCD=K80 (B)",gap="Complement (A-B)"))+
  scale_x_continuous(breaks=seq(2000,2024,4))+scale_y_continuous(labels=comma)+labs(title="Annual record counts",x=NULL,y="Deaths",colour=NULL)+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+theme(legend.position="top",legend.text=element_text(size=5.4),legend.key.width=unit(3,"mm"))
ucf_long<-annual%>%transmute(year,estimate=UCF,lo=UCF_CI_lo,hi=UCF_CI_hi,method="Crude")%>%
  bind_rows(std_primary%>%transmute(year,estimate=std_UCF,lo=CI_lo,hi=CI_hi,method="Age-sex standardized"))
p2b<-ggplot(ucf_long,aes(year,estimate,colour=method,fill=method))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.1,colour=NA)+geom_line(linewidth=.78)+
  scale_colour_manual(values=c(Crude=unname(pal["gray"]),`Age-sex standardized`=unname(pal["blue"])))+scale_fill_manual(values=c(Crude=unname(pal["gray"]),`Age-sex standardized`=unname(pal["blue"])))+
  scale_x_continuous(breaks=seq(2000,2024,4))+scale_y_continuous(limits=c(0,.70),labels=percent_format(accuracy=1))+
  labs(title="Crude and standardized UCF",subtitle="95% confidence intervals shown",x=NULL,y="Underlying-cause fraction",colour=NULL,fill=NULL)+theme(legend.position="top")
landmarks<-contrasts%>%filter(scheme=="primary_2018_2024",contrast_or_year%in%c("1999","2015","2024"))%>%mutate(contrast_or_year=factor(contrast_or_year,levels=c("1999","2015","2024")))
p2c<-ggplot(landmarks,aes(contrast_or_year,std_UCF))+geom_errorbar(aes(ymin=CI_lo,ymax=CI_hi),width=.12,linewidth=.55,colour=pal["blue"])+
  geom_point(size=2.8,shape=21,fill="white",colour=pal["blue"],stroke=.8)+geom_text(aes(label=percent(std_UCF,accuracy=.1)),vjust=-1.35,size=2.2,family="Arial")+
  scale_y_continuous(limits=c(0,.70),labels=percent_format(accuracy=1))+labs(title="Endpoint landmarks",subtitle="2015 is the descriptive maximum; not a tested peak",x=NULL,y="Standardized UCF (95% CI)")
rd_scheme<-contrasts%>%filter(contrast_or_year=="2024_minus_1999")%>%mutate(scheme_label=recode(scheme,primary_2018_2024="2018-2024 weights (primary)",sensitivity_1999="1999 weights",sensitivity_full_period="1999-2024 weights"),scheme_label=factor(scheme_label,levels=rev(c("2018-2024 weights (primary)","1999 weights","1999-2024 weights"))))
p2d<-ggplot(rd_scheme,aes(RD,scheme_label,colour=scheme))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+geom_errorbarh(aes(xmin=RD_CI_lo,xmax=RD_CI_hi),height=0,linewidth=.65)+geom_point(size=2.4)+
  scale_colour_manual(values=c(primary_2018_2024=unname(pal["blue"]),sensitivity_1999=unname(pal["orange"]),sensitivity_full_period=unname(pal["green"])),guide="none")+
  scale_x_continuous(limits=c(0,.16),labels=percent_format(accuracy=1))+labs(title="Fixed-standard sensitivity",subtitle="2024 minus 1999 risk difference",x="Probability difference (95% CI)",y=NULL)
write_source(bind_rows(
  counts_long%>%transmute(panel="A",year,series,estimate=count,lo=NA_real_,hi=NA_real_,scheme=NA_character_),
  ucf_long%>%transmute(panel="B",year,series=method,estimate,lo,hi,scheme=NA_character_),
  landmarks%>%transmute(panel="C",year=as.integer(as.character(contrast_or_year)),series="landmark standardized UCF",estimate=std_UCF,lo=CI_lo,hi=CI_hi,scheme),
  rd_scheme%>%transmute(panel="D",year=NA_integer_,series="2024-minus-1999 RD",estimate=RD,lo=RD_CI_lo,hi=RD_CI_hi,scheme)
),"Figure2_source_data.csv")
fig2<-((p2a|p2b)/(p2c|p2d))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig2,"Figure2_annual_K80_UCF",183,150)

# Figure 3 --------------------------------------------------------------------
subtype<-comp%>%filter(group=="subtype",year%in%c(1999,2024))%>%mutate(stratum=factor(stratum,levels=c("K80.0","K80.1","K80.2","K80.3","K80.4","K80.5","K80.8")),year=factor(year))
p3a<-ggplot(subtype,aes(year,composition_weight,fill=stratum))+geom_col(width=.72,colour="white",linewidth=.15)+
  scale_fill_manual(values=c("#0072B2","#56B4E9","#B8D9EA","#009E73","#E69F00","#D55E00","#6F6F6F"),drop=FALSE)+
  scale_y_continuous(labels=percent_format(accuracy=1),expand=c(0,0))+labs(title="K80 subtype composition shifted",x=NULL,y="Share of record-axis K80-coded deaths",fill="Priority subtype")+theme(legend.position="right",legend.key.height=unit(2.5,"mm"))
p3b<-ggplot(subtype,aes(UCF,stratum,colour=year,group=stratum))+geom_line(linewidth=.55,colour=pal["light"])+
  geom_errorbarh(aes(xmin=UCF_CI_lo,xmax=UCF_CI_hi),height=0,linewidth=.5,position=position_dodge(width=.28))+geom_point(size=2.1,position=position_dodge(width=.28))+
  scale_colour_manual(values=c(`1999`=unname(pal["gray"]),`2024`=unname(pal["blue"])))+scale_x_continuous(limits=c(0,1.02),labels=percent_format(accuracy=1))+
  labs(title="Subtype-specific UCF",subtitle="Priority ontology; 95% CIs",x="Underlying-cause fraction",y=NULL,colour=NULL)+theme(legend.position="top")
decomp_long<-bind_rows(
  decomp_uncertainty%>%transmute(group,bootstrap_reps,component="Composition",estimate=composition_point,lo=composition_percentile_ci_lo,hi=composition_percentile_ci_hi,method),
  decomp_uncertainty%>%transmute(group,bootstrap_reps,component="Within-stratum UCF change",estimate=within_stratum_point,lo=within_stratum_percentile_ci_lo,hi=within_stratum_percentile_ci_hi,method)
)%>%mutate(
  component=factor(component,levels=c("Composition","Within-stratum UCF change")),
  group_label=recode(group,subtype="K80 subtype",severity="Severity grouping",entity_part="Part I vs Part II",complexity="Certificate complexity"),
  group_label=factor(group_label,levels=rev(c("K80 subtype","Severity grouping","Part I vs Part II","Certificate complexity")))
)
p3c<-ggplot(decomp_long,aes(estimate,group_label,colour=component))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,linewidth=.62,position=position_dodge(width=.48))+
  geom_point(size=2.2,position=position_dodge(width=.48))+
  scale_colour_manual(values=c(Composition=unname(pal["orange"]),`Within-stratum UCF change`=unname(pal["blue"])))+
  scale_x_continuous(labels=percent_format(accuracy=1),limits=c(-.06,.19))+
  labs(title="Kitagawa decomposition with uncertainty",subtitle="Joint bootstrap percentile 95% CIs\n100,000 replicates per grouping",x="Component of crude 2024-minus-1999 UCF change",y=NULL,colour=NULL)+theme(legend.position="top")
definition_sens<-composition_sensitivity%>%
  filter(analysis%in%c("part_rule_kitagawa_total","complexity_kitagawa_total"),comparison=="1999_vs_2024")%>%
  transmute(
    family=if_else(analysis=="part_rule_kitagawa_total","Part rule","Complexity"),
    scenario=recode(scenario,Part_I_priority="Part I priority",Part_II_priority="Part II priority",Both_separate="Both separate",Both_excluded="Both excluded",record_n="record_n positions",unique_record_axis_icd="Unique ICD codes"),
    Composition=as.numeric(composition_total),
    `Within-stratum UCF`=as.numeric(selection_total),
    total_change=as.numeric(total_change)
  )
definition_sens_long<-definition_sens%>%pivot_longer(c(Composition,`Within-stratum UCF`),names_to="component",values_to="estimate")%>%
  mutate(scenario=factor(scenario,levels=rev(c("Part I priority","Part II priority","Both separate","Both excluded","record_n positions","Unique ICD codes"))))
p3d<-ggplot(definition_sens_long,aes(estimate,scenario,colour=component))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+
  geom_point(size=2.05,position=position_dodge(width=.38))+facet_grid(family~.,scales="free_y",space="free_y")+
  scale_colour_manual(values=c(Composition=unname(pal["orange"]),`Within-stratum UCF`=unname(pal["blue"])))+
  scale_x_continuous(limits=c(-.06,.17),labels=percent_format(accuracy=.1))+
  labs(title="Operational-definition sensitivity",subtitle="337 dual-Part records; 46 records differed by complexity counter",x="Component of crude 2024-minus-1999 UCF change",y=NULL,colour=NULL)+
  theme(legend.position="top",strip.background=element_rect(fill=pal["pale"],colour=NA),legend.text=element_text(size=5.1))
write_source(bind_rows(
  subtype%>%transmute(panel="A/B",group="subtype_priority",stratum=as.character(stratum),year=as.integer(as.character(year)),estimate=UCF,lo=UCF_CI_lo,hi=UCF_CI_hi,composition=composition_weight,component=NA_character_,bootstrap_reps=NA_integer_,method=NA_character_),
  decomp_long%>%transmute(panel="C",group=as.character(group_label),stratum=NA_character_,year=NA_integer_,estimate,lo,hi,composition=NA_real_,component=as.character(component),bootstrap_reps,method),
  definition_sens_long%>%transmute(panel="D",group=family,stratum=as.character(scenario),year=NA_integer_,estimate,lo=NA_real_,hi=NA_real_,composition=total_change,component,bootstrap_reps=NA_integer_,method="full-scan operational-definition sensitivity")
),"Figure3_source_data.csv")
fig3<-((p3a|p3b)/(p3c|p3d))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig3,"Figure3_composition_decomposition",183,155)

# Figure 4 --------------------------------------------------------------------
dest<-dest%>%mutate(destination=factor(destination,levels=dest_order))
p4a<-ggplot(dest,aes(year,std_probability,fill=destination))+geom_area(position="stack",colour="white",linewidth=.08)+
  scale_fill_manual(values=dest_cols,drop=FALSE,labels=function(x)gsub("_"," ",x))+scale_x_continuous(breaks=seq(2000,2024,4))+
  scale_y_continuous(labels=percent_format(accuracy=1),expand=c(0,0))+
  labs(title="Annual standardized destination composition",subtitle="Nine mutually exclusive categories; annual sum=100%",x=NULL,y="Share among gap records",fill=NULL)+
  theme(legend.position="right",legend.key.height=unit(2.4,"mm"),legend.text=element_text(size=5.2))
leading_names<-dest%>%group_by(destination)%>%summarise(m=mean(std_probability),.groups="drop")%>%arrange(desc(m))%>%slice_head(n=4)%>%pull(destination)%>%as.character()
leading<-dest%>%filter(as.character(destination)%in%leading_names)
p4b<-ggplot(leading,aes(year,std_probability,colour=destination,fill=destination))+geom_ribbon(aes(ymin=std_CI_lo,ymax=std_CI_hi),alpha=.08,colour=NA)+geom_line(linewidth=.75)+
  scale_colour_manual(values=dest_cols,labels=function(x)gsub("_"," ",x))+scale_fill_manual(values=dest_cols,labels=function(x)gsub("_"," ",x))+
  scale_x_continuous(breaks=seq(2000,2024,4))+scale_y_continuous(limits=c(0,.50),labels=percent_format(accuracy=1))+
  labs(title="Leading destination trajectories",subtitle="95% confidence intervals",x=NULL,y="Standardized probability",colour=NULL,fill=NULL)+guides(colour=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE))+theme(legend.position="top",legend.text=element_text(size=5.1),legend.key.width=unit(3.2,"mm"))
dest_con<-dest_con%>%mutate(destination=factor(destination,levels=rev(dest_order)))
dest_sim<-dest_con_simultaneous%>%mutate(destination=factor(destination,levels=rev(dest_order)),simultaneous_excludes_zero=tolower(as.character(simultaneous_excludes_zero))=="true")
p4c<-ggplot(dest_sim,aes(RD_2024_minus_1999_point,destination,colour=destination,shape=simultaneous_excludes_zero))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+
  geom_errorbarh(aes(xmin=simultaneous_95_ci_lo,xmax=simultaneous_95_ci_hi),height=0,linewidth=.65)+geom_point(size=2.1)+
  scale_colour_manual(values=dest_cols,guide="none")+scale_shape_manual(values=c(`FALSE`=21,`TRUE`=16),guide="none")+
  scale_y_discrete(labels=function(x)gsub("_"," ",x))+scale_x_continuous(labels=percent_format(accuracy=1))+
  labs(title="Destination-specific change",subtitle="2024-1999 RD; max-t simultaneous 95% CIs\nFWER controlled across 9 categories",x="Standardized probability difference",y=NULL)
other_2024<-dest_other%>%filter(year==2024,destination=="OTHER")%>%arrange(desc(raw_n))%>%slice_head(n=10)%>%mutate(ucd=factor(ucd,levels=rev(ucd)))
p4d<-ggplot(other_2024,aes(raw_n,ucd))+geom_col(fill=pal["gray"],width=.68)+geom_text(aes(label=raw_n),hjust=-.15,size=1.9,family="Arial")+
  scale_x_continuous(expand=expansion(mult=c(0,.18)))+labs(title="Residual OTHER drill-down",subtitle="Ten most frequent UCD codes in 2024",x="Gap records",y="Underlying-cause ICD-10 code")
write_source(bind_rows(
  dest%>%transmute(panel="A/B",year,destination=as.character(destination),estimate=std_probability,lo=std_CI_lo,hi=std_CI_hi,count=raw_n,code=NA_character_,interval_type="pointwise annual 95% CI",bootstrap_reps=NA_integer_,simultaneous_excludes_zero=NA),
  dest_sim%>%transmute(panel="C",year=NA_integer_,destination=as.character(destination),estimate=RD_2024_minus_1999_point,lo=simultaneous_95_ci_lo,hi=simultaneous_95_ci_hi,count=NA_real_,code=NA_character_,interval_type="FWER max-t simultaneous 95% CI across nine categories",bootstrap_reps,simultaneous_excludes_zero),
  other_2024%>%transmute(panel="D",year,destination,estimate=raw_probability,lo=NA_real_,hi=NA_real_,count=raw_n,code=as.character(ucd),interval_type=NA_character_,bootstrap_reps=NA_integer_,simultaneous_excludes_zero=NA)
),"Figure4_source_data.csv")
fig4<-((p4a|p4b)/(p4c|p4d))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig4,"Figure4_UCD_destinations",183,155)

# Figure 5 --------------------------------------------------------------------
k80_cross<-std_primary%>%transmute(disease="K80",scheme,year,std_UCF,CI_lo,CI_hi,variance,effective_weight_sum,supported_cells)
cross_all<-bind_rows(k80_cross,cross_std)%>%mutate(disease=factor(disease,levels=disease_order),disease_label=factor(disease_labels[as.character(disease)],levels=unname(disease_labels[disease_order])))
k80_con<-contrasts%>%filter(scheme=="primary_2018_2024")%>%transmute(disease="K80",contrast_or_year,std_UCF,CI_lo,CI_hi,RD,RD_CI_lo,RD_CI_hi,RR,RR_CI_lo,RR_CI_hi)
cross_con_all<-bind_rows(k80_con,cross_con)%>%filter(contrast_or_year=="2024_minus_1999")%>%mutate(disease=factor(disease,levels=rev(disease_order)),disease_label=factor(disease_labels[as.character(disease)],levels=unname(disease_labels[rev(disease_order)])))
p5a<-ggplot(cross_all,aes(year,std_UCF,colour=disease,fill=disease))+geom_ribbon(aes(ymin=CI_lo,ymax=CI_hi),alpha=.08,colour=NA)+geom_line(linewidth=.68)+
  facet_wrap(~disease_label,scales="free_y",ncol=2)+scale_colour_manual(values=disease_cols,guide="none")+scale_fill_manual(values=disease_cols,guide="none")+
  scale_x_continuous(breaks=c(2000,2008,2016,2024))+scale_y_continuous(labels=percent_format(accuracy=1))+
  labs(title="Within-disease standardized trajectories",subtitle="Each disease uses its own fixed 2018-2024 age-sex distribution of record-axis coded deaths",x=NULL,y="Standardized UCF")+
  theme(strip.background=element_rect(fill=pal["pale"],colour=NA),panel.spacing=unit(3,"mm"))
p5b<-ggplot(cross_con_all,aes(RD,disease_label,colour=disease))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+
  geom_errorbarh(aes(xmin=RD_CI_lo,xmax=RD_CI_hi),height=0,linewidth=.62)+geom_point(size=2.2)+scale_colour_manual(values=disease_cols,guide="none")+
  scale_x_continuous(labels=percent_format(accuracy=1))+labs(title="Disease-specific endpoint contrasts",subtitle="2024-1999 RD (95% CI)\nDo not rank across disease-specific targets",x="Within-disease probability difference",y=NULL)
p5c<-ggplot(cross_con_all,aes(RR,disease_label,colour=disease))+geom_vline(xintercept=1,colour=pal["gray"],linewidth=.35)+
  geom_errorbarh(aes(xmin=RR_CI_lo,xmax=RR_CI_hi),height=0,linewidth=.62)+geom_point(size=2.2)+scale_colour_manual(values=disease_cols,guide="none")+
  scale_x_continuous(trans="log10",breaks=c(.8,.9,1,1.1,1.3),labels=number_format(accuracy=.01))+
  labs(title="Disease-specific relative contrasts",subtitle="2024/1999 RR (95% CI)\nDo not rank across disease-specific targets",x="Within-disease risk ratio (log scale)",y=NULL)
p5d<-blank_panel()+labs(title="Interpretation boundary")
p5d<-add_box(p5d,.5,9.5,7.5,9.3,"SUPPORTED\nDirection of 1999-to-2024 change within each disease","#EAF6F1",pal["green"],2.2,"bold")
p5d<-add_box(p5d,.5,9.5,4.3,6.5,"NOT SUPPORTED\nCross-disease ranking of absolute UCF, RD, or RR magnitudes under disease-specific standards","#FFF0ED",pal["vermillion"],2.05,"bold")
p5d<-add_box(p5d,.5,9.5,1.2,3.2,"RESULT\nDirections differ among these six NCHS examples; no universal direction is inferred","#EAF3F8",pal["blue"],2.1,"bold")
write_source(bind_rows(
  cross_all%>%transmute(panel="A",disease=as.character(disease),disease_label,year,estimate=std_UCF,lo=CI_lo,hi=CI_hi,metric="standardized UCF",scheme="own disease fixed 2018-2024 age-sex distribution of record-axis coded deaths"),
  cross_con_all%>%transmute(panel="B",disease=as.character(disease),disease_label,year=NA_integer_,estimate=RD,lo=RD_CI_lo,hi=RD_CI_hi,metric="2024-minus-1999 within-disease RD; magnitude not cross-disease comparable",scheme="own disease 2018-2024 record-axis-coded age-sex weights"),
  cross_con_all%>%transmute(panel="C",disease=as.character(disease),disease_label,year=NA_integer_,estimate=RR,lo=RR_CI_lo,hi=RR_CI_hi,metric="2024-over-1999 within-disease RR; magnitude not cross-disease comparable",scheme="own disease 2018-2024 record-axis-coded age-sex weights")
),"Figure5_source_data.csv")
fig5<-(p5a/(p5b|p5c)/p5d)+plot_layout(heights=c(1.65,1,.72))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig5,"Figure5_cross_disease_benchmark",183,180)

# Figure 6 --------------------------------------------------------------------
axis_sens<-annual%>%transmute(year,`Primary B/A`=UCF,`A-prime sensitivity`=B_official_all_UCD_K80/A_prime_A_union_B_official)%>%pivot_longer(-year,names_to="analysis",values_to="estimate")
p6a<-ggplot(axis_sens,aes(year,estimate,colour=analysis,linetype=analysis))+geom_line(linewidth=.78)+geom_point(data=axis_sens%>%filter(year%in%c(1999,2024)),size=1.5)+scale_colour_manual(values=c(`Primary B/A`=unname(pal["blue"]),`A-prime sensitivity`=unname(pal["orange"])))+scale_linetype_manual(values=c(`Primary B/A`="solid",`A-prime sensitivity`="22"))+
  scale_x_continuous(breaks=seq(2000,2024,4))+scale_y_continuous(limits=c(0,.70),labels=percent_format(accuracy=1))+
  labs(title="Record-axis and A-prime sensitivity",subtitle="Official UCD-only orphans enter A-prime",x=NULL,y="Crude UCF",colour=NULL,linetype=NULL)+theme(legend.position="top")
rob_summary<-robust_current%>%
  mutate(
    analysis=recode(analysis_id,
      `ROBUST-01`="Primary K80",
      `ROBUST-02`="A-prime crude",
      `ROBUST-03`="Entity-axis K80",
      `ROBUST-04`="Expanded K80-K83",
      `ROBUST-05`="Age >=25 years",
      `ROBUST-06`="Age >=65 years",
      `ROBUST-07`="1999 standard",
      `ROBUST-08`="Full-period standard"),
    family=case_when(
      analysis_id=="ROBUST-01"~"Primary",
      analysis_id%in%c("ROBUST-02","ROBUST-03")~"Axis/denominator",
      analysis_id=="ROBUST-04"~"Disease boundary",
      analysis_id%in%c("ROBUST-05","ROBUST-06")~"Age restriction",
      TRUE~"Standard weights"),
    analysis=factor(analysis,levels=rev(analysis))
  )
p6b<-ggplot(rob_summary,aes(rd_2024_minus_1999,analysis,colour=family))+geom_vline(xintercept=0,colour=pal["gray"],linewidth=.35)+
  geom_errorbarh(aes(xmin=rd_ci_lower,xmax=rd_ci_upper),height=0,linewidth=.58)+geom_point(size=2.05)+
  scale_colour_manual(values=c(Primary=unname(pal["blue"]),`Axis/denominator`=unname(pal["purple"]),`Disease boundary`=unname(pal["orange"]),`Age restriction`=unname(pal["green"]),`Standard weights`=unname(pal["gray"])))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  scale_x_continuous(limits=c(0,.18),labels=percent_format(accuracy=1))+labs(title="Current-candidate sensitivity matrix",subtitle="2024-minus-1999 RD with 95% CIs; all values are locally bound",x="Probability difference",y=NULL,colour=NULL)+theme(legend.position="top",legend.text=element_text(size=5),legend.key.width=unit(2.8,"mm"))
era_df<-bind_rows(
  tibble(period="Pre_2020",years="1999-2019",n_annual_points=21L,mean_UCF=mean(std_primary$std_UCF[std_primary$year<=2019])),
  era_evidence%>%select(period,years,n_annual_points,mean_UCF)
)%>%mutate(era=recode(period,Pre_2020="Pre-2020",COVID_2020_2021="COVID 2020-2021",MedCoder_2022_2024="MedCoder 2022-2024"),era=factor(era,levels=c("Pre-2020","COVID 2020-2021","MedCoder 2022-2024")),estimate=mean_UCF)
p6c<-ggplot(era_df,aes(era,estimate,fill=era))+geom_col(width=.68)+geom_text(aes(label=percent(estimate,accuracy=.1)),vjust=-.4,size=2.1,family="Arial")+
  scale_fill_manual(values=c(`Pre-2020`=unname(pal["gray"]),`COVID 2020-2021`=unname(pal["purple"]),`MedCoder 2022-2024`=unname(pal["orange"])),guide="none")+
  scale_y_continuous(limits=c(0,.64),labels=percent_format(accuracy=1))+labs(title="Institutional windows",subtitle="Descriptive fixed-standard era estimates",x=NULL,y="Standardized UCF")+theme(axis.text.x=element_text(angle=18,hjust=1))
check_df<-tibble(check=factor(c("Primary weight support","Destination ontology","Decomposition identity","Definition sensitivity","Race Recode 40","Cross-disease weights"),levels=rev(c("Primary weight support","Destination ontology","Decomposition identity","Definition sensitivity","Race Recode 40","Cross-disease weights"))),status=c("Complete","1-135 local; yearly 100%","Residual ~0",ifelse(nrow(composition_sensitivity)>0,"Executed","Missing"),"Six groups incl. 2021","Disease-specific"),class=c("pass","pass","pass",ifelse(nrow(composition_sensitivity)>0,"pass","report"),"pass","boundary"))
p6d<-ggplot(check_df,aes(1,check,fill=class))+geom_tile(width=.88,height=.72,colour="white")+geom_text(aes(label=status),family="Arial",size=2.15)+
  scale_fill_manual(values=c(pass=alpha(unname(pal["green"]),.4),report=alpha(unname(pal["blue"]),.35),boundary=alpha(unname(pal["orange"]),.35)),guide="none")+
  scale_x_continuous(NULL,breaks=NULL)+labs(title="Measurement-boundary matrix",subtitle="Race Recode 40 is mapped at positions 489–490; 2021 is included",y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())
write_source(bind_rows(
  axis_sens%>%transmute(panel="A",analysis,year,estimate,lo=NA_real_,hi=NA_real_,note="axis/A-prime sensitivity"),
  rob_summary%>%transmute(panel="B",analysis=as.character(analysis),year=NA_integer_,estimate=rd_2024_minus_1999,lo=rd_ci_lower,hi=rd_ci_upper,note=family),
  era_df%>%transmute(panel="C",analysis=as.character(era),year=NA_integer_,estimate,lo=NA_real_,hi=NA_real_,note="descriptive era standardized UCF"),
  check_df%>%transmute(panel="D",analysis=as.character(check),year=NA_integer_,estimate=NA_real_,lo=NA_real_,hi=NA_real_,note=status)
),"Figure6_source_data.csv")
fig6<-((p6a|p6b)/(p6c|p6d))+plot_annotation(tag_levels="a")&tag_theme
save_pub(fig6,"Figure6_robustness_boundaries",183,150)

# Supplementary Figure S1 -----------------------------------------------------
pS1a<-ggplot(dest,aes(year,std_probability,colour=destination))+geom_line(linewidth=.62)+
  facet_wrap(~destination,scales="free_y",ncol=3,labeller=labeller(destination=function(x)gsub("_"," ",x)))+
  scale_colour_manual(values=dest_cols,guide="none")+scale_x_continuous(breaks=c(2000,2012,2024))+scale_y_continuous(labels=percent_format(accuracy=1))+
  labs(title="All nine destination trajectories",subtitle="Fixed 2018-2024 gap age-sex standard; panels use free y-scales",x=NULL,y="Standardized probability")
other_top<-dest_other%>%filter(destination=="OTHER")%>%group_by(ucd)%>%summarise(raw_n=sum(raw_n),.groups="drop")%>%arrange(desc(raw_n))%>%slice_head(n=20)%>%mutate(ucd=factor(ucd,levels=rev(ucd)))
pS1b<-ggplot(other_top,aes(raw_n,ucd))+geom_col(fill=pal["gray"],width=.68)+scale_x_continuous(labels=comma,expand=expansion(mult=c(0,.08)))+
  labs(title="Residual OTHER codes across 1999-2024",subtitle="Twenty most frequent underlying-cause ICD-10 codes",x="Gap records",y=NULL)
figS1<-pS1a/pS1b+plot_layout(heights=c(1.65,1))+plot_annotation(tag_levels="a")&tag_theme
write_source(bind_rows(
  dest%>%transmute(panel="A",year,group=as.character(destination),code=NA_character_,estimate=std_probability,count=raw_n),
  other_top%>%transmute(panel="B",year=NA_integer_,group="OTHER",code=as.character(ucd),estimate=NA_real_,count=raw_n)
),"Supplementary_Figure_S1_source_data.csv")
save_pub(figS1,"Supplementary_Figure_S1_destinations",183,180)

# Supplementary Figure S2 -----------------------------------------------------
race_plot<-race%>%
  arrange(regime,race,year)%>%
  mutate(regime=factor(regime,levels=c("bridged","single_race")),estimate=ifelse(!is.na(std_UCF),std_UCF,UCF),
    standardization_status=if_else(!is.na(std_UCF),
      "common pooled age-sex distribution among record-axis coded deaths within regime",standardization_status),
    estimate_type=ifelse(!is.na(std_UCF),"Standardized (common regime weights)","Crude only / not estimable"))%>%
  group_by(regime,race)%>%
  mutate(display_segment=cumsum(row_number()==1L | year!=lag(year)+1L | estimate_type!=lag(estimate_type) | is.na(estimate) | is.na(lag(estimate))))%>%
  ungroup()
pS2a<-ggplot(race_plot%>%filter(regime=="bridged"),aes(year,estimate,colour=race,linetype=estimate_type,group=interaction(race,estimate_type)))+geom_line(linewidth=.68,na.rm=TRUE)+
  scale_x_continuous(breaks=seq(2000,2020,4))+scale_y_continuous(labels=percent_format(accuracy=1))+
  labs(title="Bridged-race regime, 1999-2020",x=NULL,y="UCF",colour="Race",linetype=NULL)+theme(legend.position="right")
pS2b<-ggplot(race_plot%>%filter(regime=="single_race"),aes(year,estimate,colour=race,linetype=estimate_type,group=interaction(race,estimate_type,display_segment)))+
  geom_line(linewidth=.68,na.rm=TRUE)+geom_point(size=1.15,na.rm=TRUE)+
  scale_x_continuous(breaks=2018:2024)+scale_y_continuous(labels=percent_format(accuracy=1))+
  labs(title="Race Recode 40 regime, 2018-2024",subtitle="Official Race Recode 40 (489-490); 2021 included; no imputation",x=NULL,y="UCF",colour="Race",linetype=NULL)+theme(legend.position="right")
figS2<-pS2a/pS2b+plot_annotation(tag_levels="a")&tag_theme
write_source(race_plot,"Supplementary_Figure_S2_source_data.csv")
save_pub(figS2,"Supplementary_Figure_S2_race_regimes",183,160)

# Supplementary Figure S3 -----------------------------------------------------
diag_plot<-diag%>%filter(!is.na(deviance_per_df))%>%mutate(model_label=recode(model,main_effects="Main effects",year_x_age="Year x age",year_x_sex="Year x sex",age_x_sex="Age x sex",all_pairwise="All pairwise",saturated_year_x_age_x_sex="Saturated"),model_label=factor(model_label,levels=rev(c("Main effects","Year x age","Year x sex","Age x sex","All pairwise","Saturated"))))
pS3a<-ggplot(diag_plot,aes(deviance_per_df,model_label))+geom_vline(xintercept=1,linetype="dashed",colour=pal["gray"],linewidth=.4)+geom_point(colour=pal["blue"],size=2.2)+
  labs(title="Aggregated-binomial model diagnostics",subtitle="Deviance per residual degree of freedom",x="Deviance / df",y=NULL)
support_plot<-bind_rows(
  std%>%transmute(module=paste0("K80: ",scheme),year,support=effective_weight_sum),
  cross_std%>%transmute(module=paste0("Comparator: ",disease),year,support=effective_weight_sum),
  dest%>%group_by(year)%>%summarise(support=min(effective_weight_sum),.groups="drop")%>%mutate(module="Destinations")
)
pS3b<-ggplot(support_plot,aes(year,support,colour=module))+geom_line(linewidth=.62)+geom_hline(yintercept=1,linetype="dashed",linewidth=.35,colour=pal["gray"])+
  scale_x_continuous(breaks=seq(2000,2024,4))+scale_y_continuous(limits=c(.98,1.001),labels=percent_format(accuracy=.1))+
  labs(title="Direct-standardization support",subtitle="Primary, destination, and comparator modules",x=NULL,y="Effective weight sum",colour=NULL)+theme(legend.position="right",legend.text=element_text(size=5))
emp_delta<-emp%>%select(scheme,year,empirical=std_UCF)%>%left_join(std%>%select(scheme,year,direct=std_UCF),by=c("scheme","year"))%>%mutate(delta=empirical-direct)
pS3c<-ggplot(emp_delta,aes(year,delta,colour=scheme))+geom_hline(yintercept=0,colour=pal["gray"],linewidth=.35)+geom_line(linewidth=.68)+scale_x_continuous(breaks=seq(2000,2024,4))+
  labs(title="Empirical marginal equivalence",subtitle="Saturated year-age-sex marginal minus direct standardization",x=NULL,y="Probability difference",colour="Weight scheme")+theme(legend.position="top")
figS3<-(pS3a|pS3b)/pS3c+plot_annotation(tag_levels="a")&tag_theme
write_source(bind_rows(
  diag_plot%>%transmute(panel="A",module=as.character(model_label),year=NA_integer_,estimate=deviance_per_df,note=status),
  support_plot%>%transmute(panel="B",module,year,estimate=support,note="effective weight sum"),
  emp_delta%>%transmute(panel="C",module=scheme,year,estimate=delta,note="empirical minus direct")
),"Supplementary_Figure_S3_source_data.csv")
save_pub(figS3,"Supplementary_Figure_S3_diagnostics",183,150)

# Automatic numeric/export QA --------------------------------------------------
fig_files<-list.files(out_dir,full.names=TRUE);source_files<-list.files(source_dir,full.names=TRUE)
expected_stems<-c("Figure1_study_universe_estimands","Figure2_annual_K80_UCF","Figure3_composition_decomposition","Figure4_UCD_destinations","Figure5_cross_disease_benchmark","Figure6_robustness_boundaries","Supplementary_Figure_S1_destinations","Supplementary_Figure_S2_race_regimes","Supplementary_Figure_S3_diagnostics")
qa_files<-tibble(path=list.files(out_dir,full.names=TRUE))%>%mutate(file=basename(path),stem=tools::file_path_sans_ext(file),format=tolower(tools::file_ext(file)),bytes=file.info(path)$size,exists_nonempty=file.exists(path)&bytes>1000)
qa_matrix<-expand_grid(stem=expected_stems,format=c("png","svg","pdf","tiff"))%>%left_join(qa_files%>%select(stem,format,bytes,exists_nonempty),by=c("stem","format"))%>%mutate(pass=!is.na(bytes)&exists_nonempty)
png_dims<-lapply(file.path(out_dir,paste0(expected_stems,".png")),function(p){i<-magick::image_info(magick::image_read(p));tibble(stem=tools::file_path_sans_ext(basename(p)),png_width_px=i$width[1],png_height_px=i$height[1])})%>%bind_rows()
svg_text<-tibble(stem=expected_stems,svg_has_editable_text=vapply(file.path(out_dir,paste0(expected_stems,".svg")),function(p)any(grepl("<text",readLines(p,warn=FALSE,encoding="UTF-8"),fixed=TRUE)),logical(1)))
qa_matrix<-qa_matrix%>%left_join(png_dims,by="stem")%>%left_join(svg_text,by="stem")
dest_sum_error<-max(abs(dest%>%group_by(year)%>%summarise(s=sum(std_probability),.groups="drop")%>%pull(s)-1))
decomp_point_check<-decomp_uncertainty%>%select(group,composition_point,within_stratum_point)%>%
  inner_join(decomp%>%select(group,composition_component,selection_component),by="group")
decomp_point_max_diff<-max(abs(c(decomp_point_check$composition_point-decomp_point_check$composition_component,decomp_point_check$within_stratum_point-decomp_point_check$selection_component)))
destination_point_check<-dest_con_simultaneous%>%select(destination,RD_2024_minus_1999_point)%>%
  inner_join(dest_con%>%transmute(destination=as.character(destination),RD_frozen=RD_2024_minus_1999),by="destination")
destination_point_max_diff<-max(abs(destination_point_check$RD_2024_minus_1999_point-destination_point_check$RD_frozen))
simultaneous_excluding<-sort(as.character(dest_sim$destination[dest_sim$simultaneous_excludes_zero]))
simultaneous_flag_consistent<-all(dest_sim$simultaneous_excludes_zero==(dest_sim$simultaneous_95_ci_lo>0|dest_sim$simultaneous_95_ci_hi<0))
race_segment_gaps<-race_plot%>%filter(regime=="single_race")%>%group_by(race,estimate_type,display_segment)%>%
  summarise(max_year_gap=if(n()>1) max(diff(sort(year))) else 0,.groups="drop")%>%filter(max_year_gap>1)%>%nrow()
race_na_bridges<-race_plot%>%filter(regime=="single_race")%>%group_by(race,estimate_type,display_segment)%>%
  summarise(max_non_na_year_gap=if(sum(!is.na(estimate))>1) max(diff(year[!is.na(estimate)])) else 0,.groups="drop")%>%filter(max_non_na_year_gap>1)%>%nrow()
race_2021<-race%>%filter(regime=="single_race",year==2021)
race_2021_complete<-nrow(race_2021)==6&&setequal(as.character(race_2021$race),c("White","Black","AIAN","Asian","NHOPI","Multiple"))&&
  sum(race_2021$A)==annual$A_record_axis_K80[annual$year==2021]&&sum(race_2021$B)==annual$B_main_A_and_UCD_K80[annual$year==2021]
numeric_checks<-tibble(
  check=c("A and B frozen totals","gap identity","destination annual sums","destination nine categories","primary standard support","cross-disease support","decomposition residual","Kitagawa bootstrap reps","Kitagawa point estimates unchanged","Kitagawa bootstrap CIs finite","destination simultaneous point estimates unchanged","destination simultaneous CI flags consistent","destination max-t critical value common","2021 Race Recode 40 complete","race display segments have no gaps","NA points do not bridge adjacent years","source CSV per figure"),
  value=c(paste0(A_total,"/",B_total),max(abs(annual$gap_A_minus_B-(annual$A_record_axis_K80-annual$B_main_A_and_UCD_K80))),dest_sum_error,n_distinct(dest$destination),min(std_primary$effective_weight_sum),min(cross_std$effective_weight_sum),max(abs(decomp$decomposition_residual)),paste(sort(unique(decomp_uncertainty$bootstrap_reps)),collapse=";"),decomp_point_max_diff,sum(!is.finite(c(decomp_uncertainty$composition_percentile_ci_lo,decomp_uncertainty$composition_percentile_ci_hi,decomp_uncertainty$within_stratum_percentile_ci_lo,decomp_uncertainty$within_stratum_percentile_ci_hi))),destination_point_max_diff,paste(simultaneous_excluding,collapse=";"),paste(sort(unique(dest_sim$max_t_critical_95)),collapse=";"),sum(race_2021$A),race_segment_gaps,race_na_bridges,length(source_files)),
  pass=c(A_total==51084&&B_total==27514,all(annual$gap_A_minus_B==annual$A_record_axis_K80-annual$B_main_A_and_UCD_K80),dest_sum_error<1e-8,n_distinct(dest$destination)==9,min(std_primary$effective_weight_sum)>=.999999,min(cross_std$effective_weight_sum)>=.999999,max(abs(decomp$decomposition_residual))<1e-10,identical(sort(unique(decomp_uncertainty$bootstrap_reps)),100000L),decomp_point_max_diff<1e-14,all(is.finite(c(decomp_uncertainty$composition_percentile_ci_lo,decomp_uncertainty$composition_percentile_ci_hi,decomp_uncertainty$within_stratum_percentile_ci_lo,decomp_uncertainty$within_stratum_percentile_ci_hi))),destination_point_max_diff<1e-14,simultaneous_flag_consistent,length(unique(dest_sim$max_t_critical_95))==1,race_2021_complete,race_segment_gaps==0,race_na_bridges==0,length(source_files)==length(expected_stems))
)
utils::write.csv(as.data.frame(qa_matrix),file.path(qa_dir,"figure_qc.csv"),row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(as.data.frame(numeric_checks),file.path(qa_dir,"figure_qc_numeric.csv"),row.names=FALSE,fileEncoding="UTF-8")
if(!all(qa_matrix$pass)||!all(qa_matrix$svg_has_editable_text)||!all(numeric_checks$pass))stop("Figure QA failed")
writeLines(c(
  "# Figure-recreation QA",
  "",
  sprintf("- Expected figures: %d",length(expected_stems)),
  sprintf("- Export files: %d",nrow(qa_files)),
  sprintf("- Recreated source-data files: %d",length(source_files)),
  sprintf("- Export matrix PASS: %s",all(qa_matrix$pass)),
  sprintf("- Editable SVG text PASS: %s",all(qa_matrix$svg_has_editable_text)),
  sprintf("- Numeric checks PASS: %s",all(numeric_checks$pass))
),file.path(qa_dir,"figure_qc.md"),useBytes=TRUE)
make_contact_sheet<-function(pattern,outname,ncol){
  p<-list.files(out_dir,pattern=pattern,full.names=TRUE)
  imgs<-lapply(p,function(f){x<-magick::image_read(f);x<-magick::image_resize(x,"900x");magick::image_annotate(x,basename(f),size=18,color="black",boxcolor="white",gravity="north")})
  rows<-lapply(split(imgs,ceiling(seq_along(imgs)/ncol)),function(z)magick::image_append(do.call(c,z)))
  sheet<-magick::image_append(do.call(c,rows),stack=TRUE)
  tmp<-file.path(render_root,outname);magick::image_write(sheet,tmp)
  if(!file.copy(tmp,file.path(qa_dir,outname),overwrite=TRUE))stop(paste("Contact-sheet copy failed:",outname))
}
make_contact_sheet("^Figure[1-6].*png$","figure_qc_contact_sheet.png",2)
make_contact_sheet("^Supplementary.*png$","figure_qc_supp_contact_sheet.png",1)
hash_targets<-c(
  file.path(repo_root,"code","figures","revision_figures.R"),list.files(out_dir,full.names=TRUE),list.files(source_dir,full.names=TRUE),
  file.path(qa_dir,c("figure_qc.md","figure_qc.csv","figure_qc_numeric.csv","figure_qc_contact_sheet.png","figure_qc_supp_contact_sheet.png"))
)
hash_manifest<-tibble(
  path=sub(paste0("^",repo_root,"/"),"",gsub("\\\\","/",hash_targets)),bytes=file.info(hash_targets)$size,
  sha256=vapply(hash_targets,function(p)digest::digest(file=p,algo="sha256",serialize=FALSE),character(1))
)
utils::write.csv(as.data.frame(hash_manifest),file.path(qa_dir,"figure_qc_hashes.csv"),row.names=FALSE,fileEncoding="UTF-8")
cat(sprintf("Rendered %d figures, %d export files, %d source-data files. QA PASS.\n",length(expected_stems),nrow(qa_files),length(source_files)))
