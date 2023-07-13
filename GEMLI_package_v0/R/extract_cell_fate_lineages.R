extract_cell_fate_lineages<- function(GEMLI, selection, unique=FALSE, threshold)
{
  Lookup<-merge(as.data.frame(GEMLI[['predicted_lineage_table']]), GEMLI[['cell_type']], by="cell.ID", all=TRUE)
  if (unique){ 
    Lookup$cell.fate<-NA
    Lookup<-Lookup %>% group_by(clone.ID) %>% mutate(cell.fate=case_when((n_distinct(cell.type)==length(selection)& all(cell.type %in% selection)& is.na(clone.ID)==FALSE)~ "asym", (n_distinct(cell.type)==1& all(cell.type %in% selection) & n_distinct(cell.ID)>1& is.na(clone.ID)==FALSE)~"sym"))
  } else {
    Lookup2<-Lookup[Lookup$cell.type %in% selection, ]
    Lookup2<-Lookup2 %>% group_by(clone.ID) %>% mutate(cell.fate=case_when((n_distinct(cell.type)==length(selection)& all(cell.type %in% selection)& is.na(clone.ID)==FALSE)~ "asym", (n_distinct(cell.type)==1& all(cell.type %in% selection) & n_distinct(cell.ID)>1& is.na(clone.ID)==FALSE)~"sym"))
    Lookup2<-Lookup2[,c(1,4)]; Lookup<-merge(Lookup, Lookup2, by=c("cell.ID" ), all=TRUE)
  }
  # filter by threshold
  Lookup <- Lookup %>% group_by(cell.fate, clone.ID) %>% mutate(cnt = n()); Lookup<-Lookup%>% group_by(cell.fate, clone.ID, cell.type) %>% mutate(per= (n()/cnt)*100) 
  for(i in 1:length(selection)){Lookup <- Lookup %>% group_by(cell.fate, clone.ID) %>% mutate(cell.fate=case_when(((cell.fate=="asym") & (cell.type==selection[i] ) & (per < threshold[i]))~"filtered", TRUE~cell.fate))}
  Lookup<-Lookup %>% group_by(clone.ID) %>% mutate(cell.fate=case_when(any(cell.fate=="filtered")~NA, TRUE~cell.fate))
  Lookup$cell.fate <- paste(Lookup$cell.fate,Lookup$cell.type,sep = "_")
  GEMLI[['cell_fate_analysis']]<-Lookup[,1:4]
  return(GEMLI)
} 
