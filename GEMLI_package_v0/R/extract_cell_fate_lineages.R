extract_cell_fate_lineages<- function(GEMLI_items, selection, unique=FALSE, threshold)
{
  Lookup<-merge(as.data.frame(GEMLI_items[['predicted_lineage_table']]), GEMLI_items[['cell_type']], by="cell.ID", all=TRUE); Lookup$cell.fate<-NA
  if (unique){ 
  Lookup<-Lookup %>% group_by(clone.ID) %>% mutate(cell.fate=case_when((n_distinct(cell.type)==length(selection) & (cell.type %in% selection))~ "asym", (n_distinct(cell.type)==1& cell.type %in% selection & n_distinct(cell.ID)>1)~"sym"))} else {
  Lookup<-Lookup %>% group_by(clone.ID) %>% mutate(cell.fate=case_when((n_distinct(cell.type)>=length(selection) & (cell.type %in% selection))~"asym",(n_distinct(cell.type)>=1 & (cell.type %in% selection) & (cell.fate %in% NA))~"sym"))
  Lookup<-Lookup %>% group_by(cell.fate, clone.ID) %>% mutate(cell.fate=case_when((n_distinct(cell.ID)==1)~NA, TRUE~cell.fate))
  }
  # filter by threshold
  Lookup <- Lookup %>% group_by(cell.fate, clone.ID) %>% mutate(cnt = n()); Lookup<-Lookup%>% group_by(cell.fate, clone.ID, cell.type) %>% mutate(per= (n()/cnt)*100) 
  for(i in 1:length(selection)){Lookup <- Lookup %>% group_by(cell.fate, clone.ID) %>% mutate(cell.fate=case_when(((cell.fate=="asym") & (cell.type==selection[i] ) & (per < threshold[i]))~"filtered", TRUE~cell.fate))}
  Lookup<-Lookup %>% group_by(clone.ID) %>% mutate(cell.fate=case_when(any(cell.fate=="filtered")~NA, TRUE~cell.fate))
  Lookup$cell.fate <- paste(Lookup$cell.fate,Lookup$cell.type,sep = "_")
  GEMLI_items[['cell_fate_analysis']]<-Lookup[,1:4]
  return(GEMLI_items)
} 
