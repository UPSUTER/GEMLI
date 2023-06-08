cell_type_composition_plot <- function(GEMLI_items, ground_truth=F, cell_type_colors=F, type, intersections=NULL)
{
  base_colors = rep(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'), 100)
  if (cell_type_colors==F){
    cell.type<-unique(GEMLI_items[['cell_type']]$cell.type)
    color<-base_colors[rank(unique(GEMLI_items[['cell_type']]$cell.type))]
    GEMLI_items[['cell_type_color']] = data.frame(cell.type, color)} 
  
 if (ground_truth){cell.ID<-names(GEMLI_items[['barcodes']]); clone.ID<-unname(GEMLI_items[['barcodes']]); GT<-as.data.frame(cbind(clone.ID,cell.ID));
  Lookup<-merge(GT, GEMLI_items[['cell_type']], by="cell.ID", all=TRUE)} else {
    Lookup<-merge(as.data.frame(GEMLI_items[['predicted_lineage_table']]), GEMLI_items[['cell_type']], by="cell.ID", all=TRUE)
    }
  
  if (type == "bubble"){
    Lookup <- Lookup %>% group_by(clone.ID, cell.type) %>% summarise(cnt = n()) %>% mutate(freq = round(cnt / sum(cnt), 3)); Lookup <- reshape::cast(Lookup, clone.ID~cell.type, value="freq"); 
    base_colors = GEMLI_items[['cell_type_color']]$color[match(colnames(Lookup[,2:length(Lookup)]),GEMLI_items[['cell_type_color']]$cell.type)]
    p<-Lookup %>% gather(cell.type, percentage, -clone.ID)%>% ggplot(group=cell.type) + geom_point(aes(x = cell.type, y = clone.ID, size = percentage, col= cell.type))+ theme_classic()+ scale_colour_manual(values = base_colors)}
  
  if (type == "upsetR"){
    Lookup_list <- split(Lookup$clone.ID, Lookup$cell.type)
    p<-upset(fromList(Lookup_list), order.by = "freq", nsets = length(Lookup_list), 
          sets.x.label = "Lineages in cell type", mainbar.y.label = "Number of lineages",
          nintersects = NA, intersections = NULL, point.size=5, mb.ratio = c(0.5, 0.5), text.scale = 2,
          set_size.show = TRUE, set_size.numbers_size = 7, set_size.scale_max = length(unique(Lookup$clone.ID)))}

  if (type == "plain"){
    Lookup<-unique(Lookup[,-c(1)])
    p<-Lookup %>% group_by(clone.ID) %>% arrange(clone.ID, cell.type) %>% summarize(combi = paste0(cell.type, collapse = "__"), .groups = "drop") %>% count(combi)}
  
  return(p)
  }
