DEG_volcano_plot<-function(GEMLI_items, name1, name2){
DEG<-GEMLI_items[['DEG']]
DEG$change = ifelse(DEG$p_val_adj <= 0.05 & abs(DEG$avg_log2FC) >= 0.5, ifelse(DEG$avg_log2FC> 0.5 ,name1,name2),'Stable')
DEG$label=rownames(DEG)
plt<-ggplot(data = DEG, aes(x = avg_log2FC , y = -log10(p_val_adj), colour=change, label=label)) +
  geom_point(alpha=0.4, size=3.5)+
  xlim(c(-4.5, 4.5)) +
  scale_color_manual(values=c("#5386BD", "darkred","grey"))+
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8)+
  labs(x="log2(fold change)",y="-log10 (p-value)", title=paste0("DEG"," ",name1," vs ",name2)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank()) +
  geom_text_repel(data = subset(DEG, avg_log2FC >= 0.5 | avg_log2FC < -0.5), aes(label = label), max.overlaps = 15)
suppressWarnings(print(plt))
}
