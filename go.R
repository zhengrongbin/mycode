library(gdata) 
library(WriteXLS)

setwd('/Users/zhengrongbin/Documents/gallbladder/data/snoRNA_correlation/')

rna12_bp = read.xls('snoRNA12_bp.xls')
rna12_cc = read.xls('snoRNA12_cc.xls')
rna12_mf = read.xls('snoRNA12_mf.xls')
rna12_all = rbind(rna12_bp[,c('Pvalue', 'Term')], rna12_cc[,c('Pvalue', 'Term')], rna12_mf[,c('Pvalue', 'Term')])
rna12_all = rna12_all[order(rna12_all[,1]),]
top15 = rna12_all[1:15,]
pdf('../../plot/RNA12_GO.pdf', width = 10, height = 6)
attach(mtcars)
opar= par(no.readonly=TRUE)
par(mar=c(6,25,4,4))
ggplot(top15, aes(x = reorder(Term, log10(Pvalue)), y = -log10(Pvalue)))+geom_bar(stat="identity", fill = "#B22222")+
ylab('-log10(p.value)')+xlab("")+coord_flip()+theme_bw()+
theme(panel.border = element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), 
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"), 
    axis.title = element_text(size = 15, face = "bold"))
par(opar)
detach(mtcars)
dev.off()

rna21_bp = read.xls('snoRNA21_bp.xls')
rna21_cc = read.xls('snoRNA21_cc.xls')
rna21_mf = read.xls('snoRNA21_mf.xls')
rna21_all = rbind(rna21_bp[,c('Pvalue', 'Term')], rna21_cc[,c('Pvalue', 'Term')], rna21_mf[,c('Pvalue', 'Term')])
rna21_all = rna21_all[order(rna21_all[,1]),]
top15.21 = rna21_all[1:15,]
pdf('../../plot/RNA21_GO.pdf', width = 10, height = 6)
attach(mtcars)
opar= par(no.readonly=TRUE)
par(mar=c(6,25,4,4))
p=ggplot(top15.21, aes(x = reorder(Term, log10(Pvalue)), y = -log10(Pvalue)))+geom_bar(stat="identity", fill = "#B22222")+
  ylab('-log10(p.value)')+xlab("")+coord_flip()+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"))
print(p)
par(opar)
detach(mtcars)
dev.off()
