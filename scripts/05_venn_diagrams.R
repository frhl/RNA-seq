
files = list.files('derived/', pattern = '\\.csv', full.names = T)[c(1,3)]

# venn data
H441_SUBvsKO_SUB <- read.csv("derived//201125_H441_SUBvsKO_SUB.csv")
H441_SUBvsKO_SUB$hgnc_symbol[H441_SUBvsKO_SUB$hgnc_symbol == ''] <- NA
H441_SUBvsKO_SUB$significant <- H441_SUBvsKO_SUB$FDR < 0.05
SALIvsKO_SALI <- read.csv("derived//201125_SALIvsKO_SALI.csv" )
SALIvsKO_SALI$hgnc_symbol[SALIvsKO_SALI$hgnc_symbol == ''] <- NA
SALIvsKO_SALI$significant <- SALIvsKO_SALI$FDR < 0.05

# setup venn diagram
venn <- list(H441_SUBvsKO_SUB = as.vector(na.omit(H441_SUBvsKO_SUB$hgnc_symbol[H441_SUBvsKO_SUB$significant])),
             SALIvsKO_SALI = as.vector(na.omit(SALIvsKO_SALI$hgnc_symbol[SALIvsKO_SALI$significant])))

# remove things that are NA
library(VennDiagram)
v <- venn.diagram(venn,
                  fill = c("red", "green"),
                  alpha = c(0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  height = 480, 
                  width = 480, 
                  filename=NULL)
grid.newpage()
grid.draw(v)



