library(sveval)
library(magrittr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

sr = readSVvcf(args[1], sample.name=NULL)
lr = readSVvcf(args[2], sample.name=NULL)
region_bed=args[3]
sample=args[4]

sr$ac = 1
lr$ac = 1

ol = sveval::svevalOl(calls.gr=sr, truth.gr=lr, min.size=50, method="bipartite",
							max.ins.dist=100, min.ol=.3)

SV_overlap.df <- NULL

SV_overlap.df %<>% rbind(c("DEL", "all", length(ol$svs$DEL$FN), length(ol$svs$DEL$TP), length(ol$svs$DEL$FP)))
SV_overlap.df %<>% rbind(c("INS", "all", length(ol$svs$INS$FN), length(ol$svs$INS$TP), length(ol$svs$INS$FP)))
SV_overlap.df %<>% rbind(c("INV", "all", length(ol$svs$INV$FN), length(ol$svs$INV$TP), length(ol$svs$INV$FP)))
SV_overlap.df %<>% rbind(c("DUP", "all", length(ol$svs$DUP$FN), length(ol$svs$DUP$TP), length(ol$svs$DUP$FP)))

norepeat_ol = sveval::svevalOl(calls.gr=sr, truth.gr=lr, min.size=50, method="bipartite",
							max.ins.dist=100, min.ol=.3, bed.regions=region_bed)

SV_overlap.df %<>% rbind(c("DEL", "not_difficult", length(norepeat_ol$svs$DEL$FN), length(norepeat_ol$svs$DEL$TP), length(norepeat_ol$svs$DEL$FP)))
SV_overlap.df %<>% rbind(c("INS", "not_difficult", length(norepeat_ol$svs$INS$FN), length(norepeat_ol$svs$INS$TP), length(norepeat_ol$svs$INS$FP)))
SV_overlap.df %<>% rbind(c("INV", "not_difficult", length(norepeat_ol$svs$INV$FN), length(norepeat_ol$svs$INV$TP), length(norepeat_ol$svs$INV$FP)))
SV_overlap.df %<>% rbind(c("DUP", "not_difficult", length(norepeat_ol$svs$DUP$FN), length(norepeat_ol$svs$DUP$TP), length(norepeat_ol$svs$DUP$FP)))

SV_overlap.df %<>% as.data.frame
colnames(SV_overlap.df) <- c("SV_type", "region", "LR_unique", "Shared", "SR_unique")
SV_overlap.df %<>% mutate(LR_unique = as.integer(as.character(LR_unique)), 
						  Shared = as.integer(as.character(Shared)),
						  SR_unique = as.integer(as.character(SR_unique)))
SV_overlap.df %<>% mutate(SR_overlap = Shared / (Shared + SR_unique),
						  LR_overlap = Shared / (Shared + LR_unique),
						  sample = sample)
out.file=args[5]
write.table(SV_overlap.df, file=out.file, quote=F, row.names=F, col.names=T, sep="\t")


