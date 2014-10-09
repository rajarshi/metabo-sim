library(rcdk)
library(parallel)
library(ggplot2)
library(reshape2)

## read in data from O'Hagel paper
## SMILES generated from supplied InChI using molconvert 6.0.5
dat <- read.table('valid-smiles.txt', header=FALSE,
                  as.is=TRUE, comment='', sep='\t', quote='')
udat <- unique(dat)
mols <- parse.smiles(udat$V3)
fps <- lapply(mols, get.fingerprint, type='extended')
fps.maccs <- lapply(mols, get.fingerprint, type='maccs')
fps.circ <- lapply(mols, get.fingerprint, type='circular')

## read in a WEHI set
get.wehi.sims <- function(fname, fps=NULL, type) {
  wehi <- read.table(fname, as.is=TRUE, header=FALSE, comment='', sep=' ')
  wehi <- parse.smiles(wehi$V1)
  cat('parsed smiles\n')
  fp.wehi <- lapply(wehi, get.fingerprint, type=type)
  sims.wehi <- unlist(mclapply(fp.wehi, function(afp) {
    sims <- sapply(fps, function(x) distance(x, afp))
    return(max(sims))
  }, mc.cores=6))
  return(sims.wehi)
}

wehi <- get.wehi.sims('WEHI-RND-FAIL.sml', fps, 'extended')
wehi.m <- get.wehi.sims('WEHI-RND-FAIL.sml', fps.maccs, 'maccs')
wehi.c <- get.wehi.sims('WEHI-RND-FAIL.sml', fps.circ, 'circular')

wehi.pass <- get.wehi.sims('WEHI-RND-PASS.sml', fps, 'extended')
wehi.pass.m <- get.wehi.sims('WEHI-RND-PASS.sml', fps.maccs, 'maccs')
wehi.pass.c <- get.wehi.sims('WEHI-RND-PASS.sml', fps.circ, 'circular')

x <- rbind(data.frame(fp='MACCS (166 bit)', PAINS='Fail', y=wehi.m),
           data.frame(fp='MACCS (166 bit)', PAINS='Pass', y=wehi.pass.m),
           data.frame(fp='Path (1024 bit)', PAINS='Fail', y=wehi),
           data.frame(fp='Path (1024 bit)', PAINS='Pass', y=wehi.pass),
           data.frame(fp='ECFP6 (1024 bit)', PAINS='Fail', y=wehi.c),
           data.frame(fp='ECFP6 (1024 bit)', PAINS='Pass', y=wehi.pass.c))

ggplot(x, aes(x=PAINS, y=y))+
  geom_boxplot()+facet_wrap(~fp)+
  ylab('Tanimoto Similarity')+xlab("PAINS Status")+
  theme(strip.text.x = element_text(face='bold', size=10))

dev.copy(png, 'sim-dist.png', 600, 200);dev.off()

## Are distributions for PAINS pass/fail different?

## Lipinski Ro5 and max sim?
wehi <- read.table('WEHI-RND-PASS.sml', as.is=TRUE, header=FALSE, comment='', sep=' ')
wehi <- parse.smiles(wehi$V1)
dn <- get.desc.names()
dn <- dn[ grep('RuleOf', dn) ]
desc.wehi <- eval.desc(wehi, dn)
desc.metab <- eval.desc(mols, dn)
plot(desc.wehi[,1], wehi.pass.c)

## save.image('work.Rda')
