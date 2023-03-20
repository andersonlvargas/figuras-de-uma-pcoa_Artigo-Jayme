

# Para permanova
# para abundancia do dung:
# PARA FIGURAS DO PERMANOVA, UM NOME PARA CADA ARQUIVO!

require(vegan)     #para carregar o pacote vegan
datazoo<-read.table(file.choose(),header=T)      # para carregar os dados do zoo
str(datazoo)   # para ver um resumo dos dados
datazoo       # para ver todos os dados

factor(datazoo$treatment)->datazoo$treatment   # treatment é a variável preditora
is.factor(datazoo$treatment)
#se houver 2 variáveis preditora eu faço isso 2 vezes.

bio<-datazoo[-c(1)]   # esse comando exclui a 1ª coluna da análise, pois a  # 1ª  coluna são as variáveis preditivas. O bio é a tabela sem a 1ª coluna.
bio
sqrt(bio)->bio
adonis(bio ~ treatment ,datazoo, perm=9999, method ="bray")


#betadisper
https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/betadisper

datazoo
## Bray-Curtis distances between samples
dis <- vegdist(bio)
dis
## Calculate multivariate dispersions
mod1 <- betadisper(dis, datazoo$treatment)
mod1

## Perform test
anova(mod1)

## Permutation test for F
permutest(mod1, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod1))
plot(mod.HSD)
plot(mod1)


########################################

#figura abundance low dung sample 1 - mod1:
par(mfrow=c(2,2))
#treat=c(rep("Treatment1",12),rep("Treatment2",12),rep("Treatment3",12))
ordiplot(mod1,type="n", xlim=c(-0.3, 0.3))  # 0.3 para low dung
ordihull(mod1,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
#orditorp(mod1,display="centroids",col="red",air=0.01)
orditorp(mod1,display="centroids",col=c(rep("black",6),rep("black",6)),
   air=0.01,cex=1)
legend("topleft", "a)", bty="n") 

#figura abundance high dung sample 1 - mod2:
ordiplot(mod2,type="n", xlim=c(-0.4, 0.3))  # 0.3 para low dung
ordihull(mod2,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod2,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "b)", bty="n") 

#figura abundance low dung sample 2 - mod3:
ordiplot(mod3,type="n", xlim=c(-0.4, 0.4))  # 0.3 para low dung
ordihull(mod3,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod3,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "c)", bty="n") 

#figura abundance high dung sample 2 - mod4:
ordiplot(mod4,type="n", xlim=c(-0.4, 0.3))  # 0.3 para low dung
ordihull(mod4,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod4,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "d)", bty="n") 



# exportar figura?
tiff(filename="C:/Users/jayme/Desktop/figure.abund.juntos.tif", height=1500, width=1500, bg="white", res = 300)
par(mfrow=c(2,2))
#treat=c(rep("Treatment1",12),rep("Treatment2",12),rep("Treatment3",12))
ordiplot(mod1,type="n", xlim=c(-0.3, 0.3))  # 0.3 para low dung
ordihull(mod1,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
#orditorp(mod1,display="centroids",col="red",air=0.01)
orditorp(mod1,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "a)", bty="n") 

#figura abundance high dung sample 1 - mod2:
ordiplot(mod2,type="n", xlim=c(-0.4, 0.3))  # 0.3 para low dung
ordihull(mod2,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod2,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "b)", bty="n") 

#figura abundance low dung sample 2 - mod3:
ordiplot(mod3,type="n", xlim=c(-0.4, 0.4))  # 0.3 para low dung
ordihull(mod3,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod3,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "c)", bty="n") 

#figura abundance high dung sample 2 - mod4:
ordiplot(mod4,type="n", xlim=c(-0.4, 0.3))  # 0.3 para low dung
ordihull(mod4,groups=datazoo$treatment,draw="polygon",col="grey90",label=F)
orditorp(mod4,display="centroids",col=c(rep("black",5),rep("black",5)),
   air=0.01,cex=1)
legend("topleft", "d)", bty="n") 

dev.off()

##############  em 16/03/2022
########                   #############3
###########3##########################
require(vegan)

#para abundância. Colar os arquivos no desktop
ab1ld <- read.table("C:/Users/jayme/Desktop/dung abundance sample 1 low dung.txt", header = T)
ab1hd <- read.table("C:/Users/jayme/Desktop/dung abundance sample 1 high dung.txt", header = T)
ab2ld <- read.table("C:/Users/jayme/Desktop/dung abundance sample 2 low dung.txt", header = T)
ab2hd <- read.table("C:/Users/jayme/Desktop/dung abundance sample 2 high dung.txt", header = T)

#para biomassa:
ab1ld <- read.table("C:/Users/jayme/Desktop/dung biomass sample 1 low dung.txt", header = T)
ab1hd <- read.table("C:/Users/jayme/Desktop/dung biomass sample 1 high dung.txt", header = T)
ab2ld <- read.table("C:/Users/jayme/Desktop/dung biomass sample 2 low dung.txt", header = T)
ab2hd <- read.table("C:/Users/jayme/Desktop/dung biomass sample 2 high dung.txt", header = T)



#(1) ab1ld  , abundance, sampling 1, low dung
factor(ab1ld$treatment)->ab1ld$treatment   # treatment é a variável preditora
is.factor(ab1ld$treatment)

#(2) ab1hd
factor(ab1hd$treatment)->ab1hd$treatment   # treatment é a variável preditora
is.factor(ab1hd$treatment)

#(3) ab2ld
factor(ab2ld$treatment)->ab2ld$treatment   # treatment é a variável preditora
is.factor(ab2ld$treatment)

#(4) ab2hd
factor(ab2hd$treatment)->ab2hd$treatment   # treatment é a variável preditora
is.factor(ab2hd$treatment)

#(1)
bio1<-ab1ld[-c(1)]   # esse comando exclui a 1ª coluna da análise, pois a  # 1ª  coluna são as variáveis preditivas. O bio é a tabela sem a 1ª coluna.
bio1
sqrt(bio1)->bio1
adonis(bio1 ~ treatment ,ab1ld, perm=9999, method ="bray")

#(2)
bio2<-ab1hd[-c(1)]   # esse comando exclui a 1ª coluna da análise, pois a  # 1ª  coluna são as variáveis preditivas. O bio é a tabela sem a 1ª coluna.
bio2
sqrt(bio2)->bio2
adonis(bio2 ~ treatment ,ab1hd, perm=9999, method ="bray")

#(3)
bio3<-ab2ld[-c(1)]   # esse comando exclui a 1ª coluna da análise, pois a  # 1ª  coluna são as variáveis preditivas. O bio é a tabela sem a 1ª coluna.
bio3
sqrt(bio3)->bio3
adonis(bio3 ~ treatment ,ab2ld, perm=9999, method ="bray")

#(4)
bio4<-ab2hd[-c(1)]   # esse comando exclui a 1ª coluna da análise, pois a  # 1ª  coluna são as variáveis preditivas. O bio é a tabela sem a 1ª coluna.
bio4
sqrt(bio4)->bio4
adonis(bio4 ~ treatment ,ab2hd, perm=9999, method ="bray")


dis.ab1ld <- vegdist(bio1)
dis.ab1hd <- vegdist(bio2)
dis.ab2ld <- vegdist(bio3)
dis.ab2hd <- vegdist(bio4)

## Calculate multivariate dispersions
mod1 <- betadisper(dis.ab1ld, ab1ld$treatment)
mod2 <- betadisper(dis.ab1hd, ab1hd$treatment)
mod3 <- betadisper(dis.ab2ld, ab2ld$treatment)
mod4 <- betadisper(dis.ab2hd, ab2hd$treatment)

plot(mod4)


#esse funciona com cores: colocar os 4 grafs
tiff(filename="C:/Users/jayme/Desktop/figure.abund.juntos3.tif", height=1500, width=1500, bg="white", res = 300)
#par(mfrow=c(2,2))
par(mar=c(5.1,4,1,1), mfrow=c(2,2))
#(1)low dung, abundance, sampling 1
ordiplot(mod1,type="n", xlim=c(-0.3, 0.3), ylim=c(-0.5,0.5))  # 0.3 para low dung
ordihull(mod1,groups=ab1ld$treatment,draw="polygon",col=c("lightpink2","aquamarine","green","gold","coral","gray"),label=F,
         border=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"))
orditorp(mod1,display="centroids",col=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"),
         air=0.01,cex=0.9)
legend("topleft", "a) Low dung, sampling 1", bty="n", cex=0.7)

#(2)high dung, abundance, sampling 1
ordiplot(mod2,type="n", xlim=c(-0.3, 0.3), ylim=c(-0.5,0.5))  # 0.3 para low dung
ordihull(mod2,groups=ab1hd$treatment,draw="polygon",col=c("lightpink2","aquamarine","green","gold","coral","gray"),label=F,
         border=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"))
orditorp(mod2,display="centroids",col=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"),
         air=0.01,cex=0.9)
legend("topleft", "b) High dung, sampling 1", bty="n", cex=0.7)

#(3)
ordiplot(mod3,type="n", xlim=c(-0.3, 0.3), ylim=c(-0.5,0.5))  # 0.3 para low dung
ordihull(mod3,groups=ab2ld$treatment,draw="polygon",col=c("lightpink2","aquamarine","green","gold","coral","gray"),label=F,
         border=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"))
orditorp(mod3,display="centroids",col=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"),
         air=0.01,cex=0.9)
legend("topleft", "c) Low dung, sampling 2", bty="n", cex=0.7)


#(4)
ordiplot(mod4,type="n", xlim=c(-0.3, 0.3), ylim=c(-0.5,0.5))  # 0.3 para low dung
ordihull(mod4,groups=ab2hd$treatment,draw="polygon",col=c("lightpink2","aquamarine","green","gold","coral","gray"),label=F,
         border=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"))
orditorp(mod4,display="centroids",col=c("lightpink4","blue","darkgreen","darkorange2","darkred","gray37"),
         air=0.01,cex=0.9)
legend("topleft", "d) High dung, sampling 2", bty="n", cex=0.7)

dev.off()

