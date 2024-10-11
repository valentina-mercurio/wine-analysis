data<-read.csv("wine.csv",header=T)
data<-data.frame(data)
str(data)
head(data)
summary(data)
library(GGally)
library(ggplot2)
library(ggfortify)
library(dplyr)
std_data<-scale(data)
std_data<-data.frame(std_data)
ggpairs(std_data)

pca_princomp <- princomp(std_data, cor = TRUE)
summary(pca_princomp)

print(summary(pca_princomp, loadings = TRUE)) 

cor(data, pca_princomp$scores)

ggcorr(cbind(data, pca_princomp$scores), label = TRUE, cex = 2.5)

pca_sc <- data.frame(pca_princomp$scores)

plot(pca_princomp, main="Screeplot", type = "lines")

plot(pca_princomp, main="Regola di Kaiser", type = "lines")
abline(h=1)

library(rgl)
plot3d(pca_princomp$scores[, 1:3], size = 5)
play3d(spin3d(axis = c(1, 1, 1), rpm = 5), duration = 10)


###########CLUSTER###################

d <- dist(std_data)

data_sl <- hclust(d, method = "single")
data_cl <- hclust(d, method = "complete")
data_al <- hclust(d, method = "average")
data_ward <- hclust(d, method = "ward.D2")

op <- par(mfrow = c(2, 2))
plot(data_sl, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Legame Singolo)", xlab = "data")
plot(data_cl, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Legame completo)", xlab = "data")
plot(data_al, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Legame medio)", xlab = "data")
plot(data_ward, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Metodo di Ward)", xlab = "data")
par(op)

op <- par(mfrow = c(1, 2))
plot(data_cl, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Legame completo)", xlab = "data")
plot(data_ward, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Metodo di Ward)", xlab = "data")
par(op)

op <- par(mfrow = c(1, 1))
plot(data_cl, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Legame completo)", xlab = "data")
hc<-rect.hclust(data_cl, k=4)
par(op)

op <- par(mfrow = c(1, 1))
plot(data_ward, labels = c(1:nrow(std_data)), cex = .7,
     main = "Wine data (Metodo di Ward)", xlab = "data")
hc<-rect.hclust(data_cl, k=3)
par(op)


data_cl_m <- cutree(data_cl, k = 3)
data_cl_m

library(clusterSim)

minC <- 2
maxC <- 6
res <- sapply(minC:maxC,
              function(nc, data) index.G1(data, cutree(data_cl, k = nc)),
              data = std_data)
ggp <- ggplot(data=data.frame(x=seq_along(res)+1, y=res), mapping=aes(x=x,y=y)) + 
    geom_point() + 
    geom_line() +
    xlab("Numero di cluster") +
    ylab("Statistica pseudo-F di Calinski-Harabasz")
print(ggp)



data_ward_clust<- cutree(data_ward, k=3)
data_ward_clust
#silhouette con Ward
sil<-silhouette(data_ward_clust, d)
plot(sil)
sil


data_cl_clust<- cutree(data_cl, k=4)
data_cl_clust
#silhouette con legame completo
sil<-silhouette(data_cl_clust, d)
plot(sil)
sil

###########PROFILAZIONE CLUSTER
par(mfrow=c(1,1))
data_class<-cutree(data_ward, k=3)
plot(data_class)

df<-cbind(std_data, Cluster=data_ward_clust)
df<-as.data.frame(df)
library(dplyr)
df %>%
  mutate(Cluster = as.factor(Cluster))  %>%
  ggpairs(columns = 1:13, mapping = aes(color = Cluster))


dt2<-cbind(data, member=data_ward_clust)
dt2<-as.data.frame(dt2)
table(dt2$member)


by(data = dt2[,-(14)], INDICES = dt2$member, FUN = summary)

data_summ <- group_by(dt2, member) %>%
  summarise(Alcohol = mean(Alcohol),
            Malic_Acid= mean(Malic_Acid),
            Ash = mean(Ash),
            Ash_Alcanity= mean(Ash_Alcanity),
            Magnesium = mean(Magnesium),
            Total_Phenols = mean(Total_Phenols),
            Flavanoids = mean(Flavanoids),
            Nonflavanoid_Phenols = mean(Nonflavanoid_Phenols),
            Proanthocyanins = mean(Proanthocyanins),
            Color_Intensity = mean(Color_Intensity),
            Hue = mean(Hue),
            OD280 = mean(OD280),
            Proline = mean(Proline))
palette(rainbow(13))
data_summ

to.draw <- apply(data_summ[, -1], 2, function(x) x/max(x))

stars(to.draw, draw.segments = TRUE, scale = FALSE, key.loc = c(4.6, 2.3), locations=NULL,
      labels = c("CLUSTER 1", "CLUSTER 2", "CLUSTER 3"),
      main = "Wine data (profili cluster)", cex = .75,
      flip.labels = TRUE)

library(clusterSim)
library(ggplot2)

    ####################PAM###############


par(mfrow= c(1, 1))
library(psych)
library(cluster)
library(factoextra)


fviz_nbclust(std_data, pam, method ="silhouette")+theme_minimal()

fviz_nbclust(std_data, pam, method ="wss")+theme_minimal()

fviz_nbclust(std_data, pam, method ="gap_stat")+theme_minimal()

pamResult <-pam(std_data, k = 3)
pamResult
data$cluster = pamResult$cluster
head(data)
pamResult$medoids
pamResult$clustering
fviz_cluster(pamResult, 
             palette("default"),
             ellipse.type ="convex",
             repel =TRUE,
             ggtheme =theme_minimal())