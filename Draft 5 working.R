#eda 
library("faraway")
protein_train <- read.csv("protein-train.csv")

# check for NAs
na.list = sapply(protein_train, function(x){sum(is.na(x))})
na.list[na.list>0]

# correlation with accuracy

cor.list = sapply(protein_train[-1], function(x){cor(protein_train$accuracy,x)})

cor.tab  = data.frame(var = names(cor.list), r=cor.list,row.names = NULL)

cor.tab2 = subset(cor.tab, abs(r)>0.35)

vars = cor.tab2$var

vars = as.character(vars)
protein = protein_train[ c("accuracy",as.character(vars))]
#print(protein)

# summary statistics
summary(protein)



counts = protein%>% dplyr::count(carbonylC_aliph1HC_medlong)

prot.model = lm(accuracy~., data = protein)
summary(prot.model)

library(ggplot2)
library(gridExtra)

counts = list()
hplots = list()

for(i in 1:20){
  
  counts[[i]] = protein %>% dplyr::count(get(vars[i]))
  colnames(counts[[i]]) = c(vars[i],"n")
  p      = ggplot2::ggplot(data = counts[[i]], aes_string(vars[i])) + 
    ggplot2::geom_bar(aes(weight = n)) +  
    labs(title= vars[i],y="", x = "") + 
    theme(title = element_text(size = 9), plot.margin=unit(c(0,-0.5,0,-0.5),"cm"))
  
  
  hplots[[i]] = p
}

# arrange the plots in a grid
do.call("grid.arrange", c(hplots, ncol=5, nrow = 4))

counts.1 = protein %>% dplyr::count(bbC_bbC_medshort)
ggplot2::ggplot(data = counts.1, aes(bbC_bbC_medshort)) + 
  ggplot2::geom_bar(aes(weight = n))

sort(unique(protein$bbC_bbC_medshort))


# Scatter plots
splots = list()

for(i in 1:20){
  
  p   = ggplot2::ggplot(data = protein, aes(y = accuracy)) + 
    ggplot2::geom_point(aes_string(vars[i])) +  
    labs(title= vars[i],y="", x = "") + 
    theme(title = element_text(size = 9), plot.margin=unit(c(0,-0.5,0,-0.5),"cm"))
  
  
  splots[[i]] = p
}
do.call("grid.arrange", c(splots, ncol=5, nrow = 4))


prot.model = lm(accuracy~., data = protein)
summary(prot.model)


# Multicollinearity

vfit = vif(prot.model)
max(vfit)


# Homoscedasticity

library(lmtest)
lmtest::bptest(prot.model)


#residual plot
library(car)
residualPlots(model = prot.model)

# residula plot with fitted values of accuracy

plot(x = prot.model$fitted.values, y = prot.model$residuals,
     ylab = "residuals", xlab = "fitted")
abline(h =mean(prot.model$residuals), col = "red")

# update model

prot.model2 = lm(accuracy~.-aliph3HC_bbProN_vlong-aliph2HC_aliph3HC_medshort-aliph1HC_hydroxylO_long-aliph2HC_aliph3HC_medlong, data = protein)
summary(prot.model2)
#removed due to insignificant p-values

#residual plot for prot.model2

residualPlots(model = prot.model2)

# residula plot with fitted values of accuracy

plot(x = prot.model2$fitted.values, y = prot.model2$residuals,
     ylab = "residuals", xlab = "fitted")
abline(h =mean(prot.model2$residuals), col = "red")


# Multicollinearity for prot.model2

vfit2 = vif(prot.model2)
max(vfit2)


# Homoscedasticity for prot.model2

lmtest::bptest(prot.model2)

