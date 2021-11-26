library(ghyp)
library(fitdistrplus)
library(pastecs)
library(psych)
library(corrplot)
library(VineCopula)
library(copula)
library(corrgram)

## data load
speedlog <- read.csv('C:\\Users\\g\\Desktop\\博论书写\\论文3\\TTV\\data\\case2_ttv_data.csv')
pairs(speedlog) # 散点图

speedlog <- as.matrix(speedlog)

## -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* statistical analysis-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
## statistics
summary(speedlog)
stat.desc(speedlog, norm = TRUE)
describe(speedlog)

## method 1： correlations plot
corpearson <- cor(speedlog, method = 'pearson')
corrplot(corpearson, method = 'number', number.cex = 0.8, diag = TRUE,tl.cex = 0.8)
corrplot(corpearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')

corspearman <- cor(speedlog, method = 'spearman')
corrplot(corspearman, method = 'number', number.cex = 0.8, diag = TRUE,tl.cex = 0.8)
corrplot(corspearman, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')

corkendall <- cor(speedlog[1:62,], method = 'kendall')
corrplot(corkendall, method = 'number', number.cex = 0.8, diag = TRUE,tl.cex = 0.8)
corrplot(corkendall, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')
write.csv(round(corkendall,3),'C:\\Users\\g\\Desktop\\博论书写\\论文3\\TTV\\data\\case2_cor1.csv', row.names = FALSE)

## method 2： correlations plot
# pts:scatters; ellipse:平滑曲线和置信椭圆 shade:阴影; 
# pie; bar; cor; conf; fill
corrgram(speedlog[1:62,], 
         cor.method='pearson', 
         order = FALSE, 
         lower.panel = panel.cor, 
         upper.panel = panel.pts,
         text.panel = panel.txt,
         col.regions = colorRampPalette(c("darkgoldenrod4", "burlywood1",
                                          "darkkhaki", "darkgreen")),
         diag.panel = NULL) 

## -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* distribution fitting and testing -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

dat <- speedlog[1:62,12]
#dat_r <- speedlog[62:92,12]

## GH distribution
fg <- fit.ghypuv(dat, lambda = 0.5, alpha.bar = 0.5, mu = mean(dat), sigma = mad(dat), gamma = 0.5)
summary(fg)

## T distribution
ft <- fit.tuv(dat, nu = 20)
#ft2 <- fitdistr(dat,'t',df=20)
summary(ft)
ks.test(dat,'pt',df=2)

## Gaussian distribution
fn <- fit.gaussuv(dat, na.rm = T, save.data = T)
summary(fn)
ks.test(dat_r,"pnorm",-0.0853,	0.1902)

## -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* Copula fitting and GoF test -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

M1 <- as.matrix(speedlog[1:62,])
u1 <- pobs(M1) # Pseudo observations, convert to rank,and then normalization

param_len <- ncol(u1)*(ncol(u1)-1)/2
df <- 10
## nornalcopula
nc.ml <- fitCopula(normalCopula(dim=ncol(u1), dispstr="un"), 
                   u1, 
                   method="ml", 
                   optim.method="Nelder-Mead",
                   estimate.variance=TRUE,
                   start=c(rep(0,param_len))) # start = [n(n-1)/2]个参数
summary(nc.ml)
gofCopula(normalCopula(dim = ncol(u1)), u1, N = 50, ties=TRUE)

## t_Copula
tc.mpl <- fitCopula(tCopula(dim=ncol(u1), dispstr="un"),
                    u1, 
                    method="ml", 
                    estimate.variance=TRUE,
                    start=c(rep(0,param_len),df)) # start = (rho[1:3], df)
summary(tc.mpl)
gofCopula(tCopula(dim=ncol(u1), df=5.124), method= 'SnC', u1, N = 50, ties=TRUE)

## frankcopula
fc.mpl <- fitCopula(frankCopula(dim=ncol(u1)),u1, method="ml")
summary(fc.mpl)
gofCopula(frankCopula(dim = ncol(u1)), u1, N = 1000, ties=TRUE)

## claytoncopula
cc.mpl <- fitCopula(claytonCopula(dim=ncol(u1)),u1, method="ml")
summary(cc.mpl)
gofCopula(claytonCopula(dim = ncol(u1)), u1, N = 50, ties=TRUE)

## gumbelcopula
gc.mpl <- fitCopula(gumbelCopula(dim=ncol(u1)),u1, method="ml")
summary(gc.mpl)
gofCopula(gumbelCopula(dim = ncol(u1)), u1, N = 50, ties=TRUE)

## joecopula
jc.mpl <- fitCopula(joeCopula(dim=ncol(u1)),u1, method="ml")
summary(jc.mpl)
gofCopula(joeCopula(dim = ncol(u1)), u1, N = 50, ties=TRUE)


summary(nc.ml)
summary(tc.mpl)
summary(fc.mpl)
summary(cc.mpl)
summary(gc.mpl)
summary(jc.mpl)

## -*-*-*-*-*-*-*-*-*-*-* Correlation matrix from optim copula -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

## ralation matrix 
cor_cop <- matrix(c(rep(1,144)),ncol = 12)

# rewrite the up-triangular matrix
param_ind <- 1
for (i in 1:11) { # row
  for (j in (i+1):12) { # col
    cor_cop[i,j] <- tc.mpl@copula@parameters[param_ind]
    param_ind <- param_ind+1
  }
}
# rewrite the low-triangular matrix
param_ind <- 1
for (i in 1:11) { # col
  for (j in (i+1):12) {  # row
    cor_cop[j,i] <- tc.mpl@copula@parameters[param_ind]
    param_ind <- param_ind+1
  }
}

write.csv(round(cor_cop,3),'C:\\Users\\g\\Desktop\\博论书写\\论文3\\TTV\\data\\case2_cor2.csv', row.names = FALSE)
corrplot(cor_cop, method = 'number', number.cex = 0.8, diag = TRUE,tl.cex = 0.8)
corrplot(cor_cop, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')
