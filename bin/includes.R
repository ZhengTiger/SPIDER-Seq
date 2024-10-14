library(rgl)
library(wholebrain)

path.bin <- 'F:/1.课题/1.神经环路/20230713/Figure/bin'
path.matrices <-'F:/1.课题/1.神经环路/20230713/Figure/data/3D_data'
allen.annot.path <- paste(path.matrices , 'allen-mesh', sep = '/')

load(paste(path.matrices ,'atlasspots.RData',sep='/'))
load(paste(path.matrices , 'vivid-colors.RData', sep='/'))

source(paste(path.bin,'allenAnnotationsFunctions.R',sep='/'))