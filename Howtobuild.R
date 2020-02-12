#How to build and install package

getwd()

#install roxygen2
library(roxygen2)
#set namespaces
roxygen2::roxygenise()

#get compiler
install.packages(devtools)
library(devtools)

build()

install()
