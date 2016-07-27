all.R.files <- list.files(path = '.', pattern = 'R$|Rmd$', 
                          recursive = TRUE, full.names = TRUE)

all.objects <- ls()

scan.file <- function(x){
  test.files <- function(y){
    aa <- scan(y, what = 'character')
    any(regexpr(x, aa, fixed=TRUE)>1, na.rm=TRUE)
  }
  
  size <- as.numeric(object.size(eval(as.symbol(x))))
  counts <- as.numeric(sapply(all.R.files, test.files))
  
  output <- data.frame(size = size, 
                       t(counts),
                       sum = sum(counts))
  colnames(output)[2:21] <- all.R.files
  output
}

file.test <- sapply(all.objects, scan.file)

write.csv(t(file.test), 'variables.csv')
