check_vars <- function(x = '.') {
  R_files <- list.files(path = x, pattern = 'R$|Rmd$', 
                            recursive = TRUE, full.names = TRUE)
  
  all_objects <- ls(.GlobalEnv)
  
  
  scan.file <- function(x){
    test.files <- function(y){
      file_text <- readLines(y)
      any(regexpr(x, file_text, fixed = TRUE) > 0, 
          na.rm = TRUE)
    }
    
    size <- as.numeric(object.size(eval(as.symbol(x))))
    counts <- as.numeric(sapply(R_files, test.files))
    
    output <- data.frame(size = size, 
                         t(counts),
                         sum = sum(counts))
    colnames(output)[-c(1, ncol(output))] <- R_files
    output
  }
  
  file.test <- sapply(all_objects, scan.file)
  
  write.csv(t(file.test), 'variables.csv')
  
  file.test

}