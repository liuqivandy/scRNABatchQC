lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

getFtpFilenames <- function(curl) {
  filenames<-getURL(curl, dirlistonly=TRUE)
  tmpcon <- textConnection(filenames, "r")
  b <- read.table(tmpcon)
  close(tmpcon)
  filenames<-as.character(b[, ncol(b)])
}

deleteFileAndMd5<-function(localfile){
  if(file.exists(localfile)){
    unlink(localfile)
    
    md5file = paste0(localfile, ".md5")
    if(file.exists(md5file)){
      unlink(md5file)
    }
  }
}

deleteFilesAndMd5<-function(localfiles){
  for(localfile in localfiles){
    deleteFileAndMd5(localfile)
  }
}

##' arrayExpressDownload
##'
##' The function downloads array express raw data from EBI ftp server based on 
##' datasets user provided. Once the compressed raw data is downloaded, 
##' individual target file will be extracted from compressed raw data. The 
##' dataset/count table will be returned.
##'
##' @param datasets the dataset names, for example: c("E-TABM-43", "E-TABM-158") 
##' @param targetDir the target directory to store the datasets
##' @param filePattern the file pattern of the expected data file extracted 
##' from gzipped file
##' @param unzip the path to the command to be used in unzip function
##' @param overwrite If TRUE, overwrite existing files, otherwise ignore such 
##' files.
##' @return a data frame containing dataset and how many expected data files in 
##' that dataset 
##' @importFrom RCurl getURL
##' @export
##' @examples 
##' #download three datasets from ArrayExpress website
##' rootDir<-paste0(dirname(tempdir()), "/DupChecker")
##' dir.create(rootDir, showWarnings = FALSE)
##' datatable<-arrayExpressDownload(datasets = c("E-MEXP-3872"), targetDir=rootDir, filePattern="cel$")
arrayExpressDownload<-function(datasets, 
                               targetDir = getwd(), 
                               filePattern=NULL, 
                               unzip="internal", 
                               overwrite=FALSE){
  targetDir<-gsub("[/\\]$","",targetDir)
  counts<-data.frame(dataset = datasets, count = rep(0, length(datasets)))
  for(dataset in datasets){
    dname<-unlist(strsplit(dataset, '-'))[2]
    if(is.na(dname)){
      stop(paste0("Wrong array express dataset name : ", dataset))
    }
    
    subdir<-file.path(targetDir, dataset)
    message(paste0("Dataset ", dataset, " is downloading to ", subdir, " ..."))
    
    dir.create(subdir, showWarnings = FALSE)
    
    curl<-paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/"
                 , dname, "/", dataset, "/")
    
    filenames<-getFtpFilenames(curl)
    
    for(filename in filenames){
      link<-paste0(curl, filename)
      localfile<-paste0(subdir, "/", filename)
      
      if(overwrite){
        deleteFileAndMd5(localfile)
      }
      
      if(!file.exists(localfile)){
        message("Downloading file ", link, " ...")
        download.file(link, localfile, method="auto", mode="wb")
        Sys.sleep(2)
      }
      
      if(grepl(".zip$", localfile)){
        expectFiles<-paste0(subdir, "/", unzip(localfile, unzip=unzip, list=TRUE)$Name)
        if(overwrite){
          deleteFilesAndMd5(expectFiles)
        }
        
        if(all(file.exists(expectFiles))){
          next
        }
        
        message("Decompressing file ", localfile, " ...")
        unzip(localfile, overwrite=overwrite, exdir=subdir, setTimes=TRUE)
        Sys.sleep(2)
      }
    }
    
    celfiles<-list.files(subdir, filePattern, ignore.case=TRUE)
    counts$count[counts$dataset==dataset] = length(celfiles)
    message("Dataset ", dataset, " has been downloaded.")
  }
  return(counts)
  
}

##' geoDownload
##'
##' The function downloads GEO raw data from ncbi ftp server based on datasets 
##' user provided. Once the compressed raw data is downloaded, 
##' individual gzipped target file will be extracted from compressed raw data, 
##' and individual target file will be extracted from corresponding 
##' gzipped file. The dataset/count table will be returned.
##'
##' @param datasets the GEO dataset names, for example: c("GSE14333") 
##' @param targetDir the target directory to store the datasets
##' @param filePattern the file pattern of the expected data file may or may not 
##'        extracted from gzipped file, for example: "cel$" for AffyMetrix
##'        CEL files. Default is NULL.
##' @param tar the path to the command to be used in untar function
##' @param overwrite If TRUE, overwrite existing files, otherwise ignore such 
##'        files. The equivalent of unzip -o.
##' @return a data frame containing dataset and how many target files in 
##'         that dataset 
##' @importFrom RCurl getURL
##' @importFrom R.utils gunzip
##' @export
##' @examples 
##' #download three datasets from GEO website
##' rootDir<-paste0(dirname(tempdir()), "/DupChecker")
##' dir.create(rootDir, showWarnings = FALSE)
##' datatable<-geoDownload(datasets = c("GSE1478"), targetDir=rootDir, filePattern="cel$")
geoDownload<-function(datasets, 
                      targetDir = getwd(), 
                      filePattern=NULL, 
                      tar="internal",
                      overwrite=FALSE){
  targetDir<-gsub("[/\\]$","",targetDir)
  counts<-data.frame(dataset = datasets, count = rep(0, length(datasets)))
  for(dataset in datasets){
    subdir<-file.path(targetDir, dataset)
    message(paste0("Dataset ", dataset, " is downloading to ", subdir, " ..."))

    dir.create(subdir, showWarnings = FALSE)
    
    curl<-paste0("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/"
                 , dataset, "/")

    filenames<-getFtpFilenames(curl)

    for(filename in filenames){
      link<-paste0(curl, filename)
      localfile<-paste0(subdir, "/", filename)

      if(overwrite){
        deleteFileAndMd5(localfile)
      }
      
      if(!file.exists(localfile)){
        message("Downloading file ", link, " ...")
        download.file(link, localfile, method="auto", mode="wb")
        Sys.sleep(2)
      }
      
      if(grepl(".tar$", localfile)){
        expectFiles<-paste0(subdir, "/", untar(localfile, tar=tar, list=TRUE))
        unzippedFiles<-gsub("\\.gz$","",expectFiles)
        if(overwrite){
          deleteFilesAndMd5(expectFiles)
          deleteFilesAndMd5(unzippedFiles)
        }
        
        if(all(file.exists(unzippedFiles))){
          next
        }
        
        message("Untar file ", localfile, " ...")
        untar(localfile, exdir=subdir, tar=tar)
        Sys.sleep(2)
        
        gzfiles<-expectFiles[grepl("\\.gz$", expectFiles)]
        for(file in gzfiles){
          message("Decompressing file ", file, " ...")
          tryCatch(gunzip(file, overwrite = TRUE, remove=TRUE),
                   error=function(w){
                     warning(paste0("de-compress file failed for ", file, " : ", w))
                   })
        }
      }
    }
    
    celfiles<-list.files(subdir, filePattern, ignore.case=TRUE)
    counts$count[counts$dataset==dataset] = length(celfiles)
    message("Dataset ", dataset, " has been downloaded.")
  }
  return(counts)
}

##' buildFileTable
##' 
##' The function build file table in the subdirectories under root directories 
##' user provided. The result table contains two columns, dataset and filename
##' 
##' @param rootDir the root of directories whose sub directories contains file 
##'        waiting for validation. It can be vector of directories, or just 
##'        one directory 
##' @param subDirPattern the pattern of sub directory name. Default is NULL.
##' @param filePattern the pattern of file waiting for validation. For example, 
##'        "cel$" for AffyMetrix CEL file only. Default is NULL. 
##' @param ignoreExtensions the extensions of file that will be ignored. 
##'        Default is c("tar", "md5").
##' @param ignore.case ignore the case difference when list files from sub 
##'        directory using filePattern
##' @return a data frame containing full file name and its corresponding 
##'         dataset, which will be used at validateFile
##' @importFrom tools file_ext
##' @export
##' @examples 
##' rootDir<-paste0(dirname(tempdir()), "/DupChecker")
##' datafile<-buildFileTable(rootDir=rootDir, filePattern="cel$")
##' #or
##' datafile<-buildFileTable(rootDir=c(paste0(rootDir, 
##'                          c("/E-MEXP-3872", "/GSE1478") )), filePattern="cel$")
buildFileTable<-function(rootDir, 
                         subDirPattern = NULL, 
                         filePattern = NULL, 
                         ignoreExtensions = c("tar", "md5"),
                         ignore.case=TRUE){
  rootDir<-gsub("[/\\]$","",rootDir)

  lst<-list()
  
  if(is.character(rootDir)){
    dirs<-c(rootDir)    
  }else{
    dirs<-rootDir
  }

  for(dir in dirs){
    subdirs<-list.dirs(dir, recursive=FALSE)
    if(!is.null(subDirPattern)){
      basedirs<-basename(subdirs)
      acceptSubDirs<-grepl(subDirPattern, basedirs, perl=TRUE)
      subdirs<-subset(subdirs, acceptSubDirs)
    }
    
    subdirs<-c(dir, subdirs)
    for(subdir in subdirs){
      basedir<-basename(subdir)
      files<-list.files(subdir, filePattern, full.names=TRUE, no..=TRUE, 
                        ignore.case=ignore.case)
      files<-files[!file.info(files)$isdir]
      for(file in files){
        ext = file_ext(file)
        if(ext %in% ignoreExtensions){
          next
        }
        
        lst<-lappend(lst, c(basedir, file))
      }
    }
  }
  
  result = data.frame(dataset = character(length(lst)), 
                      file = character(length(lst)), stringsAsFactors=FALSE)
  if(length(lst) > 0){
    for(i in c(1:length(lst))){
      result$dataset[i]<-lst[[i]][1]
      result$file[i]<-lst[[i]][2]
    }
  }
  
  return (result)
}

##' validateFile
##' 
##' The function calculate MD5 fingerprint for each file in table and then 
##' check to see if any two files have same MD5 fingerprint. The files with 
##' same fingerprint will be treated as duplication. The function will return 
##' a table contains all duplicated files and datasets.
##' 
##' @param fileTable a table with column name "dataset" and "file", 
##'        here column "file" should contain full name of file.
##' @param saveMd5File if calculated MD5 fingerprint should be save to 
##'        local file
##' @return a list contains two tables. One is the table contains three 
##'         columns: "dataset", "file" and "md5". Another one is the 
##'         duplication table whose row indicates MD5 fingerprint and 
##'         whose column indicates dataset, table cell indicates the 
##'         corresponding filename.
##' @importFrom tools md5sum
##' @export
##' @examples 
##' rootDir<-paste0(dirname(tempdir()), "/DupChecker")
##' datafile<-buildFileTable(rootDir=rootDir)
##' if(nrow(datafile) > 0){
##'   result<-validateFile(datafile)
##'   if(result$hasdup){
##'     duptable<-result$duptable
##'     write.csv(duptable, file="duptable.csv")
##'   }
##' }
validateFile<-function(fileTable, saveMd5File=TRUE){
  message("calculate and validate files, make sure the directory is readable.")
  filemd5<-apply(fileTable, 1, function(x){
    celfile<-as.character(x["file"])
    if(!file.exists(celfile)){
      stop(paste0("File not exists : ", celfile))
    }
    md5file<-paste0(celfile, ".md5")
    if(file.exists(md5file)){
      md5<-readLines(md5file, 1, warn=FALSE)
    }else{
      md5<-""
    }
    
    if(nchar(md5) == 0){
      message("Calculating md5 for ", celfile, "...")
      md5<-md5sum(celfile)
      if(saveMd5File){
        fileConn<-file(md5file)
        writeLines(c(md5), fileConn)
        close(fileConn)
      }
    }
    
    md5
  })
  
  oldtable<-data.frame(fileTable)
  oldtable$md5<-filemd5
  
  md5table<-table(filemd5)
  dupmd5<-names(md5table[md5table > 1])
  
  if(length(dupmd5) > 0){
    dup<-oldtable[oldtable$md5 %in% dupmd5,]
    warning(paste0(nrow(dup), " entries out of total ", nrow(oldtable), 
                   " entries are duplicated at least once."))
    
    dupdatasets<-unique(as.character(dup$dataset))
    dupdatasets<-dupdatasets[order(dupdatasets)]
    x<-"GSE14333"
    dupdsnames<-unlist(lapply(dupdatasets, function(x){
      totalcount<-table(oldtable$dataset==x)["TRUE"]
      dupcount<-table(dup$dataset==x)["TRUE"]
      paste0(x, "[", dupcount, "/", totalcount, "]")
    }))
    
    restable<-matrix(rep("", length(dupmd5) * length(dupdatasets)), 
                     nrow=length(dupmd5), ncol=length(dupdatasets))
    rownames(restable)<-dupmd5
    colnames(restable)<-dupdsnames
    
    for(i in c(1:length(dupmd5))){
      md5<-dupmd5[i]
      ds<-oldtable[oldtable$md5==md5,]
      for(j in c(1:length(dupdatasets))){
        if(dupdatasets[j] %in% ds$dataset){
          restable[i,j]<-paste0(basename(ds[ds$dataset==dupdatasets[j],]$file),
                                collapse=";")
        }    
      }      
    }
    
    result<-list(filetable=oldtable, duptable=restable, hasdup=TRUE)
  }else{
    result<-list(filetable=oldtable, hasdup=FALSE)
  }
  class(result)<-"FileDup"
  result
}
