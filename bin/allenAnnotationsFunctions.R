#------------------ Librairies ------------------  
library(rjson)



#------------------ Ontology handling ------------------ 

get.id.from.acronym <- function(acronym){
  
  rec.func <- function(acronym, json){
    
    if(is.null(json))
      return()
    
    if(acronym == json$acronym)
      return(json$id)
    
    return(unlist(lapply(json$children, function(x){return(rec.func(acronym, x))})))
  }
  
  #Loading json data
  json_data <- fromJSON(file= paste(path.matrices,'ontology.json',sep='/'))
  json <- json_data$msg[[1]]
  
  return(rec.func(acronym, json))
}



#------------------ 3D handling ------------------ 

mesh3d.allen.annot.from.id <- function(id, max.nrows = 349696, bregma = c(5200, 650, 5700)/1000, no.normals = F){ #5400, 650, 5700 In mm, AP/DV/ML
  
  fpath <- sprintf('%s/%d.obj', allen.annot.path, id)
  
  vert.start <- 0
  vert.end <- 0
  norm.start <- 0
  norm.end <- 0
  face.start <- 0
  face.end <- 0
  i <- 0
  
  con <- file(fpath, 'r')
  while(TRUE){
    
    i <- i + 1
    line = readLines(con, n = 1)
    
    if (length(line) == 0)
      break
    
    if(startsWith(line, 'v ') & vert.start == 0)
      vert.start <- i
    
    if(vert.start != 0 & vert.end == 0 & !startsWith(line, 'v '))
      vert.end <- i
    
    if(startsWith(line, 'vn') & norm.start == 0)
      norm.start <- i
    
    if(norm.start !=0 & norm.end == 0 & !startsWith(line, 'vn'))
      norm.end <- i
    
    if(startsWith(line, 'f') & face.start == 0)
      face.start <- i
    
    if(face.start !=0 & face.end == 0 & !startsWith(line, 'f'))
      face.end <- i
    
  }
  
  if(face.end == 0)
    face.end <- i - 1
  
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (vert.start - 1))
  line = readLines(con, n = vert.end - vert.start)
  vertices <- strsplit(line, ' ')
  vertices.ul <- unlist(vertices)
  df.vertices <- data.frame(x = as.numeric(vertices.ul[seq(from = 2, to = 4*length(vertices), by = 4)]),
                            y = as.numeric(vertices.ul[seq(from = 3, to = 4*length(vertices), by = 4)]),
                            z = as.numeric(vertices.ul[seq(from = 4, to = 4*length(vertices), by = 4)]))
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (norm.start - 1))
  line = readLines(con, n = norm.end - norm.start)
  norm <- strsplit(line, ' ')
  norm.ul <- unlist(vertices)
  df.norm <-  data.frame(   x = as.numeric(norm.ul[seq(from = 2, to = 4*length(norm), by = 4)]),
                            y = as.numeric(norm.ul[seq(from = 3, to = 4*length(norm), by = 4)]),
                            z = as.numeric(norm.ul[seq(from = 4, to = 4*length(norm), by = 4)]))
  close(con)
  
  con <- file(fpath, 'r')
  line = readLines(con, n = (face.start - 1))
  line = readLines(con, n = face.end - face.start)
  face <- strsplit(line, ' ')
  face <- lapply(face, function(x){return(unlist(strsplit(x, '//')))})
  face.ul <- unlist(face)
  
  df.face <-  data.frame(   x = as.numeric(face.ul[seq(from = 2, to = 7*length(face), by = 7)]),
                            y = as.numeric(face.ul[seq(from = 4, to = 7*length(face), by = 7)]),
                            z = as.numeric(face.ul[seq(from = 6, to = 7*length(face), by = 7)]))
  close(con)
  
  if(!no.normals){
    mesh <- tmesh3d(vertices = rbind(t(df.vertices),1),
                    indices = t(df.face),
                    normals = t(df.norm)[c(3,2,1),])
  }else{
    mesh <- tmesh3d(vertices = rbind(t(df.vertices),1),
                    indices = t(df.face))
  }
  
  mesh$vb[1:3,] <- mesh$vb[1:3,] / 1000 
  
  mesh$vb[1,] <- -mesh$vb[1,] + bregma[1]
  mesh$vb[2,] <- -mesh$vb[2,] + 0.05 #Maybe not the +0.05
  mesh$vb[3,] <- -mesh$vb[3,] + bregma[3]
  
  tmp <- mesh$vb[1,]  
  mesh$vb[1,] <- mesh$vb[3,]
  mesh$vb[3,] <- tmp
  
  return(mesh)
  
}



mesh3d.show.outline <- function(col = 'lightgray', alpha = 0.1){
  shade3d(mesh3d.allen.annot.from.id(get.id.from.acronym('root')), col = col, alpha = alpha)
}

mesh3d.new.window <- function(show.outline = TRUE){
  
  rgl.open()
  rgl.bg(color = 'white')
  par3d(windowRect = c(0, 0, 1920, 1080))
  
  if(show.outline)
    mesh3d.show.outline()
}




