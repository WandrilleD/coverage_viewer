## utilitary functions for the coverage viewer
### ---
### author: "Wandrille Duchemin"
### date: 2024-02-15
### ---
##




## a function to add transcript annotation to the plot
add_transcript_annotations = function( fig , TX_exons , TX_links , TX_names){
  t <- list(
    family = "sans serif",
    size = 12,
    color = toRGB("black"))
  if(nrow(TX_exons) > 0){
    for( ex in rownames( TX_exons ) ) {
      fig = fig %>% add_segments(x = TX_exons[ex,'start'], 
                                 xend = TX_exons[ex,'end'], 
                                 y = TX_exons[ex,'y']*yscale, 
                                 yend = TX_exons[ex,'y']*yscale,
                                 size = I(10), color=I("grey50"),
                                 visible=T)
    }
    for( ex in rownames( TX_links ) ) {
      fig = fig %>% add_segments(x = TX_links[ex,'start'], 
                                 xend = TX_links[ex,'end'], 
                                 y = TX_links[ex,'y']*yscale, 
                                 yend = TX_links[ex,'y']*yscale,
                                 size = I(2), color=I("grey50"),
                                 visible=T)
    }
    for(ex in TX_names){
      fig = fig %>% add_text(x = TX_names[ex,'start']-25, 
                             y = TX_names[ex,'y']*yscale, 
                             text=TX_names[ex,'symbol'],
                             textposition = "middle left",
                             textfont=t, visible=T)
    }
    
  }
  
  fig <- fig %>% rangeslider(autorange =FALSE , range = c(window_start , window_end)) 
  return(fig)  
}

## a helper function we can use to smoothen coverage profiles
rolling_average_subsample = function( x , n  ){
  if( n == 1){return(list(indices = 1:length(x),
                          rolling_average = x))}
  
  cx <- c(0,cumsum(as.numeric(x)))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  rsum_OS = c(rep(NA , floor( n/2 )-1)  , rsum)
  S = seq(floor( n/2 ) , length(rsum) , n)
  return( list( indices = S , rolling_average = rsum_OS[S] ) )
}



## takes:
##   chr start stop
##   refseq_tx < AH object
##   orgdb < AH object
## returns: 
##  list
##    exons ->  data.frame of exons start, end, ID and y-position
##    links ->  data.frame of exons link start, end, ID and y-position
##    names -> data.frame of transcript name ID and position
load_transcript_tracks = function(  chromosome, start , stop , refseq_tx , orgdb , key="REFSEQ"   ){
  
  ss_track = subsetByOverlaps(refseq_tx , GenomicRanges::GRanges(chromosome, 
                                                                 IRanges::IRanges(c(start),c(stop))))
  TX_exons = data.frame(start= c() , end = c() , tx = c() , y = c()  )
  TX_links = data.frame(start= c() , end = c() , tx = c() , y = c()  )
  TX_names = data.frame(start= c() , tx = c() , y = c()  )
  if( length(ss_track$name) > 0 ){
    for( i in 1:length(ss_track$name)){
      
      tx = ss_track[i]
      
      tx_start = tx@ranges@start
      tx_exon = as.data.frame( tx$blocks )
      tx_exon$start = tx_exon$start + tx_start
      tx_exon$end   = tx_exon$end + tx_start
      tx_exon$tx = tx$name
      tx_exon$y = i
      
      tx_link = data.frame( start = tx_exon$end[-nrow(tx_exon)], 
                            end=tx_exon$start[-1]  )
      if(nrow(tx_link)>0){
        tx_link$tx = tx$name
        tx_link$y = i
      }
      
      TX_exons = rbind( TX_exons , tx_exon )
      TX_links = rbind( TX_links , tx_link )
      
      TX_names = rbind(TX_names , 
                       data.frame( start = tx_start , tx= tx$name , y = i ) )
    }
    
    ## refseq id to gene symbol
    SYMBOL = AnnotationDbi::mapIds(orgdb, TX_names$tx, "SYMBOL", key , multiVals='first')
    
    TX_names$symbol = SYMBOL[ TX_names$tx ]
    TX_names$symbol[ is.na(TX_names$symbol) ] = TX_names$tx[is.na(TX_names$symbol)] 
    
  }
  return( list( exons =  TX_exons , links = TX_links , names = TX_names ) )
}


## takes:
##    named vector of coverage files
##    chromosome, start, stop
##    desired number of points / subsample factor
## returns : coverageTracks (rows are pos, columns are coverage tracks)
loadCoverageTracks = function( chromosome, start, stop , coverageFiles , subsample_factor = 10  ){
  BWSp = BigWigSelection(ranges = GenomicRanges::GRanges( chromosome, IRanges::IRanges(c(start),c(stop))), 
                         colnames = "score")
  
  tmp = rolling_average_subsample(window_start:window_end , subsample_factor)
  
  coverageTracks = data.frame(pos=window_start + tmp$indices)
  
  for( sampleID in names(coverageFiles)){
    track = import(coverageFiles[sampleID],
                   format='bw',
                   as = 'NumericList',
                   selection= BWSp )
    
    tmp = rolling_average_subsample(track[[1]] , subsample_factor)
    
    df_tmp = data.frame( tmp$rolling_average )
    names(df_tmp) = sampleID
    coverageTracks = cbind( coverageTracks, df_tmp )
    
  }
  
  print(paste("number of sampled points:",length(unique(coverageTracks$pos))))  
  return(coverageTracks)
}


## takes: 
##  sample table
##  color_column
##  dash_column
##
## return : color_dash_table
make_color_dashes_table = function( sample_table, color_column , dashes_column ){
  
  color_dashes_table = data.frame( row.names = sample_table$sampleID )
  color_dashes_table$color = 'grey'
  color_dashes_table$dash = 'solid'
  
  if( !is.null(color_column) ){
    C = unique( sample_table[,color_column] )
    if(length(C)>12){
      print(paste("ERROR: more than 12 categories in column ",color_column,", which used for color specification:",))
      print(C)
      print("You will have to design your color palette on your own :-(")
    }
    palette = brewer.pal( length(C) , "Set3")  
    names(palette) = C
    color_dashes_table$color = palette[ sample_table[,color_column] ]
  }
  if( !is.null(dashes_column) ){
    C = unique( sample_table[,dashes_column] )
    dashes = c('solid','dash','dashdot','dot')
    if(length(C)>length( dashes )){
      print(paste("ERROR: more than 4 categories in column ",color_column,", which used for dashes specification:",))
      print(C)
      print("You will have to design your color palette on your own :-(")
    }
    dashes = dashes[1:length(C)]
    names(dashes) = C
    color_dashes_table$dash = dashes[ sample_table[,dashes_column] ]
  }
  
  return(color_dashes_table)
}


get_breaks = function(a,e,brk){ (as.integer( a/brk ):as.integer( e/brk ))*brk}


## adding buttons to make certain group of tracks appear/dissapear
## takes:
##     trace_index_list : named list of track indexes. 
##     Ntracks : total number of tracks
##
## returns:
##    list with plotly buttons info. 
##          contains 1 button for each element in the input list, set with the element name as label, and whose action toggle the visibility of the associated track indexes
##           + 1 button "ALL" which toggles all track visibility
add_toggle_buttons = function( trace_index_list , Ntracks  ){
  
  
  updatemenus3 <- list(
    list(
      active = -1,
      type = 'buttons',
      buttons = list(
        list(
          label = "ALL",
          method = "restyle",
          args = list("visible", as.list( rep(TRUE,Ntracks)  ) , 1:Ntracks ),
          args = list("visible", as.list( rep(FALSE,Ntracks) ) , 1:Ntracks )
        )
      )
    )
  )
  
  n=2
  for( column in categoryButton_columns ){
    for( item in unique( sample_table[,column] )  ){
      updatemenus3[[1]]$buttons[[n]] = list( label = paste(item,"toggle"), method = "restyle",
                                             args  = list("visible", 
                                                          as.list( rep(FALSE , length(trace_index_list[[item]]) ) ) , 
                                                          trace_index_list[[item]] )  ,
                                             args2 = list("visible", 
                                                          as.list( rep(TRUE , length(trace_index_list[[item]]) ) ) , 
                                                          trace_index_list[[item]] ) ) 
      n = n+1    
    }
  }
  return(updatemenus3)
}
