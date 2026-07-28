
contFactor = 50

path = "/home/projects/ulitsky/shared/data/human/TRIP/LP_run1//"
fileListFile = "LP.txt"

myArgs = commandArgs(trailingOnly = T)
if (length(myArgs)>0)
{
  path=myArgs[1]
  fileListFile=myArgs[2]
}

tabs = list()
beds = list()
fileList = read.table(paste0(path,fileListFile))
for (name in fileList$V1)
{
  print(name)
  tabs[[name]]=read.table(paste0(path,name,".filtered_barcodes_poses.txt"),sep="\t")
  beds[[name]]=read.table(paste0(path,name,".filtered_barcodes.bed"),sep="\t")
  if (ncol(tabs[[name]]) != 11) stop(paste(name, "has", ncol(tabs[[name]]), "columns, expected 11"))
}

allBarcodes = list()
for (group in unique(fileList$V2))
{
  # Find the maximum count for each barcode
  allBarcodes[[group]] = c()
  for (name in fileList[fileList$V2==group,"V1"])
  {
    print(name)
    agg=aggregate(tabs[[name]][,10],by=tabs[[name]][,c(1,2,9)],FUN=max)
    allBarcodes[[group]][paste0(agg$V1,":",agg$V2,":",agg$V9)]=
      apply(data.matrix(cbind(agg$x,allBarcodes[[group]][paste0(agg$V1,":",agg$V2,":",agg$V9)])),1,max,na.rm=T)
  }
  print(paste(length(allBarcodes[[group]])," barcodes recoded"))
}
# Remove contaminants
for (name in fileList$V1)
{
  print(name)
  myGroup = fileList[fileList$V1==name,"V2"]
  allOtherBarcodes = c()
  for (otherGroup in setdiff(unique(fileList$V2),myGroup))
  {
    if (is.null(allOtherBarcodes[names(allBarcodes[[otherGroup]])]))
    {
      allOtherBarcodes=allBarcodes[[otherGroup]]
    } else {
      allOtherBarcodes[names(allBarcodes[[otherGroup]])]=
        apply(data.matrix(cbind(allBarcodes[[otherGroup]],allOtherBarcodes[names(allBarcodes[[otherGroup]])])),1,max,na.rm=T)
    }
  }
  if (length(allOtherBarcodes)>0)
  {
    head(allOtherBarcodes)
    myFilter = tabs[[name]][,10]>=(allOtherBarcodes[match(paste0(tabs[[name]]$V1,":",tabs[[name]]$V2,":",tabs[[name]]$V9),names(allOtherBarcodes))]/contFactor)
  } else {
    myFilter=rep(NA,nrow(tabs[[name]]))
  }
  print(paste(format(length(myFilter),big.mark=","),"barcodes considered"))
  print(paste(format(sum(is.na(myFilter)),big.mark=",")," not observed outside group"))
  print(paste(format(sum(myFilter,na.rm=T),big.mark=","),"to be retained"))
  print(paste(format(sum(!myFilter,na.rm=T),big.mark=","),"to be removed as contaminants"))
  head(myFilter[which(!myFilter)])
  write.table(subset(tabs[[name]],is.na(myFilter) | myFilter),paste0(path,name,".retained_barcodes_counts.txt"),col.names=F,row.names=F,quote=F,sep="\t")
  write.table(subset(beds[[name]],is.na(myFilter) | myFilter),paste0(path,name,".retained_barcodes.bed"),col.names=F,row.names=F,quote=F,sep="\t")
}

