#to do for GC
# distribution dump from Q for the quarter
#all other input files with explanation (cleaned lay summary, and pop list?)
#output report files
#include script as reference if needed

#read in data files
lay <- read.csv('../data/NHGRI_lay_summary_20210119.csv') #cleaned version of the lay summary text for each order (this has been manually curated by team for readability and error correction)
detail <- read.csv('../data/NHGRI_Shipping_History_Detail_110120_013121_BK.csv') #QueueLite distribution data for NHGRI for the quarter (11/01/20-01/31/21), retreived from Q-Depot/Distribution/Shipping History, filtered search by ship date (11/01/20-01/31/21) and by collection type (HG | NHGRI Sample Repository), BK version
#Q0422135942.xls is raw QueueLite dump for same time period and collection
pops <- read.csv('../data/pop_list_20210219.csv') #list of relevant population sample designations

#loop through lay, column 4 to pull out order number, formatted in excel
orderIdLay <- NULL
for(i in 1:nrow(lay)){
	proj <- lay[i,4]
	tmp <- strsplit(as.character(proj), "[(]")[[1]][2]
	tmpb <- strsplit(as.character(tmp), "[)]")[[1]][1]
	orderIdLay <- c(orderIdLay, as.character(tmpb))	
}

#map cleaned lay summaries to detail$Order ID
#make individual vectors of relevant variables(e.g., institution and country)
#for each product (CC or DNA panel) count the #ordered (detail$Indv..Aliquot.Qty.)
#cycle through each population sample
ordersDet <- NULL
ordersLay <- NULL
investigatorLast <-NULL
investigatorFirst <- NULL
names <- NULL
institutions <- NULL
countries <- NULL
prodTypes <- NULL
ordered <- NULL
rintents <- NULL
lays <- NULL
pops <- NULL
prodNum <- NULL

for(i in 1:length(orderIdLay)){
	order <- orderIdLay[i]
	ordDet <- detail[grep(order, detail$Order.Id), ]
	prodtmp <- unique(ordDet$Product)
	poptmp <- unique(ordDet$Diag.Desc)
	nametmp <- unique(as.character(ordDet$Name))
	if(length(nametmp)>1){nametmp <- nametmp[-(which(nametmp==""))]}
	instmp <- unique(as.character(ordDet$Institution))
	counttmp <- unique(as.character(ordDet$Country))
	ordertmp <- unique(as.character(ordDet$Order.Id))
	rintenttmp <- unique(as.character(ordDet$RIntent))
	for(j in 1:length(poptmp)){
			popDet <- ordDet[which(ordDet$Diag.Desc==poptmp[j]),]
			for(k in 1:length(prodtmp)){
			prodDet <- popDet[which(popDet$Product==prodtmp[k]),]
			ordersDet <- c(ordersDet, ordertmp)
			ordersLay <- c(ordersLay, as.character(orderIdLay[i]))
			names <- c(names, nametmp)
			institutions <- c(institutions, instmp)
			countries <- c(countries, counttmp)
			prodTypes <- c(prodTypes, as.character(prodtmp[k]))
			ordered <- c(ordered, as.character(sum(prodDet$Quantity.Ordered)))
			rintents <- c(rintents, rintenttmp)
			lays <- c(lays, as.character(lay[i,1]))
			pops <- c(pops, as.character(poptmp[j]))

		}	
	}
}


#remove new lines from research intents
rinttmp <- rintents
rintentsFor <- gsub("[\r\n]", "", rinttmp)

#create summary table
tableFor <- cbind(ordersDet, ordersLay, pops, names, institutions, countries, prodTypes, ordered, rintentsFor, lays)
colnames(tableFor) <- c("order det", "order lay", "population", "Investigator Name", "institution", "Country", "Product", "# Ordered", "Research Intent", "Lay Summary")


#since looping through each pop and each prod, some have quantity of 0 - need to filter these out
remove <- which(tableFor[,8]=="0")


write.table(tableFor[-(remove),], "combined_comm_report_info_20210219.tsv", quote=FALSE,  sep='\t', row.names=F)
