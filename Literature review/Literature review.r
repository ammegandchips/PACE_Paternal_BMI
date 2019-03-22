############################################################################################
#
# Literature search for PACE paternal BMI
#
############################################################################################

############################################################################################
##Load libraries, setwd etc
############################################################################################

setwd("//ads.bris.ac.uk/filestore/MyFiles/Staff17/gs8094/Documents/work/ieu/pace/pat_bmi/literature_search/")

library(RISmed)

############################################################################################
##Define Search terms
############################################################################################

#Paternal related terms
paternal_terms<-c(
"paternal",
"father")
 
#Epigenetics related terms  
epigenetic_terms<-c(
  "methylation")
  
#Obesity/weight related terms
BMI_terms<-c(
  "obesity",
  "overweight",
  "BMI",
  "body mass index",
  "fat mass")

#NOT terms (imprinted syndromes)
not_terms <- c(
"Prader*",
"Angelman*",
"beckwith*")

#get intersects between these terms
searchterms <- as.vector(outer(paternal_terms,epigenetic_terms,paste,sep=" AND "))
searchterms <- as.vector(outer(searchterms,BMI_terms,paste,sep=" AND "))
searchterms <- paste(searchterms,"NOT (",paste(not_terms,collapse=" OR "),")")

print(searchterms)

############################################################################################
## Extract search results from pubmed for each search term and save results object
############################################################################################

search.function <- function(x){
res <- EUtilsSummary(x)
fetch <- EUtilsGet(res)
  out<-as.data.frame(cbind(
    as.character(rep(Query(fetch), times=length(PMID(fetch)))),
    as.character(PMID(fetch)),
    as.character(ELocationID(fetch)),
    as.character(YearPubmed(fetch)),
    as.character(Title(fetch)),
    as.character(lapply(Author(fetch), `[`, 1, 1)),
    as.character(ArticleTitle(fetch)),
    as.character(AbstractText(fetch)),
    as.character(Language(fetch)),
    as.character(PublicationType(fetch))))
  out
}

searchresults <- lapply(searchterms,search.function)
searchresults <- do.call(rbind,searchresults)
searchresults[] <- lapply(searchresults, as.character)
colnames(searchresults)<-c("search query", "PMID", "DOI", "year_pubmed", "journal_name",
                           "first_author", "title", "abstract", "language", "publication_type" )

############################################################################################
## Clean up results object and save clean output
############################################################################################

nrow(searchresults)

#remove any duplicate rows
searchresults<-searchresults[duplicated(searchresults$PMID)==F,]
nrow(searchresults)

#remove any unwanted article types
searchresults$publication_type<-as.character(searchresults$publication_type)

  #keep Journal Articles
  searchresults <- searchresults[grepl("Journal Article", searchresults$publication_type),]
  nrow(searchresults)
  
  #remove reviews
  searchresults <-searchresults[!grepl("Review", searchresults$publication_type),]
  nrow(searchresults)
  
#remove any papers not in English language
searchresults<-searchresults[searchresults$language=="eng",]
nrow(searchresults)

#add extra columns for data extraction
  searchresults <- cbind(searchresults,
	data.frame(Inclusion = rep(NA,nrow(searchresults)),
	reason.for.exclusion = rep(NA,nrow(searchresults)), 
	study.design= rep(NA,nrow(searchresults)), 
	global.ewas.candidate = rep(NA,nrow(searchresults)), 
	exposure = rep(NA,nrow(searchresults)), 
	outcome = rep(NA,nrow(searchresults)), 
	tissue = rep(NA,nrow(searchresults)), 
	n = rep(NA,nrow(searchresults)), 
	n.cases = rep(NA,nrow(searchresults))))
    

#save final object
write.csv(searchresults, paste("./pubmedsearch_", Sys.Date(), "_clean.csv", sep="" ),row.names=FALSE,na="")
