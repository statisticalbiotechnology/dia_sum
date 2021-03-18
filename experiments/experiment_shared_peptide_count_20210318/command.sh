
 curl -s https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz | gunzip | gawk '/^>/{print a;a="";next}{a=a $0}END{print a}' | sed 's/K/K\n/g' | sed 's/R/R\n/g' | gawk 'length($0)>7 {print}' | sort | uniq -c | gawk '{print $1}' | sort -n | uniq -c
