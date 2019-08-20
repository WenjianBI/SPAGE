
# Package Update History

## SPAGE_0.1.4: 
* Add an option 'G.Model' to specify genetic model of "Add" (for additive), "Dom" (for dominant), or "Rec" (for recessive). 
* Subject IDs are required to check if genotype, environmental exposure and phenotypes share the same order. (Suggested by Alexander Rix) 

## SPAGE_0.1.5:
* Update the package in case of errors from SPAtest:::Saddle_Prob(). The errors seem only happen for some rare variants. We ouput NA instead of pval-SPA. (Thank Lars G. Fritsche for pointing this out)
