require(phylostratr)
require(magrittr)
require(dplyr)
require(readr)
require(knitr)
require(taxizedb)
require(ape)

focal_taxid <- '4932'

weights <-c(
  '1355160'= 0,
  '245562'= 0,
  '1094981'= 0,
  '284813'= 1.1,
  '766039'= 1.1,
  '237561'= 1.1
)

strata <-
  uniprot_strata(focal_taxid) %>%
  uniprot_strata %>%
  strata_apply(f=diverse_subtree, n=5, weights=uniprot_weight_by_ref()) %>%
  use_recommended_prokaryotes %>%
  uniprot_fill_strata()

saccharomyces <- Strata(
  tree = ape::read.tree(system.file('extdata', 'yeast', 'tree', package='phylostratr')),
  data = list(faa=list(
    # The user is responsible for collecting these data
    Saccharomyces_cerevisiae = 'yeast/cerevisiae.faa',
    Saccharomyces_paradoxus = 'yeast/paradoxus.faa',
    Saccharomyces_mikatae = 'yeast/mikatae.faa',
    Saccharomyces_kudriavzevii = 'yeast/kudriavzevii.faa',
    Saccharomyces_arboricola = 'yeast/arboricola.faa',
    Saccharomyces_eubayanus = 'yeast/eubayanus.faa',
    Saccharomyces_uvarum = 'yeast/uvarum.faa'
  )),
  focal_species = 'Saccharomyces_cerevisiae'
) %>% strata_convert(target='tip', to='id')
strata <- replace_branch(strata, y=saccharomyces, node='4930')

strata <- strata_blast(strata, blast_args=list(nthreads=2)) %>% strata_besthits
results <- merge_besthits(strata)
phlyage <- stratify(result)
hits <- do.call(rbind,strata2@data$besthit)
plot_heatmaps(hits[1:100,],"heatmap.pdf",tree=strata2@tree,focal_id = '3847')
