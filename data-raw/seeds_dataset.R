## code to prepare `seeds_dataset` dataset goes here

seeds <- readr::read_table('seeds_dataset.txt',
                           col_names = c('area',
                                         'perimeter',
                                         'compactness',
                                         'length',
                                         'width',
                                         'asymmetry',
                                         'groove',
                                         'variety' 
                           )
)
usethis::use_data(seeds_dataset, overwrite = TRUE)
