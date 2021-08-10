
feglm_cases()[19:20,]
K = 18
fml_fixest = feglm_cases()[K,]$fml_fixest
fml_stats = feglm_cases()[K,]$fml_stats
my_family = feglm_cases()[K,]$my_family
my_offset = feglm_cases()[K,]$my_offset
my_weight = feglm_cases()[K,]$my_weight

rm(fml_fixest, fml_stats, my_family, my_offset, my_weight)


