## code to prepare `DATASET` dataset goes here
load("data-raw/fish.RData")
fulton = data
fulton$day = "Fri" 
fulton[fulton$mon == 1, "day"] = "Mon"
fulton[fulton$tues == 1, "day"] = "Tue"
fulton[fulton$wed == 1, "day"] = "Wed"
fulton[fulton$thur == 1, "day"] = "Thu"
fulton[, c("mon", "tues", "wed", "thur")] <- NULL

fulton <- fulton[, c("t", "day", "avgprc", "totqty", "speed2", "wave2", "speed3", "wave3", "prca", "prcw", "qtya", "qtyw")]
fulton <- setNames(fulton, c("t", "day", "price", "qty", "speed2", "wave2", "speed3", "wave3", "price_asian", "price_white", "qty_asian", "qty_white"))

# head(fulton)
save(fulton, file = "data/fulton.RData")
