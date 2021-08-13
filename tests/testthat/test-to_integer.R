####
#### To Integer ####
####

chunk("TO_INTEGER")

base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$z <- sample(5, 150, TRUE)

# Normal
m <- to_integer(base$species)
test(length(unique(m)), 3)

m <- to_integer(base$species, base$z)
test(length(unique(m)), 15)

# with NA
base$species_na <- base$species
base$species_na[base$species == "setosa"] <- NA

m <- to_integer(base$species_na, base$z)
test(length(unique(m)), 11)

m <- to_integer(base$species_na, base$z, add_items = TRUE, items.list = TRUE)
test(length(m$items), 10)
