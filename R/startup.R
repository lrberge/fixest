#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Apr 23 15:33:05 2022
# ~: Startup related functions
#----------------------------------------------#



#' Permanently removes the fixest package startup message
#'
#' Package startup messages can be very annoying, although sometimes they can be necessary. Use this function to prevent \code{fixest}'s package startup message from popping when loading. This will be specific to your current project.
#'
#' @param x Logical, no default. If \code{FALSE}, the package startup message is removed.
#'
#' @details
#' Note that this function is introduced to cope with the first \code{fixest} startup message (in version 0.9.0).
#'
#' This function works only with R >= 4.0.0. There are no startup messages for R < 4.0.0.
#'
fixest_startup_msg = function(x){

    check_arg(x, "logical scalar mbt")

    config_update("fixest_startup_msg", x)

}

initialize_startup_msg = function(startup_msg){
    # When new versions of the package are installed => we reset the display of the startup message
    # we need to keep track of the versions for which this default has been set

    # NOTA:
    # - the variable fixest_version is written when the user uses fixest_startup_msg()
    # - if this function returns TRUE, then it forces the msg to pop

    # NOTA:
    # - one problem is that I check the version using a local variable
    # specific to a project.
    # - this means that when one creates a new project, the message will necessarily pop!
    # - so I MUST turn off the message for newly created projects.
    # otherwise it would be so annoying.
    # - still => this is a problem if the person uses fixest for the first time
    # -> the project can be deemed old, while in fact fixest was never used
    # so startup messages weren't necessary (bc it would break nothing in the existing code)
    # -> new way: I look at the R files to check whether fixest is used:
    # - if TRUE: startup message
    # - if FALSE: nothing
    # - that's quite costly, but should happen only the very first time the package is attached

    # Note that we must return the value of 'fixest_startup_msg' since these are
    # updated only at session restart (and hence are not directly accessible)

    # message("fixest_startup_msg")

    if(getRversion() < "4.0.0"){
        # No startup message for version < 4.0
        # because there's no way to monitor the messages
        return(FALSE)
    }

    if(is_Rmarkdown()){
        # Never in Rmarkdown: too ugly
        return(FALSE)
    }

    if(is.null(find_project_path())){
        return(FALSE)
    }

    # message("getting version")

    previous_version = config_get("fixest_version")
    is_corrupt_version = !is.null(previous_version) && !is_pkg_version(previous_version)

    # message("version is ", previous_version)

    if(is.null(previous_version)){
        # compatibility with previous versions
        # message("trying to get version from renviron")
        previous_version = renvir_get("fixest_version")
    }

    current_version = fixest_version()

    if(!is_pkg_version(current_version)){
        # If we're here, it's a bug: this should NEVER happen
        return(FALSE)
    }

    if(!is_pkg_version(previous_version)){
        # We first update the version
        # message("updating the version")
        config_update("fixest_version", current_version)

        # message("Is fixest used? ", is_fixest_used())

        # Is it a new project? Or was fixest simply never used before?
        if(!is_corrupt_version && is_fixest_used()){
            # => message
            # Since I register versions since 0.9.0, this means that the
            # version of fixest used was anterior => all msgs should pop

            config_update("fixest_startup_msg", TRUE)
            return(TRUE)
        } else {
            # fixest was never used or the version was corrupt
            # => we don't show any message since it will not break any existing code
            config_update("fixest_startup_msg", FALSE)
            return(FALSE)
        }

        # message("updating done ")

    } else if(!identical(previous_version, current_version)){

        if(version2num(current_version) < version2num(previous_version)){
            # Can happen in projects shared in the cloud
            # In that case, we don't touch the startup message

            msg = paste0("The current project used 'fixest' version ", previous_version, ", but the current version is only ", current_version, ". Maybe update the package?")
            packageStartupMessage(fit_screen(msg, 1))

        } else {

            # A) we update the version
            config_update("fixest_version", current_version)

            # B) we reset the value of fixest_startup_msg
            #    only if the previous_version is anterior to the version that introduced the
            #    message (means the message SHOULD pop since it would be the first time)

            max_version_msg = names(startup_msg)[1]

            if(version2num(previous_version) < version2num(max_version_msg)){
                # You force a startup message even if it was turned off in a previous version

                # use case:
                # - v0.9.0: startup message, user uses fixest_startup_msg(FALSE)
                # - v0.10.0: new breaking changes, you want to inform the user even if he had set
                # fixest_startup_msg(FALSE) in v0.9.0
                #

                config_update("fixest_startup_msg", previous_version)
                return(previous_version)

            } else {
                # The previous version is already posterior to the last message
                # => no startup message any more

                config_update("fixest_startup_msg", FALSE)
                return(FALSE)
            }
        }
    }

    # If null, we'll get the value thanks to renvir_get("fixest_startup_msg")
    # but in some instances, it may be corrupt, so we fix it
    res = config_get("fixest_startup_msg")
    if(is.null(res)){
        # corrupt situation (can occur in dev)
        config_update("fixest_startup_msg", FALSE)
        return(FALSE)
    }

    return(res)
}

version2num = function(x){
    sum(as.numeric(strsplit(x, "\\.")[[1]]) * c(1e6, 1e3, 1))
}

fixest_version = function(){
    as.character(packageVersion("fixest"))
}

is_pkg_version = function(x){
    length(x) == 1 && is.character(x) && length(strsplit(x, "\\.")[[1]]) == 3
}

is_fixest_used = function(){
    # To return TRUE:
    # - fixest in the files
    # - + file saved > 7 days
    #
    # - if fixest but file saved < 7 days, very likely a new project

    # Only level 1 recursivity
    files = list.files(pattern = "\\.(r|R)$")
    dirs = c("./", list.dirs(recursive = FALSE))
    sub_files = unlist(lapply(dirs, list.files, pattern = "\\.(r|R)$", full.names = TRUE))
    file_extra = if(file.exists(".Rprofile")) ".Rprofile" else NULL

    files = c(files, file_extra, sub_files)
    files = files[!dir.exists(files)]

    if(length(files) == 0) return(FALSE)

    big_text = lapply(files, readLines, warn = FALSE)

    # we get the files that have fixest in them
    id_fixest = which(sapply(big_text, function(x) any(grepl("fixest", x, fixed = TRUE))))

    fixest_files = files[id_fixest]
    if(length(fixest_files) == 0) return(FALSE)

    now = Sys.time()

    for(f in fixest_files){
        f_created = file.mtime(f)
        if("POSIXt" %in% class(f_created)){
            d = as.numeric(difftime(now, f_created, units = "days"))
            if(d > 7){
                return(TRUE)
            }
        }
    }

    return(FALSE)
}

renvir_get = function(key){
    # Get the values of envir variables
    # we also evaluate them

    value_raw = Sys.getenv(key)

    if(value_raw == ""){
        return(NULL)
    }

    # Any default value should be able to be evaluated "as such"
    value_clean = gsub("__%%;;", "\n", value_raw)
    value_clean = gsub("&quot;", '"', value_clean)
    value_clean = gsub("&apos;", "'", value_clean)

    value = eval(str2lang(value_clean))

    return(value)
}

find_project_path = function(force = FALSE){
    # finds the root directory
    # we just look up the search path to find the root
    # Only works for Rstudio projects!

    past_path = "init"
    path = normalizePath(".", "/")

    is_found = FALSE
    i = 1
    nmax = 10
    while(past_path != path && i <= nmax){
        i = i + 1
        if(length(list.files(path, pattern = "Rproj$")) > 0){
            is_found = TRUE
            break
        } else {
            past_path = path
            path = dirname(path)
        }
    }

    proj_path = NULL
    if(is_found){
        proj_path = path
    }

    if(force && is.null(proj_path)){
        proj_path = normalizePath(".", "/")
    }

    proj_path
}

find_renviron = function(path = NULL){
    # Simply attaches .Renviron to the project path

    if(is.null(path)){
        proj_path = find_project_path()
        if(is.null(proj_path)) return(NULL)
    } else {
        if(!dir.exists(path)){
            if(file.exists(path)){
                path = dirname(path)
            } else {
                stop_up("The path provided in 'save' does not exist.", .up = 2)
            }
        }

        proj_path = path
    }

    file.path(proj_path, ".Renviron")
}

renvir_update = function(key, value){
    # Updates the .Renviron file
    # asks permission to the user => avoids messing up their workspace!
    # I was thinking to add an argument path, given by the user... but in fact no
    # the .Renviron works only at the Rstudio project level so making the user think
    # that giving a path for saving would help is misleading, since the .Renviron from
    # that path very likely wouldn't be loaded

    check_arg(key, "character scalar mbt")
    check_arg(value, "NULL mbt")

    renv_path = find_renviron()

    if(is.null(renv_path)){
        message("The 'save' feature only works with Rstudio projects. The root directory of the Rstudio project could not be found: settings cannot be saved at the project level, sorry.")
        return(NULL)
    }

    message("To save the settings at the project level 'fixest' needs to update the '.Renviron' file, currently located at:\n\n ", renv_path, "\n\n If the path indeed leads to your current project, do you give persmission? ")

    consent = readline("ok/y/yes to consent:")
    consent = tolower(trimws(consent))

    if(!consent %in% c("ok", "y", "ye", "yes")){
        message("aborting save")
        return(NULL)
    }

    if(file.exists(renv_path)){
        file = file(renv_path, "r", encoding = "UTF-8")

        renvir_raw = readLines(file)

        close(file)
    } else {
        renvir_raw = ""
    }

    all_keys = trimws(gsub("=.*", "", renvir_raw))

    do_write = TRUE
    if(is.null(value)){

        line_to_drop = all_keys == key
        if(any(line_to_drop)){
            renvir_raw = renvir_raw[!line_to_drop]
        } else {
            do_write = TRUE
        }

    } else {

        # we need to do some extra legwork... => sys env don't do quotes
        value_text = paste(deparse(value, width.cutoff = 500), collapse = "\n")
        value_text = gsub("\n", "__%%;;", value_text)
        value_text = gsub("\"", "&quot;", value_text)
        value_text = gsub("'", "&apos;", value_text)

        key_line = all_keys == key
        renvir_raw = c(renvir_raw[!key_line], paste0(key, " = ", value_text))
    }

    if(do_write){
        file = file(renv_path, "w", encoding = "UTF-8")

        renvir_raw = writeLines(renvir_raw, file)

        close(file)
    }


}

find_config_path = function(){

    if(getRversion() < "4.0.0"){
        return(NULL)
    }

    dir = tools::R_user_dir("fixest", "config")

    # We create the directory if needed
    if(!dir.exists(dir)){
        dir.create(dir, recursive = TRUE)
    }

    dir = normalizePath(dir, "/")

    file.path(dir, "fixest_config.csv")
}


config_update = function(key, value){

    if(getRversion() < "4.0.0"){
        return(NULL)
    }

    path = find_config_path()
    proj = find_project_path(force = TRUE)

    if(file.exists(path)){
        data = read.csv(path)
    } else {
        data = data.frame(proj = proj, fixest_version = fixest_version(), stringsAsFactors = FALSE)
    }

    if(!key %in% names(data)){
        data[[key]] = NA_character_
    }

    if(!proj %in% data$proj){
        row = data[1, , drop = FALSE]
        for(i in 1:ncol(row)) row[1, i] = NA
        row[1, 1] = proj
        data = rbind(data, row)
    }

    i = which(data$proj %in% proj)

    data[["fixest_version"]][i] = fixest_version()

    if(is.null(value)) value = "NULL"
    data[[key]][i] = as.character(value)

    write.csv(data, path, row.names = FALSE)
}

config_get = function(key){

    path = find_config_path()

    if(is.null(path) || !file.exists(path)){
        return(NULL)
    }

    data = read.csv(path)

    proj = find_project_path(force = TRUE)

    if(!proj %in% data$proj){
        return(NULL)
    }

    i = which(data$proj %in% proj)

    value = data[[key]][i]

    if(is.character(value) && value %in% c("NULL", "TRUE", "FALSE")){
        value = str2lang(value)
    }

    value
}



