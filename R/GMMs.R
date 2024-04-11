#' Download GMMs file
#' @description Pulls the Raes lab gut metabolic modules file (v1.0.7). Will not overwrite an existing file.
#' @param filepath The filepath where the GMMs file should be downloaded
#' @return Nothing
#' @export
#'
pull_GMMs_file <- function(filepath, version="v1.07"){
    if(!file.exists(filepath)){
        gmm.url <- paste0("https://raw.githubusercontent.com/raeslab/GMMs/master/GMMs.",
                          version, ".txt")
        download.file(gmm.url,
                      filepath)
    } else {
        warning("File already exists.")
    }
}


#' Adds KEGG line to the GMM matrix
#' @description Parser for each line of the GMM file
#' @param GMM.matrix The matrix with GMM data
#' @param line A string containing the line to be parsed and added to the GMM matrix
#' @param mod.name A string containing the name of the module for which current line is being parsed
#' @param mod.number A string containing the ID of the module for which current line is being parsed
#' @return The updated matrix of GMMs, with the new line parsed and added
#'
add_K_line_to_matrix <- function(GMM.matrix, line, mod.name, mod.number){
    # Splits line on either a comma or tab
    split.line <- stringr::str_split(line, "[,\t]")[[1]]
    new.GMM.matrix <- GMM.matrix

    for(entry in split.line){
        new.GMM.matrix <- rbind(new.GMM.matrix, c(entry, mod.name, mod.number))
    }
    return(new.GMM.matrix)

}

#' Gets the GMM matrix in a tabular format
#' @description Download the GMM data (if needed) and parses this file to return a table with the GMMs and corresponding KOs
#' @param filepath The filepath at which the GMMs file should be downloaded or found.
#' @param version The version of GMMs to be used. As of April '24, the newest version is v1.07. This is used in pulling the data from the Raes lab GMMs github, so old versions that are not on the master branch may cause the filepath to not work.
#' @param cleanup A boolean denoting if the filepath should be removed at the end of this function.
#' @return A tabular version of the GMMs as a dataframe
#' @export
#'
get_GMM_matrix <- function(filepath="GMMs.txt", version="v1.07", cleanup=TRUE){
    pull_GMMs_file(filepath=filepath, version=version)

    GMM.fileconts <- readLines(filepath)

    GMM.matrix <- matrix(data=c("KEGG","Module","Module ID"), nrow=1, ncol=3)
    for(i in 1:length(GMM.fileconts)){
        line <- GMM.fileconts[i]
        # if new module, get its name and number
        if(grepl("^MF\\d{4}", line)){
            split.line <- stringr::str_split(line, "\t")[[1]]
            mod.number <- split.line[1]
            mod.name <- split.line[2]
        }
        # If KO line, add that info to the matrix
        else if(grepl("^K\\d{5}", line)){
            GMM.matrix <- add_K_line_to_matrix(GMM.matrix, line, mod.name, mod.number)
        }
    }

    # tidy
    colnames(GMM.matrix) <- GMM.matrix[1,]
    GMM.matrix <- GMM.matrix[2:nrow(GMM.matrix),]
    GMM.matrix <- as.data.frame(GMM.matrix)

    GMM.matrix <- GMM.matrix[!is.na(GMM.matrix$KEGG),]
    GMM.matrix <- GMM.matrix[GMM.matrix$KEGG!="",]

    if(cleanup){
        message(paste("Deleting", filepath))
        file.remove(filepath)
    }
    return(GMM.matrix)
}
