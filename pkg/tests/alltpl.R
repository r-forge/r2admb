library(R2admb)
## find ~/admb/admb-trunk/examples -name "*.tpl" -exec cp {} inputs \;
res <- list()
tplfiles <- list.files("inputs",pattern=".tpl")
tpllist <- gsub("\\.tpl$","",tplfiles)
skipfiles <- "glmmadmb"
tpllist <- setdiff(tpllist,skipfiles)
for (i in seq_along(tpllist)) {
    cat(tpllist[[i]],"\n")
    res[[i]] <- R2admb:::read_tpl(file.path("inputs",tpllist[[i]]))
}
# glmmadmb is the only one that fails
## don't know if these are all *right*, but they don't crash ...
