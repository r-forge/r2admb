library(R2admb)
## find ~/admb/admb-trunk/examples -name "*.tpl" -exec cp {} inputs \;
system('find ~/admb/admb-trunk/examples -name "*.par" -exec cp {} inputs \\;')
res <- list()
## TPL files
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

## PARS
parfiles <- list.files("inputs",pattern=".par")
parlist <- gsub("\\.par$","",parfiles)
## skipfiles <- "glmmadmb"
parlist <- setdiff(parlist,skipfiles)
for (i in seq_along(parlist)) {
    cat(parlist[[i]],"\n")
    res[[i]] <- R2admb:::read_pars(file.path("inputs",parlist[[i]]))
}
