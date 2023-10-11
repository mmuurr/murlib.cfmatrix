devtools::document()

desc::desc_del("Imports")
list(
    "Imports" = c(
        renv::dependencies("R/")$Package
    ),
    "Depends" = c(),
    "Suggests" = c(),
    "Enhances" = c(),
    "LinkingTo" = c()
) |>
    purrr::iwalk(function(packages, type) {
        purrr::walk(packages, \(package) {
            usethis::use_package(package, type = type, min_version = TRUE)
        })
    })

usethis::use_testthat()
usethis::use_package_doc()
