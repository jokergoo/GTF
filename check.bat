
R CMD REMOVE GTF
Rscript -e "roxygen2::roxygenise()"

cd ..

R CMD build GTF
R CMD check --as-cran --timings GTF_0.0.1.tar.gz
R CMD INSTALL GTF_0.0.1.tar.gz