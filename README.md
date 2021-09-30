# mdtry

Source file .rmd in is a chapter of a bookdown doc. Note that there is no independent .yaml

Chapter was compiled using rmarkdown::render(output_format="github_document"). there's lots of associated options 
(https://rmarkdown.rstudio.com/github_document_format.html), but I used all defaults

The file displays decently here on github, but of course the figure and reference tags for bookdown are bulloxed,
just the same as if we did a non bookdown compile with knitter

I then uploaded compiled file and associated figure directory here. Did a global search and replace on the 
figure folder to list the correct location (https://raw.githubusercontent.com/kcudding/mdtry/main/sample_files/figure-gfm/).

And then just copied to EAforum. 

Essentially did nothing at all except figure folder location correction.
