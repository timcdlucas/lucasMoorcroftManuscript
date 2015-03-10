# Repo for "A generalised random encounter model for estimating animal density with remote sensor data"

This is a repo with our manuscript and code.

The fully formatted paper can be found, open access, here [http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12346/abstract](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12346/abstract)


## Navigation

### Branches

There are two main branches.
[master](https://github.com/timcdlucas/lucasMoorcroftManuscript/tree/master) is our general use branch and is a mess.
[postPeerReview](https://github.com/timcdlucas/lucasMoorcroftManuscript/tree/postPeerReview) is a tidy branch with just the files needed for the manuscript.


### postPeerReview

[lucas_et_al_mainms_*date*.tex](lucas_et_al_mainms_2015-01-20.tex) contains the main manuscript.

[lucas_et_al_mainms_*date*.pdf](lucas_et_al_mainms_2015-01-20.pdf) is a compiled, unformatted pdf of the manuscript.

[lucas_et_al_tex_suppl_file_2015-01-20.bib](lucas_et_al_tex_suppl_file_2015-01-20.bib) and [lucas_et_al_tex_suppl_file_2015-01-20.bst](lucas_et_al_tex_suppl_file_2015-01-20.bst) are the latex reference and style file.

[imgs/](imgs/) contains all the .pdf figures for the paper.


[supplementary-material/](supplementary-material/) contains:
- [REM-methods.tex](supplementary-material/REM-methods.tex), [lucas_et_al_supplementarymaterial_2015-01-20.tex](supplementary-material/lucas_et_al_supplementarymaterial_2015-01-20.tex) and a directory of .tex files,  [latexFiles/](https://github.com/timcdlucas/lucasMoorcroftManuscript/tree/postPeerReview/supplementary-material/latexFiles), that make up the text of the supplementary material.
- [imgs/](https://github.com/timcdlucas/lucasMoorcroftManuscript/tree/postPeerReview/imgs) contains the pdf images needed for the supplementary material. 

Code:
The code supplement include with the paper includes
- [lucas_et_al_S3.py](supplementary-material/lucas_et_al_S3.py) a python script which performs the analytical maths for the paper. It also provides checks that the model specifications are correct. 
- [lucas_et_al_S4.R](supplementary-material/lucas_et_al_S4.R) is a simple script that implements the gREM in R.
- [supplementary-material/testSupplementaryRscript.R](testSupplementaryRscript.R) contains some unit tests for the functions in [lucas_et_al_S4.R](supplementary-material/lucas_et_al_S4.R).




### master

This branch is basically a mess. It follows the same structure as `postPeerReview` but without any of the tidying. 






 
