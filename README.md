[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8042339.svg)](https://doi.org/10.5281/zenodo.8042339)

# Stability paper --- Figures

Supporting repository for the stability manuscript. The reporsitory contains code to generate figures and correspdoning statistics where needed. 

## Figures and code
Currently, generated figures are located in the ```output``` folder.
To regenerate the figures, first step would be to downloaded necessary data with the helper script.
### Downloading the data
```get_data.sh``` script is here to access the data used in the work. Example usage: 
 ```bash
 git clone https://github.com/melonheader/stability---figures.git
 cd stability---figures
 bash get_data.sh -o ./
 ```
### Regenerating the figures
Running the entirety of R notebook ```./src/code.RMD``` will save all figures in the ```./output```folder. 
Short descriptions for the outpout could be found above code-chunks in the same notebook.<br>
Opening the repository as an Rproject is recommended.

### Reproducible environment
Currently, the state of the project's libraries is saved with [```renv```](https://rstudio.github.io/renv/articles/renv.html#reproducibility) R package. 
To restore the packages with the same versions used in the initial analysis one can run:
 ```R
 renv::restore()
 ```
in the beginning of the ./src/code.RMD notebook.
