# Silnice v Číslech - Data
Data processing for silnicevcislech.cz

Processing is implemented using Snakemake. 

Hatch is set up for operation (Hatch needs to be installed first):
* To run the workflow: `hatch run snakemake`
* To run the test suite: `hatch run tests:all`
* To check conventions `hatch run conventions:check`
* To format the code using `black`: `hatch run conventions:format`
* To sort the imports using `isort`: `hatch run conventions:imports`

## Crash data from Policie ČR
The workflow requires crash data from Policie ČR as input. Download them here: https://www.policie.cz/clanek/statistika-nehodovosti-900835.aspx

Place the archives for each year in their directories such as `workflow/data/2016/archive.rar`. Currently the archives need to be renamed to `archive.rar`. Processing supports data format from 2016 onward and files for each year have to be present. 

After the workflow succeeds, the results are placed in `workflow/results.json`.