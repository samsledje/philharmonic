# Philharmonic Dev
In the future, we hope to create a config.yml file that will be passed a single argument to a shell or python script, and the script would take care of running each of the components for you. For now, this is what we have. 

### Steps for running Philharmonic in it's current state

1. install the dependencies

2. Aquire data (2 ways):

- Either Download coral data from the google drive. Running this command will download the prerequisite DSD files for Pdam and SymbC1: ```gdown https://drive.google.com/drive/folders/1hPrkHyTSUak0jnpkZZ_JloUN5eDLvO9S -O ./data --folder```

- OR Run dscript to generate a DSD file. You will also need a go.obo file, and "GO map" mapping protein names to annotations. You can see examples of these by just running the gdown command above. Running DSD might look something like:
```fastDSD -c --converge -t 0.5 --outfile pdam_dscript_distances 50M_20201106_Pdam_predictions.csv.positive```


3. Run the initial Philharmonic clustering algorithm. 
If you used gdown and want to run the clustering algorithm on the SymbC1 data, run the following command:
```python3 clustering.py --net-name SymbC1 --interaction-file ./data/SymbC1_predictions_positive.tsv --dsd-file ./data/SymbC1_dscript_distances.DSD1.tsv --results-dir results/ ```.
As a note, this may need to be done on a high-performance computer. Depending on the size of the dataset, my laptop was not always able to run the recursive clustering. 

4. Run ReCIPE to add back proteins. Navigate to the `ReCIPE` folder then run: 
```python recipe.py --network ../data/SymbC1_predictions_positive.tsv  --cluster-filepath ../results/SymbC1clusters.csv --lr .1 -cthresh 0.75 --outfile SymbC1.recipe_clusters.json --max 20 --metric degree```.
As another note, this process could be improved, and the recipe script would ideally not be run in a separate folder.

5. remove proteins added back more than 15 times and also convert the results from JSON format to CSV 
```python results_convert.py -icf ./results/SymbC1clusters.csv -rrf ./results/SymbC1.recipe_clusters.json --sep , --outfile ./results/SymbC1.recipe_clusters.csv```
Again, in an ideal world, this code would be integrated into the recipe script, so the additional results_convert script wouldn't need to be run.

6. Run the remaining cells in the file `Philharmonic.ipynb` to generate output files and perform analysis on the clusters.
