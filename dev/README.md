# 

to run Philharmonic in dev stage (before creation of config.yml file)

1. install the dependencies:

2. Aquire data (2 ways):

- Download the data from the google drive. Running this command will download the DSD files for Pdam and SymbC1: ```gdown https://drive.google.com/drive/folders/1hPrkHyTSUak0jnpkZZ_JloUN5eDLvO9S -O ./data --folder```

- Run dscript to generate a DSD file:
```fastDSD -c --converge -t 0.5 --outfile pdam_dscript_distances 50M_20201106_Pdam_predictions.csv.positive```

4. Run the initial Philharmonic clustering algorithm. 
If you used gdown and want to run the clustering algorithm on the SymbC1 data, run the following command:
```python3 clustering.py --net-name SymbC1 --interaction-file ./data/SymbC1_predictions_positive.tsv --dsd-file ./data/SymbC1_dscript_distances.DSD1.tsv --results-dir results/ ```

