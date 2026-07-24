import pandas as pd

contaminants_zip3 = pd.read_csv(
    "../output/ucmr5_zip3_contaminant.csv", dtype={"Zip3": str}
)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
print(contaminants_zip3.head(n=50))
