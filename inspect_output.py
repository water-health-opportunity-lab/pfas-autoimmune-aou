import pandas as pd
contaminants_zip3 = pd.read_csv("../output/ucmr5_analysis_zip3.csv", dtype={"Zip3": str})
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
print(contaminants_zip3.head(n=150))
