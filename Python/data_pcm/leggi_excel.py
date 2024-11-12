import pandas as pd
import os
import sys


df = pd.read_excel('C:\Volume_D\Programming\Tesi\Tesi_2\Python\data_pcm\IMB_A_2022.xlsx')
print(df)

df.shape
df.drop(df.tail(100).index,inplace=True) # drop last n rows
df.shape
df["PREZZO UNITARIO 2024"] = pd.to_numeric(df['PREZZO UNITARIO 2024'])
df