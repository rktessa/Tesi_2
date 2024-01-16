# Behaviour indicator
'''Examples of plotting data collecting behaviour for some different countries. Used indicators:
- wash_hands_24h_7orMore, 2020-06-27 to 2021-06-18
- contact, 2020-04-23 to 2021-06-18'''
import pandas as pd
import numpy as np
import datetime as DT
import requests
import json 

from matplotlib import pyplot as plt

#matplotlib inline  
#widget
# Contact in the last 24 hours
response_contact_IT = requests.get("https://covidmap.umd.edu/api/resources?indicator=contact&type=smoothed&country=Italy&daterange=20200423-20210618").text
response_contact_FR = requests.get("https://covidmap.umd.edu/api/resources?indicator=contact&type=smoothed&country=France&daterange=20200423-20210618").text
# Wash a lot the hands
response_hands = requests.get("https://covidmap.umd.edu/api/resources?indicator=contact&type=smoothed&country=Italy&daterange=20200627-20210618").text


jsonData_IT = json.loads(response_contact_IT)
jsonData_FR = json.loads(response_contact_FR)

jsonData1 = json.loads(response_hands)

df_contact_IT = pd.DataFrame.from_dict(jsonData_IT['data'])
df_contact_FR = pd.DataFrame.from_dict(jsonData_FR['data'])

df_hands_IT = pd.DataFrame.from_dict(jsonData1['data'])

#df_contact_IT
df_hands_IT
df_contact_IT["survey_date"] = pd.to_datetime(df_contact_IT['survey_date'])
df_contact_FR["survey_date"] = pd.to_datetime(df_contact_FR['survey_date'])

for i in range(3):
    ax = df_contact_IT.plot(x="survey_date", y="smoothed_dc",label="Italy")
    df_contact_FR.plot(x="survey_date", y="smoothed_dc", ax=ax,label="France")
    ax.set_ylabel("% direct contact (> 1 minute) with others ")



df_hands_IT["survey_date"] = pd.to_datetime(df_hands_IT['survey_date'])

ax = df_hands_IT.plot(x="survey_date", y="smoothed_dc",label="Italy")
#df_contact_FR.plot(x="survey_date", y="smoothed_dc", ax=ax,label="France")
ax.set_ylabel("% washed their hands 7+ times in the last 24 hours")

plt.show()